## Nanopore Structure Featurizer

`nanopore-structure-featurizer` 是一个面向蛋白纳米孔结构的特征提取工程。
当前版本以 `PDB/mmCIF` 结构文件为输入，围绕孔道几何中心和残基理化属性，生成：

- chain-level 结构摘要
- residue-level 特征表
- 候选孔道内壁残基表
- nanopore-level 聚合结构特征表

工程当前已经包含一个示例输入：`5JZT` biological assembly。

## 当前流程

主入口是 [scripts/run_pipeline.py](scripts/run_pipeline.py)。

完整流程包含 4 个步骤：

1. 解析结构并生成链级统计表
2. 提取标准氨基酸残基的基础特征
3. 根据径向距离阈值筛选候选孔道内壁残基
4. 聚合得到 nanopore-level 结构特征

当前筛选逻辑是一个几何近似版本，核心规则为：

`radial_distance <= inner_radius_threshold`

这适合先完成可运行的数据流程，但它不等同于专业孔道识别工具给出的严格 pore-lining residue 定义。

## 目录结构

```text
nanopore-structure-featurizer/
├─ config/
│  ├─ default.yaml
│  └─ nanopores.yaml
├─ data/
│  ├─ raw/
│  └─ processed/
├─ outputs/
│  └─ figures/
├─ scripts/
│  ├─ 00_check_config.py
│  ├─ 01_parse_structure.py
│  ├─ 02_extract_residue_features.py
│  ├─ 03_select_inner_residues.py
│  ├─ 04_build_nanopore_features.py
│  ├─ 05_threshold_sensitivity.py
│  ├─ 06_visualize_inner_residues.py
│  └─ run_pipeline.py
├─ src/
│  └─ npstructfeat/
│     ├─ io.py
│     ├─ parser.py
│     ├─ features.py
│     ├─ geometry.py
│     ├─ pore.py
│     ├─ residue_props.py
│     └─ utils.py
├─ tests/
├─ requirements.txt
└─ README.md
```

## 环境依赖

依赖见 [requirements.txt](requirements.txt)：

- `biopython`
- `pandas`
- `numpy`
- `pyyaml`
- `matplotlib`

建议使用 Python 3.10+。

安装依赖：

```bash
pip install -r requirements.txt
```

## 输入数据

默认配置文件是 [config/default.yaml](config/default.yaml)。

当前默认输入：

- `pdb_id`: `5JZT`
- `nanopore_id`: `aerolysin_WT`
- `structure_file`: `data/raw/5JZT/5JZT_assembly.cif`
- `file_format`: `cif`
- `assembly_type`: `biological_assembly`

## 配置说明

`config/default.yaml` 主要包含以下配置段：

### `project`

- `name`
- `version`

### `input`

- `pdb_id`
- `nanopore_id`
- `structure_file`
- `file_format`
- `assembly_type`

### `structure`

- `use_model_index`: 选择第几个 model，默认 `0`
- `use_standard_residues_only`: 是否只保留标准氨基酸
- `use_ca_only`: 是否只保留含 `CA` 原子的残基

### `pore`

- `axis_mode`: 当前版本仅支持 `z_axis`
- `center_mode`: 当前支持 `xy_mean` 和 `xy_median`
- `inner_radius_threshold`: 候选孔道内壁筛选阈值，单位 Å

### `output`

定义中间结果和处理结果的输出目录。

### `save`

控制是否保存各阶段结果。

## 运行方式

### 运行完整流程

```bash
python scripts/run_pipeline.py
```

或显式指定配置文件：

```bash
python scripts/run_pipeline.py --config config/default.yaml
```

### 按步骤运行

```bash
python scripts/00_check_config.py
python scripts/01_parse_structure.py
python scripts/02_extract_residue_features.py
python scripts/03_select_inner_residues.py
python scripts/04_build_nanopore_features.py
```

### 分析和可视化脚本

```bash
python scripts/05_threshold_sensitivity.py
python scripts/06_visualize_inner_residues.py
```

## 输出文件

基于当前默认配置，工程会在 `data/processed/` 下生成：

### 1. chain summary

路径示例：

`data/processed/chain_summary/5JZT_chain_summary.csv`

包含：

- `chain_id`
- `total_residue_count`
- `standard_residue_count`
- `nonstandard_residue_count`
- `atom_count`
- `ca_missing_count`
- `first_residue_number`
- `last_residue_number`

### 2. residue features

路径示例：

`data/processed/residue_features/5JZT_residue_features.csv`

包含：

- `chain_id`
- `residue_number`
- `insertion_code`
- `residue_name`
- `x`, `y`, `z`
- `has_ca`
- `charge`
- `hydrophobicity`
- `is_aromatic`
- `is_polar`

### 3. inner candidate residues

路径示例：

`data/processed/inner_residues/5JZT_inner_candidate_residues.csv`

在 residue-level 特征基础上增加：

- `center_x`
- `center_y`
- `radial_distance`
- `theta_rad`
- `theta_deg`
- `z_norm`
- `is_inner_candidate`

### 4. nanopore structure features

路径示例：

`data/processed/nanopore_features/5JZT_nanopore_structure_features.csv`

包含聚合后的孔道级特征，例如：

- 候选内壁残基数量与比例
- 每条链的候选残基数量统计
- `z` 方向覆盖范围
- 径向距离统计
- 净电荷与平均电荷
- 平均疏水性
- 芳香族和极性残基占比
- 候选残基类型数

## 核心模块

- [src/npstructfeat/io.py](src/npstructfeat/io.py): 配置读取、路径解析、输出目录准备
- [src/npstructfeat/parser.py](src/npstructfeat/parser.py): 结构加载、标准残基判断、链级摘要
- [src/npstructfeat/features.py](src/npstructfeat/features.py): 残基级特征提取与 nanopore-level 聚合
- [src/npstructfeat/geometry.py](src/npstructfeat/geometry.py): 几何中心、径向距离、角度、z 归一化
- [src/npstructfeat/pore.py](src/npstructfeat/pore.py): 候选孔道内壁残基筛选
- [src/npstructfeat/residue_props.py](src/npstructfeat/residue_props.py): 氨基酸理化属性定义

## 当前假设与限制

当前版本依赖以下假设：

- 结构已经基本沿 `z` 轴对齐
- 孔道中心可以用所有残基 `CA` 坐标在 `x-y` 平面上的均值或中位数近似
- 候选孔道内壁残基可通过径向距离阈值粗筛
- 默认关注标准氨基酸残基

因此，这一版更适合：

- 建立可复用的数据处理流程
- 快速生成结构特征表
- 做阈值敏感性分析和初步探索

如果后续需要更严格的孔道识别，可以考虑引入更稳健的主轴估计方法或外部专业工具结果。
