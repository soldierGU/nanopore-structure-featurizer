# Nanopore Structure Featurizer

`nanopore-structure-featurizer` 用于从蛋白纳米孔三维结构中提取可分析、可建模的结构特征。

工程当前支持两类输入：

1. WT 结构：`PDB/mmCIF`
2. 基于 FoldX 建模得到的 mutant 结构：`PDB`

当前仓库已经验证过以下场景：

- WT: `5JZT`
- 单突变: `T232K`
- 双突变: `T232K/K238Q`

同样的流程也可以用于新的单点突变，例如 `T274S`。

## 1. 工程做什么

对于一个纳米孔结构，工程会依次生成：

1. chain-level summary
2. residue-level features
3. candidate inner residues
4. nanopore-level structure features
5. WT-mutant delta features

主流程如下：

```text
PDB/mmCIF / modeled PDB
        |
chain summary
        |
residue features
        |
candidate inner residues
        |
nanopore-level features
        |
WT-mutant delta features
```

## 2. 目录结构

```text
nanopore-structure-featurizer/
├─ config/
├─ data/
│  ├─ raw/
│  ├─ foldx/
│  │  ├─ input/
│  │  └─ mutants/
│  ├─ modeled/
│  └─ processed/
├─ outputs/
├─ scripts/
└─ src/npstructfeat/
```

## 3. 依赖与运行环境

安装依赖：

```powershell
pip install -r requirements.txt
```

建议使用 Python 3.10+。

由于工程采用 `src/` 布局，在 Windows PowerShell 下建议先执行：

```powershell
cd E:\nanopore-structure-featurizer
chcp 65001
$env:PYTHONPATH="src"
```

## 4. 当前配置文件

### 4.1 WT 默认配置

`config/default.yaml`

用于 WT `5JZT` biological assembly。

### 4.2 已有 mutant 配置

- `config/5JZT_T232K_foldx.yaml`
- `config/5JZT_T232K_K238Q_foldx.yaml`

这类配置用于 FoldX 建模后的 mutant PDB，再跑后续特征提取流程。

## 5. 基础用法

### 5.1 运行 WT 全流程

```powershell
python scripts/run_pipeline.py --config config/default.yaml
```

如果你想逐步执行：

```powershell
python scripts/00_check_config.py
python scripts/01_parse_structure.py
python scripts/02_extract_residue_features.py
python scripts/03_select_inner_residues.py
python scripts/04_build_nanopore_features.py
```

### 5.2 运行分析和可视化

```powershell
python scripts/05_threshold_sensitivity.py
python scripts/06_visualize_inner_residues.py
```

## 6. 这份工程是否兼容新的突变

兼容。

### 6.1 新的单点突变

例如：

- `T274S`
- `K238Q`
- `R220A`

当前工程已经支持新的单点突变，不需要为单点突变再改 Python 代码。

### 6.2 新的多点突变

例如：

- `T232K/K238Q`
- `T274S/K238Q`

当前工程已经支持把多点突变字符串拆分后传给：

- `07_check_mutation_site.py`
- `09_run_foldx_buildmodel.py`

同时也支持为多突变生成规范化文件名和目录名，例如：

```text
T232K/K238Q -> T232K_K238Q
```

## 7. 新突变的标准执行流程

下面这套流程适用于：

- 新的单点突变，例如 `T274S`
- 新的多点突变，例如 `T232K/K238Q`

### 第 1 步：检查 WT 结构中的突变位点

先确认 WT 结构里目标位点存在，而且残基编号和残基类型都匹配。

#### 单点突变示例：`T274S`

```powershell
python scripts/07_check_mutation_site.py --config config/default.yaml --mutation T274S
```

#### 多点突变示例：`T232K/K238Q`

```powershell
python scripts/07_check_mutation_site.py --config config/default.yaml --mutation T232K/K238Q
```

输出结果保存在：

```text
data/processed/mutation_sites/
```

你应该重点检查这些字段：

- `mutation_id`
- `chain_id`
- `residue_number`
- `expected_residue`
- `observed_residue`
- `is_site_found`
- `is_expected_match`
- `is_inner_candidate`

只有当位点检查通过后，才建议继续做 FoldX 建模。

### 第 2 步：准备 FoldX 输入 PDB

```powershell
python scripts/08_prepare_foldx_input.py --config config/default.yaml --mutation T274S --chains A B C D E F G
```

或者：

```powershell
python scripts/08_prepare_foldx_input.py --config config/default.yaml --mutation T232K/K238Q --chains A B C D E F G
```

这一步会生成：

```text
data/foldx/input/5JZT_WT/5JZT_assembly.pdb
```

以及 FoldX 输入准备报告。

### 第 3 步：运行 FoldX RepairPDB

进入 WT FoldX 输入目录：

```powershell
cd E:\nanopore-structure-featurizer\data\foldx\input\5JZT_WT
```

运行：

```powershell
E:\nanopore-structure-featurizer\tools\foldx\foldx_1_20270131.exe --command=RepairPDB --pdb=5JZT_assembly.pdb
```

会生成：

```text
data/foldx/input/5JZT_WT/5JZT_assembly_Repair.pdb
```

#### 什么时候可以复用已有的 `5JZT_assembly_Repair.pdb`

如果下面这件事没有变，就可以直接复用，不必重复跑：

- `data/foldx/input/5JZT_WT/5JZT_assembly.pdb` 的内容没有被重新生成或手工修改过

换句话说，**只要 WT FoldX 输入 PDB 还是同一个文件内容，已有的 `RepairPDB` 结果就可以继续用。**

### 第 4 步：先准备 BuildModel 工作目录

建议先准备工作目录和 `individual_list.txt`，不要直接运行 FoldX：

#### 单点突变示例：`T274S`

```powershell
cd E:\nanopore-structure-featurizer
$env:PYTHONPATH="src"
python scripts/09_run_foldx_buildmodel.py --config config/default.yaml --mutation T274S --chains A B C D E F G --prepare-only --overwrite-individual-list
```

#### 多点突变示例：`T232K/K238Q`

```powershell
cd E:\nanopore-structure-featurizer
$env:PYTHONPATH="src"
python scripts/09_run_foldx_buildmodel.py --config config/default.yaml --mutation T232K/K238Q --chains A B C D E F G --prepare-only --overwrite-individual-list
```

这一步会生成：

```text
data/foldx/mutants/<PDB_ID>_<MUTATION_LABEL>/individual_list.txt
```

你应该先检查 `individual_list.txt` 内容是否符合预期。

#### 单点示例 `T274S`

```text
TA274S,TB274S,TC274S,TD274S,TE274S,TF274S,TG274S;
```

#### 多点示例 `T232K/K238Q`

```text
TA232K,TB232K,TC232K,TD232K,TE232K,TF232K,TG232K,KA238Q,KB238Q,KC238Q,KD238Q,KE238Q,KF238Q,KG238Q;
```

### 第 5 步：正式运行 FoldX BuildModel

确认 `individual_list.txt` 没问题后执行：

#### 单点突变 `T274S`

```powershell
python scripts/09_run_foldx_buildmodel.py --config config/default.yaml --mutation T274S --chains A B C D E F G
```

#### 多点突变 `T232K/K238Q`

```powershell
python scripts/09_run_foldx_buildmodel.py --config config/default.yaml --mutation T232K/K238Q --chains A B C D E F G
```

这一步会：

1. 调用 FoldX `BuildModel`
2. 收集 mutant PDB
3. 验证每条链上目标位点是否真的变成目标残基

输出主要位于：

```text
data/foldx/mutants/
data/modeled/
```

### 第 6 步：为新的 mutant 新建 YAML 配置

每做出一个新的 mutant，都建议新建一个独立 YAML。

原因：

1. mutant PDB 是新的输入结构
2. 需要单独的 `pdb_id`
3. 需要单独的 `nanopore_id`
4. 需要把 `structure_file` 指向新的 modeled PDB

#### 单点突变 `T274S` 的 YAML 应该怎么写

建议新建：

```text
config/5JZT_T274S_foldx.yaml
```

关键字段示例：

```yaml
input:
  pdb_id: 5JZT_T274S
  nanopore_id: aerolysin_T274S
  structure_file: data/modeled/5JZT_T274S/5JZT_T274S_model.pdb
  file_format: pdb
  assembly_type: modeled_biological_assembly
  structure_source_type: template_based_mutant_model
  template_pdb_id: 5JZT
  modeling_method: FoldX
  mutation_list:
    - T274S
```

#### 多点突变 `T232K/K238Q` 的 YAML

本仓库已提供：

```text
config/5JZT_T232K_K238Q_foldx.yaml
```

### 第 7 步：对 mutant 重新跑结构特征提取流程

#### 单点突变 `T274S`

```powershell
python scripts/run_pipeline.py --config config/5JZT_T274S_foldx.yaml
```

#### 多点突变 `T232K/K238Q`

```powershell
python scripts/run_pipeline.py --config config/5JZT_T232K_K238Q_foldx.yaml
```

这一步会生成 mutant 自己的：

- chain summary
- residue features
- inner candidate residues
- nanopore features

### 第 8 步：比较 WT 和 mutant

#### 单点突变 `T274S`

```powershell
python scripts/10_compare_wt_mutant_features.py --wt-feature data/processed/nanopore_features/5JZT_nanopore_structure_features.csv --mutant-feature data/processed/nanopore_features/5JZT_T274S_nanopore_structure_features.csv --wt-id aerolysin_WT --mutant-id aerolysin_T274S --template-pdb-id 5JZT --mutation T274S --modeling-method FoldX
```

#### 多点突变 `T232K/K238Q`

```powershell
python scripts/10_compare_wt_mutant_features.py --wt-feature data/processed/nanopore_features/5JZT_nanopore_structure_features.csv --mutant-feature data/processed/nanopore_features/5JZT_T232K_K238Q_nanopore_structure_features.csv --wt-id aerolysin_WT --mutant-id aerolysin_T232K_K238Q --template-pdb-id 5JZT --mutation T232K/K238Q --modeling-method FoldX
```

输出在：

```text
data/processed/delta_features/
```

## 8. 处理新的单点突变 `T274S` 的完整命令示例

```powershell
cd E:\nanopore-structure-featurizer
chcp 65001
$env:PYTHONPATH="src"

python scripts/07_check_mutation_site.py --config config/default.yaml --mutation T274S
python scripts/08_prepare_foldx_input.py --config config/default.yaml --mutation T274S --chains A B C D E F G
```

如果 WT 的 `5JZT_assembly_Repair.pdb` 已经存在且对应的 WT 输入 PDB 没变，可以直接跳过 RepairPDB。

然后继续：

```powershell
python scripts/09_run_foldx_buildmodel.py --config config/default.yaml --mutation T274S --chains A B C D E F G --prepare-only --overwrite-individual-list
python scripts/09_run_foldx_buildmodel.py --config config/default.yaml --mutation T274S --chains A B C D E F G
```

新建 `config/5JZT_T274S_foldx.yaml` 后执行：

```powershell
python scripts/run_pipeline.py --config config/5JZT_T274S_foldx.yaml
python scripts/10_compare_wt_mutant_features.py --wt-feature data/processed/nanopore_features/5JZT_nanopore_structure_features.csv --mutant-feature data/processed/nanopore_features/5JZT_T274S_nanopore_structure_features.csv --wt-id aerolysin_WT --mutant-id aerolysin_T274S --template-pdb-id 5JZT --mutation T274S --modeling-method FoldX
```

## 9. 结果文件说明

### WT / mutant 结构特征

```text
data/processed/chain_summary/
data/processed/residue_features/
data/processed/inner_residues/
data/processed/nanopore_features/
```

### 突变位点检查

```text
data/processed/mutation_sites/
```

### FoldX 工作目录与模型

```text
data/foldx/input/
data/foldx/mutants/
data/modeled/
```

### WT-mutant 差异

```text
data/processed/delta_features/
```

## 10. 当前方法的边界

当前工程的孔道残基筛选是几何近似版本，核心规则是：

```text
radial_distance <= inner_radius_threshold
```

因此：

1. 这里的 `candidate inner residues` 不等于严格定义的 pore-lining residues
2. 如果你看到某些突变对全局统计影响很小，不一定说明突变没有作用
3. 多数情况下，还需要进一步做局部位点邻域分析

## 11. 建议的使用习惯

对于每一个新的 mutant，建议固定遵循下面的节奏：

1. 先做位点检查
2. 再准备 FoldX 输入
3. 再做 BuildModel
4. 为 mutant 建一个独立 YAML
5. 再跑后续特征流程
6. 最后做 WT-mutant 对比

这样最容易排查问题，也最不容易把不同 mutant 的结果混在一起。
