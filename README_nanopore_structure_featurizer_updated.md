# Nanopore Structure Featurizer

`nanopore-structure-featurizer` is a structure feature extraction project for protein nanopores. It takes `PDB/mmCIF` structure files as input and gradually converts three-dimensional nanopore structures into chain-level, residue-level, candidate pore-region-level, and nanopore-level feature tables.

The current project uses the `5JZT` wild-type aerolysin biological assembly as the initial example and has been extended to support FoldX-based template mutant modeling, using the `T232K` aerolysin mutant as the first case.

---

## 1. Project Goal

The goal of this project is to transform nanopore protein structures into structured, interpretable, and model-ready feature tables.

The overall workflow is:

```text
Nanopore PDB/mmCIF structure
        ↓
chain-level summary
        ↓
residue-level features
        ↓
candidate pore-region residues
        ↓
nanopore-level structure features
        ↓
WT-mutant delta features
```

These structural features can later be integrated with analyte sequence features and experimental condition features for nanopore sensing response modeling:

```text
Analyte features + Nanopore structural features + Experimental conditions
→ Blockade time / Residual current
```

---

## 2. Current Main Functions

The current project supports:

1. Loading `PDB/mmCIF` structure files;
2. Parsing biological assemblies;
3. Generating chain-level structure summaries;
4. Extracting Cα coordinates and physicochemical properties of standard amino acid residues;
5. Selecting candidate pore-lining residues using a radial-distance criterion;
6. Aggregating candidate pore-region residues into nanopore-level structural features;
7. Performing radial-distance threshold sensitivity analysis;
8. Visualizing candidate inner residue selection results;
9. Checking mutation sites before mutant modeling;
10. Preparing FoldX-compatible input PDB files;
11. Building template-based mutant models using FoldX;
12. Extracting structural features from FoldX-modeled mutants;
13. Comparing WT and mutant structural features.

---

## 3. Repository Structure

```text
nanopore-structure-featurizer/
├── config/
│   ├── default.yaml
│   ├── nanopores.yaml
│   └── 5JZT_T232K_foldx.yaml
│
├── data/
│   ├── raw/
│   │   └── 5JZT/
│   │       └── 5JZT_assembly.cif
│   │
│   ├── foldx/
│   │   ├── input/
│   │   └── mutants/
│   │
│   ├── modeled/
│   │   └── 5JZT_T232K/
│   │       └── 5JZT_T232K_model.pdb
│   │
│   └── processed/
│       ├── chain_summary/
│       ├── residue_features/
│       ├── inner_residues/
│       ├── nanopore_features/
│       ├── mutation_sites/
│       └── delta_features/
│
├── outputs/
│   └── figures/
│
├── scripts/
│   ├── 00_check_config.py
│   ├── 01_parse_structure.py
│   ├── 02_extract_residue_features.py
│   ├── 03_select_inner_residues.py
│   ├── 04_build_nanopore_features.py
│   ├── 05_threshold_sensitivity.py
│   ├── 06_visualize_inner_residues.py
│   ├── 07_check_mutation_site.py
│   ├── 08_prepare_foldx_input.py
│   ├── 09_run_foldx_buildmodel.py
│   ├── 10_compare_wt_mutant_features.py
│   └── run_pipeline.py
│
├── src/
│   └── npstructfeat/
│       ├── io.py
│       ├── utils.py
│       ├── parser.py
│       ├── residue_props.py
│       ├── features.py
│       ├── geometry.py
│       ├── pore.py
│       ├── mutation.py
│       ├── foldx.py
│       └── compare.py
│
├── tests/
├── requirements.txt
├── .gitignore
└── README.md
```

---

## 4. Installation

Python 3.10+ is recommended.

Install dependencies:

```bash
pip install -r requirements.txt
```

Main dependencies:

```text
biopython
pandas
numpy
pyyaml
matplotlib
```

If the package `npstructfeat` cannot be found when running scripts, set `PYTHONPATH` temporarily in PowerShell:

```powershell
$env:PYTHONPATH="src"
```

---

## 5. Input Data

The default input structure is:

```text
data/raw/5JZT/5JZT_assembly.cif
```

Default input configuration:

```yaml
input:
  pdb_id: 5JZT
  nanopore_id: aerolysin_WT
  structure_file: data/raw/5JZT/5JZT_assembly.cif
  file_format: cif
  assembly_type: biological_assembly
```

The current example structure is the `5JZT` wild-type aerolysin biological assembly.

---

## 6. Configuration

The default configuration file is:

```text
config/default.yaml
```

Main configuration sections:

```text
project
input
structure
pore
output
save
foldx
```

The current pore-region selection uses a geometry-based approximation:

```yaml
pore:
  axis_mode: z_axis
  center_mode: xy_mean
  inner_radius_threshold: 20.0
```

The candidate pore-lining residue selection rule is:

```text
radial_distance <= inner_radius_threshold
```

where:

```text
radial_distance = radial distance from residue Cα to the estimated pore axis
```

This is a geometry-based coarse screening method and is not equivalent to rigorous pore-lining residue definitions from tools such as HOLE, CHAP, or MOLEonline.

---

## 7. Basic Pipeline

Run the full pipeline:

```bash
python scripts/run_pipeline.py
```

Or explicitly specify the configuration file:

```bash
python scripts/run_pipeline.py --config config/default.yaml
```

The full pipeline includes four steps:

1. Parse the structure and generate a chain-level summary;
2. Extract residue-level structural and physicochemical features;
3. Select candidate pore-region residues using a radial-distance threshold;
4. Aggregate candidate pore-region residues into nanopore-level structure features.

The scripts can also be run step by step:

```bash
python scripts/00_check_config.py
python scripts/01_parse_structure.py
python scripts/02_extract_residue_features.py
python scripts/03_select_inner_residues.py
python scripts/04_build_nanopore_features.py
```

---

## 8. Output Tables

With the default `5JZT` configuration, the project generates the following major output tables.

### 8.1 Chain Summary

Example path:

```text
data/processed/chain_summary/5JZT_chain_summary.csv
```

Each row corresponds to one chain and is used to check whether the biological assembly is complete.

Main fields:

```text
chain_id
total_residue_count
standard_residue_count
nonstandard_residue_count
atom_count
ca_missing_count
first_residue_number
last_residue_number
```

For the current `5JZT` example, the structure contains 7 chains, with 423 standard residues per chain.

---

### 8.2 Residue Features

Example path:

```text
data/processed/residue_features/5JZT_residue_features.csv
```

Each row corresponds to one standard amino acid residue.

Main fields:

```text
chain_id
residue_number
insertion_code
residue_name
x, y, z
has_ca
charge
hydrophobicity
is_aromatic
is_polar
```

This table can be used for residue graph construction, knowledge graph construction, spatial neighborhood analysis, and local structural feature extraction.

---

### 8.3 Candidate Inner Residues

Example path:

```text
data/processed/inner_residues/5JZT_inner_candidate_residues.csv
```

This table contains candidate pore-region residues selected using the radial-distance criterion.

Additional fields include:

```text
center_x
center_y
radial_distance
theta_rad
theta_deg
z_norm
is_inner_candidate
inner_radius_threshold
axis_mode
center_mode
```

The current default threshold is:

```text
20 Å
```

---

### 8.4 Nanopore Structure Features

Example path:

```text
data/processed/nanopore_features/5JZT_nanopore_structure_features.csv
```

This table aggregates candidate pore-region residues into one row of nanopore-level structural features.

Main features include:

```text
inner_candidate_residue_count
inner_candidate_residue_ratio
chain_count
pore_region_length_approx_A
radial_distance_mean_A
radial_distance_std_A
inner_net_charge
inner_positive_residue_count
inner_negative_residue_count
inner_mean_hydrophobicity
inner_aromatic_ratio
inner_polar_ratio
inner_residue_type_count
```

This table can be used as a nanopore structural prior for downstream tabular modeling.

---

## 9. Threshold Sensitivity Analysis

Run:

```bash
python scripts/05_threshold_sensitivity.py
```

Output:

```text
data/processed/nanopore_features/5JZT_threshold_sensitivity.csv
```

Tested thresholds:

```text
15 Å
18 Å
20 Å
22 Å
25 Å
```

The current results indicate that the 15–20 Å interval shows relatively smooth changes in candidate residue count and physicochemical properties, while thresholds of 22 Å or larger begin to include more hydrophobic and less polar residues. Therefore, 20 Å is currently used as the default threshold, with 18 Å retained as a stricter comparison setting.

---

## 10. Visualization

Run:

```bash
python scripts/06_visualize_inner_residues.py
```

Outputs:

```text
outputs/figures/5JZT_radial_distance_hist.png
outputs/figures/5JZT_z_radial_scatter.png
outputs/figures/5JZT_threshold_sensitivity_count.png
outputs/figures/5JZT_threshold_sensitivity_properties.png
```

These figures are used to inspect:

1. Whether the 20 Å threshold is a reasonable cutoff;
2. Whether candidate residues continuously cover the pore region along the z-axis;
3. Whether candidate residue counts change smoothly with threshold;
4. Whether physicochemical properties remain stable near the 20 Å threshold.

---

## 11. FoldX-based T232K Mutant Modeling

This project has been extended to support FoldX-based template mutant modeling. The current example is the `T232K` mutant of wild-type aerolysin based on the `5JZT` biological assembly.

### 11.1 Motivation

For some nanopore mutants, experimentally determined PDB structures may not be available. In such cases, arbitrary mutant PDB files should not be treated as experimentally resolved structures. They should be clearly labeled as:

```text
template-based computational mutant model
```

The current workflow is:

```text
5JZT wild-type aerolysin
        ↓
mutation site validation
        ↓
FoldX input preparation
        ↓
FoldX RepairPDB
        ↓
FoldX BuildModel
        ↓
5JZT_T232K_model.pdb
        ↓
nanopore structure feature extraction
        ↓
WT-mutant delta feature analysis
```

---

### 11.2 Step 1: Mutation Site Validation

Run:

```bash
python scripts/07_check_mutation_site.py --mutation T232K
```

Output:

```text
data/processed/mutation_sites/5JZT_T232K_site_check.csv
```

This step checks:

1. Whether residue number 232 exists in chains A-G;
2. Whether the observed residue is THR;
3. Whether the site belongs to candidate inner residues;
4. The `radial_distance`, `z`, and `z_norm` of the site.

Current result:

```text
Residue number 232 is THR in chains A-G.
All seven THR232 residues are included in the candidate inner residue set.
THR232 radial_distance ≈ 11.7–12.3 Å.
THR232 z_norm ≈ 0.632.
```

This indicates that T232K numbering is consistent with the current PDB numbering and that T232 is located in the candidate pore region under the current geometric definition.

---

### 11.3 Step 2: FoldX Input Preparation

Run:

```bash
python scripts/08_prepare_foldx_input.py --mutation T232K
```

Outputs:

```text
data/foldx/input/5JZT_WT/5JZT_assembly.pdb
data/foldx/input/5JZT_WT/5JZT_foldx_input_report.csv
```

This step converts:

```text
data/raw/5JZT/5JZT_assembly.cif
```

into a FoldX-compatible PDB file while retaining chains A-G and standard amino acid residues.

---

### 11.4 Step 3: FoldX RepairPDB

FoldX is external licensed software and is not included in this repository. Users should download FoldX separately and place it locally, for example:

```text
tools/foldx/foldx_1_20270131.exe
```

The `tools/foldx/` directory should be ignored by Git and should not be uploaded to GitHub.

Run RepairPDB:

```bash
cd data/foldx/input/5JZT_WT
path/to/foldx --command=RepairPDB --pdb=5JZT_assembly.pdb
```

Expected output:

```text
data/foldx/input/5JZT_WT/5JZT_assembly_Repair.pdb
```

`RepairPDB` repairs the wild-type structure to make it more suitable for subsequent FoldX mutation modeling.

---

### 11.5 Step 4: FoldX BuildModel

T232K should be applied to all seven chains of the aerolysin homo-heptamer.

FoldX mutation list:

```text
TA232K,TB232K,TC232K,TD232K,TE232K,TF232K,TG232K;
```

Run:

```bash
python scripts/09_run_foldx_buildmodel.py --mutation T232K
```

Outputs:

```text
data/foldx/mutants/5JZT_T232K/
data/modeled/5JZT_T232K/5JZT_T232K_model.pdb
data/foldx/mutants/5JZT_T232K/T232K_mutant_validation.csv
```

Current validation result:

```text
Residue number 232 in chains A-G has been mutated from THR to LYS.
```

---

### 11.6 Step 5: Mutant Feature Extraction

Create a separate configuration file:

```text
config/5JZT_T232K_foldx.yaml
```

This configuration should point to:

```text
data/modeled/5JZT_T232K/5JZT_T232K_model.pdb
```

Run:

```bash
python scripts/run_pipeline.py --config config/5JZT_T232K_foldx.yaml
```

Outputs:

```text
data/processed/chain_summary/5JZT_T232K_chain_summary.csv
data/processed/residue_features/5JZT_T232K_residue_features.csv
data/processed/inner_residues/5JZT_T232K_inner_candidate_residues.csv
data/processed/nanopore_features/5JZT_T232K_nanopore_structure_features.csv
```

---

### 11.7 Step 6: WT-mutant Delta Feature Analysis

Run:

```bash
python scripts/10_compare_wt_mutant_features.py
```

Output:

```text
data/processed/delta_features/5JZT_WT_vs_T232K_delta_features.csv
```

Current key differences:

| Feature | WT | T232K | Delta |
|---|---:|---:|---:|
| `inner_net_charge` | 7 | 14 | +7 |
| `inner_positive_residue_count` | 49 | 56 | +7 |
| `inner_negative_residue_count` | 42 | 42 | 0 |
| `inner_neutral_residue_count` | 429 | 422 | -7 |
| `inner_mean_hydrophobicity` | -0.7233 | -0.7663 | -0.0431 |
| `inner_polar_ratio` | 0.6327 | 0.6327 | 0 |
| `inner_aromatic_ratio` | 0.0942 | 0.0942 | 0 |
| `inner_candidate_residue_count` | 520 | 520 | 0 |
| `radial_distance_mean_A` | 13.7264 | 13.7264 | 0 |
| `pore_region_length_approx_A` | 91.9670 | 91.9670 | 0 |

Interpretation:

```text
T232K replaces seven neutral THR232 residues with seven positively charged LYS232 residues.
Therefore, the candidate pore region gains +7 net charge and +7 positively charged residues.
Because the current geometric selection is based on Cα coordinates, the FoldX point mutation does not change the candidate inner residue set.
Therefore, inner_candidate_residue_count and coarse radial-distance geometry remain unchanged.
```

---

## 12. Core Modules

```text
src/npstructfeat/io.py
    Configuration loading, path handling, and output directory preparation.

src/npstructfeat/parser.py
    Structure loading, standard residue detection, and chain-level summary.

src/npstructfeat/residue_props.py
    Amino acid physicochemical property definitions.

src/npstructfeat/features.py
    Residue-level feature extraction and nanopore-level aggregation.

src/npstructfeat/geometry.py
    Geometric center estimation, radial distance, angle, and z-normalization.

src/npstructfeat/pore.py
    Candidate pore-region residue selection.

src/npstructfeat/mutation.py
    Mutation site checking.

src/npstructfeat/foldx.py
    FoldX input preparation, BuildModel execution, and mutant validation.

src/npstructfeat/compare.py
    WT-mutant delta feature analysis.
```

---

## 13. Methodological Notes and Limitations

This project should be understood as a structural feature engineering workflow rather than a strict structural biology or molecular dynamics simulation pipeline.

Current assumptions:

1. The structure is approximately aligned along the z-axis;
2. The pore center can be approximated using the mean or median of residue Cα coordinates in the x-y plane;
3. Candidate pore-region residues can be approximated using a Cα radial-distance threshold;
4. Residue charges are simplified;
5. FoldX mutant models are template-based computational models, not experimentally determined structures.

The FoldX T232K model should be described as:

```text
template-based FoldX mutant model
```

not as:

```text
experimentally resolved T232K structure
```

The current method does not explicitly model:

```text
side-chain exposure
true pore-lining accessibility
membrane environment
water and ion distributions
applied voltage
dynamic conformational changes
pH-dependent protonation states
```

Potential future improvements include:

```text
HOLE / CHAP / MOLEonline
PCA-based pore axis estimation
side-chain orientation analysis
local mutation-region features
pore-profile features along the z-axis
Rosetta local relax
MD simulation with membrane/ions/electric field
```

---

## 14. Suggested Next Steps

Recommended future directions:

1. Add `pore-profile` features along the pore axis;
2. Build local neighborhood features around T232;
3. Compare additional mutants such as `K238Q` and `T232K/K238Q`;
4. Introduce higher-resolution aerolysin WT templates for comparison;
5. Integrate structural features with analyte ESM embeddings, experimental conditions, blockade time, and residual current data;
6. Consider KG/GNN modeling only after the tabular baseline and data schema become stable.

---

## 15. Summary

This project has completed a structural feature extraction workflow from wild-type aerolysin to a FoldX-modeled T232K mutant.

Current completed workflow:

```text
5JZT WT biological assembly
        ↓
WT feature extraction
        ↓
T232K mutation site validation
        ↓
FoldX input preparation
        ↓
FoldX RepairPDB
        ↓
FoldX BuildModel
        ↓
T232K mutant validation
        ↓
T232K feature extraction
        ↓
WT-mutant delta feature analysis
```
test 20260522 on git push
The current results indicate that, under the current geometric feature definition, T232K mainly changes the charge and hydrophilicity of the candidate pore region while leaving the coarse Cα-based pore geometry unchanged.
