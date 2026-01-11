# Batch Molecular Docking

[![README-English](https://img.shields.io/badge/README-English-blue)](README_EN.md)
[![README-%E4%B8%AD%E6%96%87](https://img.shields.io/badge/README-%E4%B8%AD%E6%96%87-blue)](README.md)

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![OpenBabel](https://img.shields.io/badge/OpenBabel-3.x-green.svg)](https://openbabel.org/)
[![AutoDock%20Vina](https://img.shields.io/badge/AutoDock%20Vina-1.2%2B-orange.svg)](https://vina.scripps.edu/)

An end-to-end toolkit for: receptor/ligand preparation → batch docking → results summarization → visualization.

- Basic pipeline entry: `batch_docking.py` (config: `docking_config.json`)
- Advanced pipeline entry: `batch_docking_advanced.py` (config: `docking_config_advanced.json`)

Plots are exported as both PNG and PDF by default, and plot text is English by default.

---

## 1. Environment setup (macOS / Windows / WSL)

General principle:
- Install Python dependencies via `requirements.txt`
- Install external tools (Vina/OpenBabel) via system package managers or official binaries

### 1.1 macOS (Recommended: Conda + Homebrew)

1) Create a Conda environment and install Python deps:

```bash
conda create -n docking_env python=3.11 -y
conda activate docking_env
cd /path/to/BatchDocking
pip install -r requirements.txt
```

2) Install external tools:

```bash
brew install open-babel vina
```

3) (Optional) Install PDBFixer + OpenMM (structure fixing / adding hydrogens):

```bash
conda install -n docking_env -c conda-forge pdbfixer openmm -y
```

4) Verify:

```bash
vina --help | head -n 2
obabel -V
python -c "import pandas,numpy,rdkit,gemmi; print('OK')"
```

### 1.2 Windows (Native)

1) In Anaconda Prompt:

```bash
conda create -n docking_env python=3.11 -y
conda activate docking_env
cd C:\\path\\to\\BatchDocking
pip install -r requirements.txt
```

2) Install OpenBabel (one of the following):
- Recommended: `conda install -n docking_env -c conda-forge openbabel`
- Or install OpenBabel for Windows and add `obabel.exe` to PATH

3) Install AutoDock Vina:
- Download from https://vina.scripps.edu/
- Add `vina.exe` to PATH

### 1.3 Windows WSL (Recommended: Linux-style)

Inside WSL, follow Linux steps:

```bash
conda create -n docking_env python=3.11 -y
conda activate docking_env
cd /path/to/BatchDocking
pip install -r requirements.txt

# Ubuntu/Debian example:
sudo apt-get update
sudo apt-get install -y openbabel
# Vina: download official binary and add to PATH
```

---

## 2. Dependencies and external tools

### 2.1 Python dependencies

Install from `requirements.txt`:

```bash
pip install -r requirements.txt
```

Package roles (high-level):
- `pandas/numpy`: result tables and matrix ops
- `matplotlib/seaborn`: heatmap and top hits plots (PNG+PDF)
- `requests`: download structures from RCSB (CSV mode)
- `rdkit-pypi`: 3D generation and ligand enumeration (advanced pipeline)
- `meeko`: optional PDBQT preparation pipeline
- `gemmi`: mmCIF/CIF parsing & conversion (also useful as fallback)
- `biopython`: supplementary structure utilities

### 2.2 External tools (required)

The scripts call these commands directly (must be in PATH):
- `vina` (AutoDock Vina)
- `obabel` (OpenBabel)

---

## 3. Project structure and outputs

Typical layout (example: `work_dir=./hif1a`):

```
hif1a/
  receptors/
    raw/
    prepared/
  ligands/
    raw/
    prepared/
  results/
    <Receptor[_ID]>_vs_<Ligand>/
      pocket_*/
      docked_poses.pdbqt
      best_complex.pdb
      docking_log.txt
      summary.txt
    docking_results.csv
    docking_heatmap.png/.pdf
    top_docking_results.png/.pdf
```

---

## 4. Basic pipeline: batch_docking.py

### 4.1 What it does

Designed for robustness and reproducibility:
- Receptor preparation: mmCIF/CIF/PDB → PDBQT
- Ligand preparation: SMILES → 3D → PDBQT
- Batch docking via Vina: per pair output folder + logs + complex PDB
- Summary: `docking_results.csv`
- Visualization: heatmap + top hits (PNG + PDF)

### 4.2 Run

```bash
python batch_docking.py
```

By default it reads `docking_config.json`.

### 4.3 docking_config.json key fields

| Field | Meaning | Example |
|---|---|---|
| `work_dir` | Working directory | `./hif1a` |
| `receptor_input_mode` | Receptor input mode | `csv` / `local` |
| `ligand_input_mode` | Ligand input mode | `csv` / `local` |
| `receptor_csv` | Receptor CSV (csv mode) | `hk2测试.csv` |
| `ligand_csv` | Ligand CSV (csv mode) | `瑞香素测试.csv` |
| `receptor_local_path` | Local receptor path (local mode) | `/path/to/receptors/raw` |
| `ligand_local_path` | Local ligand path (local mode) | `/path/to/ligands/raw` |
| `box_center` | Docking box center | `auto` or `[x,y,z]` |
| `box_size` | Docking box size (Å) | `[40,40,40]` |
| `exhaustiveness` | Search exhaustiveness | higher = slower |
| `num_modes` | Number of output poses | 9 |
| `cpu` | CPU cores | 4 |

CSV column conventions (recommended):
- Receptor CSV: at least `ID` (PDB ID), recommended `PRO` (protein name)
- Ligand CSV: at least `SMILES`, recommended `CHE` (ligand name)

---

## 5. Advanced pipeline: batch_docking_advanced.py

### 5.1 What’s added

On top of the basic pipeline:

1) Ligand enumeration
- Generate multiple tautomers/conformers for the same SMILES
- Automatically keep the best variant per receptor–ligand pair

2) Multi-pocket docking (optional)
- Predict multiple candidate pockets using different algorithms
- Dock each pocket and keep the best affinity

3) Storage naming vs plotting labels
- Output folders are stored as `ProteinName + ID` to avoid collisions
- Plots show only `ProteinName` as receptor label (same names are aggregated by best affinity)

### 5.2 Run

```bash
python batch_docking_advanced.py
```

Config loading order: prefer `docking_config_advanced.json`, fall back to `docking_config.json`.

### 5.3 docking_config_advanced.json key fields

Ligand enumeration:

| Field | Meaning | Typical |
|---|---|---|
| `ligand_enumeration_enable` | Enable enumeration | `true` |
| `ligand_enum_max_tautomers` | Max tautomers | 4–12 |
| `ligand_enum_max_conformers` | Max conformers per tautomer | 5–20 |
| `ligand_enum_max_variants` | Max variants kept per ligand | 12–30 |
| `ligand_enum_random_seed` | RNG seed | fixed for reproducibility |
| `ligand_ph` | pH | 7.4 |

Multi-pocket docking:

| Field | Meaning | Typical |
|---|---|---|
| `pocket_enable` | Enable multi-pocket | default `false` |
| `pocket_algorithms` | Algorithms | `cocrystal/surface_random/receptor_center` |
| `pocket_max_pockets` | Max pockets per receptor | 3–8 |
| `pocket_box_size` | Box size for pocket docking (Å) | ≥ `box_size` |

Key tuning params for `surface_random`:
- `pocket_random_points`: random sampling points (bigger = slower)
- `pocket_top_points`: kept candidate points
- `pocket_dist_min/pocket_dist_max`: distance window to protein surface (Å)
- `pocket_cluster_radius`: clustering radius for de-duplication (Å)

---

## 6. Recommended workflow

1) Start with 1 receptor × 1 ligand to validate environment
2) Scale up gradually
3) Enable multi-pocket only when needed (`pocket_enable: true`)

---

## 7. FAQ

- `vina` not found: ensure it’s in PATH (`where vina` on Windows; `command -v vina` on macOS/Linux)
- `obabel` not found: ensure it’s in PATH (`where obabel` / `command -v obabel`)
- Too slow: reduce `exhaustiveness` or decrease `pocket_max_pockets` / `pocket_random_points`

