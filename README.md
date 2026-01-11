# 批量分子对接脚本（Batch Molecular Docking）

[![README-中文](https://img.shields.io/badge/README-%E4%B8%AD%E6%96%87-blue)](README.md)
[![README-English](https://img.shields.io/badge/README-English-blue)](README_EN.md)

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](https://www.python.org/)
[![OpenBabel](https://img.shields.io/badge/OpenBabel-3.x-green.svg)](https://openbabel.org/)
[![AutoDock%20Vina](https://img.shields.io/badge/AutoDock%20Vina-1.2%2B-orange.svg)](https://vina.scripps.edu/)

这是一套“受体/配体准备 → 批量对接 → 结果汇总 → 可视化导出”的一站式脚本。

- 基础版入口：`batch_docking.py`（配置：`docking_config.json`）
- 高级版入口：`batch_docking_advanced.py`（配置：`docking_config_advanced.json`）

输出图表默认保存为 PNG + PDF 双格式，且图中文字默认英文。

---

## 1. 三种系统环境配置

总体原则：
- Python 依赖使用 `requirements.txt`
- 外部工具（Vina/OpenBabel）使用系统包管理器或官方二进制

### 1.1 macOS（推荐：Conda + Homebrew）

1）创建环境并安装 Python 依赖：

```bash
conda create -n docking_env python=3.11 -y
conda activate docking_env
cd /你的/项目路径/批量对接脚本
pip install -r requirements.txt
```

2）安装外部工具：

```bash
brew install open-babel vina
```

3）（可选）安装 PDBFixer + OpenMM（用于自动补残基/补原子/加氢）：

```bash
conda install -n docking_env -c conda-forge pdbfixer openmm -y
```

4）验证安装：

```bash
vina --help | head -n 2
obabel -V
python -c "import pandas,numpy,rdkit,gemmi; print('OK')"
```

### 1.2 Windows（原生：Anaconda Prompt + 外部工具）

1）在 Anaconda Prompt 中创建环境并安装 Python 依赖：

```bash
conda create -n docking_env python=3.11 -y
conda activate docking_env
cd C:\\你的\\项目路径\\批量对接脚本
pip install -r requirements.txt
```

2）安装 OpenBabel（二选一）：
- 推荐：`conda install -n docking_env -c conda-forge openbabel`
- 或者：安装 OpenBabel Windows 安装包，并把 `obabel.exe` 加入 PATH

3）安装 AutoDock Vina：
- 从官网下载安装包：https://vina.scripps.edu/
- 将 `vina.exe` 所在目录加入 PATH

### 1.3 Windows WSL（推荐：按 Linux 方式）

WSL 内按 Linux 步骤安装即可（建议使用 Miniconda/Conda）：

```bash
conda create -n docking_env python=3.11 -y
conda activate docking_env
cd /你的/项目路径/批量对接脚本
pip install -r requirements.txt

# Ubuntu/Debian 示例：
sudo apt-get update
sudo apt-get install -y openbabel
# Vina 建议从官网下载二进制并加入 PATH
```

---

## 2. 依赖安装与工具安装

### 2.1 Python 依赖（requirements.txt）

仓库根目录提供 `requirements.txt`，直接安装：

```bash
pip install -r requirements.txt
```

依赖用途概览：
- `pandas/numpy`：批量结果整理、矩阵化
- `matplotlib/seaborn`：热图与 Top 结果图（PNG+PDF）
- `requests`：从 RCSB 下载结构（CSV 模式）
- `rdkit-pypi`：配体三维构象与枚举（高级版核心）
- `meeko`：可选的 PDBQT 处理链路
- `gemmi`：mmCIF/CIF 解析与转换（也可作为 CIF→PDB 备用）
- `biopython`：结构处理能力补充

### 2.2 外部工具（必须）

脚本会直接调用以下命令，请确保在 PATH 中可用：
- `vina`（AutoDock Vina）
- `obabel`（OpenBabel）

---

## 3. 项目结构与输出

典型目录（以 `work_dir=./hif1a` 为例）：

```
hif1a/
  receptors/
    raw/        # 原始受体文件（csv下载或local输入）
    prepared/   # 处理后的受体（含 pdbqt/pdb 等）
  ligands/
    raw/        # 原始配体（csv/SMILES生成或local输入）
    prepared/   # 处理后的配体（含 pdbqt；高级版可能有多变体）
  results/
    <Receptor[_ID]>_vs_<Ligand>/
      pocket_*/               # 多口袋时存在（高级版）
      docked_poses.pdbqt
      best_complex.pdb
      docking_log.txt
      summary.txt
    docking_results.csv
    docking_heatmap.png/.pdf
    top_docking_results.png/.pdf
```

---

## 4. 基础版：batch_docking.py

### 4.1 功能概述

基础版目标是“稳、可复现、少参数”，适合：
- 受体/配体数量中等
- 不需要对同一配体做互变体/构象枚举
- 对接盒参数相对固定（或用 `auto` 取受体几何中心）

核心能力：
- 受体准备：mmCIF/CIF/PDB → PDBQT
- 配体准备：SMILES → 3D → PDBQT
- Vina 批量对接：每个受体-配体独立目录、保存日志与复合物
- 结果汇总：输出 `docking_results.csv`
- 可视化：热图 + Top 结果图（PNG + PDF）

### 4.2 运行方式

```bash
python batch_docking.py
```

脚本默认读取 `docking_config.json`。

### 4.3 docking_config.json 关键配置说明

常用字段：

| 字段 | 说明 | 示例 |
|---|---|---|
| `work_dir` | 工作目录 | `./hif1a` |
| `receptor_input_mode` | 受体输入模式 | `csv` / `local` |
| `ligand_input_mode` | 配体输入模式 | `csv` / `local` |
| `receptor_csv` | 受体CSV（csv模式） | `受体测试.csv` |
| `ligand_csv` | 配体CSV（csv模式） | `配体测试.csv` |
| `receptor_local_path` | 本地受体路径（local模式） | `/path/to/receptors/raw` |
| `ligand_local_path` | 本地配体路径（local模式） | `/path/to/ligands/raw` |
| `box_center` | 对接盒中心 | `auto` 或 `[x,y,z]` |
| `box_size` | 对接盒尺寸(Å) | `[40,40,40]` |
| `exhaustiveness` | 搜索强度 | 8~64（越大越慢） |
| `num_modes` | 输出构象数 | 9 |
| `energy_range` | 能量窗口 | 3 |
| `cpu` | 并行核数 | 4 |

CSV 列名约定（建议遵守，避免脚本无法识别）：
- 受体 CSV：至少包含 `ID`（PDB ID），建议包含 `PRO`（蛋白名）
- 配体 CSV：至少包含 `SMILES`，建议包含 `CHE`（配体名）

---

## 5. 高级版：batch_docking_advanced.py

### 5.1 功能概述

高级版在基础版上增加两类“更接近实战”的能力：

1）配体枚举（Ligand Enumeration）
- 对同一 SMILES 自动生成互变异构体与构象
- 对每个受体-配体组合，自动选择最优变体作为最终结果

2）多口袋对接（Multi-pocket Docking，可选）
- 对同一受体用多种算法预测多个候选口袋
- 对每个口袋分别对接，自动取最优 Affinity

3）结果命名与绘图聚合
- 结果目录按 `蛋白名 + ID` 落盘，避免同名受体互相覆盖
- 绘图（热图/Top）仅用“蛋白名”作为受体标签展示（同名受体按最优 Affinity 聚合）

### 5.2 运行方式

```bash
python batch_docking_advanced.py
```

读取配置规则：优先 `docking_config_advanced.json`，若不存在则回退 `docking_config.json`。

### 5.3 docking_config_advanced.json 关键配置说明

#### 5.3.1 配体枚举相关

| 字段 | 说明 | 建议 |
|---|---|---|
| `ligand_enumeration_enable` | 是否启用枚举 | `true` |
| `ligand_enum_max_tautomers` | 最大互变体数 | 4~12 |
| `ligand_enum_max_conformers` | 每互变体最大构象数 | 5~20 |
| `ligand_enum_max_variants` | 每配体最多保留变体数 | 12~30 |
| `ligand_enum_random_seed` | 随机种子 | 固定可复现 |
| `ligand_ph` | pH | 7.4 |

#### 5.3.2 多口袋对接相关（可选）

| 字段 | 说明 | 建议 |
|---|---|---|
| `pocket_enable` | 是否启用多口袋 | 默认 `false` |
| `pocket_algorithms` | 口袋算法组合 | `cocrystal/surface_random/receptor_center` |
| `pocket_max_pockets` | 每受体最多口袋数 | 3~8 |
| `pocket_box_size` | 口袋对接盒尺寸(Å) | 一般 ≥ `box_size` |

`surface_random` 算法常调参数：
- `pocket_random_points`：随机采样点数（越大越慢）
- `pocket_top_points`：保留高分候选点数
- `pocket_dist_min/pocket_dist_max`：候选点与蛋白表面距离窗口（Å）
- `pocket_cluster_radius`：候选中心去重半径（Å）

---

## 6. 推荐工作流（更稳、更好排错）

1）先用 1 个受体 × 1 个配体跑通流程（确认 Vina/OpenBabel/依赖都正常）
2）再逐步放大规模（受体/配体数量）
3）需要多口袋时再打开：`pocket_enable: true`
4）受体同名但不同 PDB：
- 目录落盘：`PRO_ID_vs_Ligand/...`
- 绘图展示：只显示 `PRO`

---

## 7. FAQ

- `vina` 找不到：确认 PATH（macOS/Linux：`command -v vina`；Windows：`where vina`）
- `obabel` 找不到：确认 PATH（macOS/Linux：`command -v obabel`；Windows：`where obabel`）
- 运行慢：降低 `exhaustiveness` 或减少 `pocket_max_pockets` / `pocket_random_points`

