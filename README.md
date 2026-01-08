# 批量分子对接脚本 (Batch Molecular Docking)

[![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)](https://www.python.org/)
[![OpenBabel](https://img.shields.io/badge/OpenBabel-3.0+-green.svg)](https://openbabel.org/)
[![Vina](https://img.shields.io/badge/AutoDock_Vina-1.2+-orange.svg)](https://vina.scripps.edu/)

## 对接环境快速创建（推荐从这里开始）

下面给出在 macOS / Linux / Windows 上快速搭建可运行本脚本的对接环境的推荐流程。

### 1️⃣ macOS / Linux：Conda + requirements.txt（推荐）

1. 安装 Miniconda / Anaconda（如已安装可跳过）：
   - 从 https://docs.conda.io/en/latest/miniconda.html 下载对应安装包
   - 按官网提示完成安装

2. 创建并激活独立环境（避免污染 base 环境）：

   ```bash
   conda create -n docking_env python=3.11
   conda activate docking_env
   ```

3. 安装 Python 依赖（使用项目自带 requirements.txt）：

   ```bash
   cd /你的/项目路径/批量对接脚本
   pip install -r requirements.txt
   ```

4. 安装 OpenBabel 与 AutoDock Vina（外部二进制）：
   - macOS：

     ```bash
     brew install open-babel vina
     ```

   - Linux (Ubuntu/Debian)：

     ```bash
     sudo apt-get install openbabel
     # Vina 请从官网 https://vina.scripps.edu/ 下载对应二进制并加入 PATH
     ```

5. （可选）安装 PDBFixer + OpenMM 用于补残基/补侧链：

   ```bash
   conda install -n docking_env -c conda-forge pdbfixer openmm
   ```

完成以上步骤后，在该环境中运行：

```bash
conda activate docking_env
cd /你的/项目路径/批量对接脚本
python batch_docking.py
```

### 2️⃣ Windows：Anaconda + 命令行工具

Windows 上建议在 **Anaconda Prompt** 中完成以下步骤：

1. 安装 Anaconda/Miniconda：
   - 从 https://www.anaconda.com/download 或 Miniconda 官网下载 Windows 版本安装包
   - 安装时勾选“将 Anaconda 添加到系统 PATH”或使用 Anaconda Prompt 运行命令

2. 创建并激活对接专用环境：

   ```bash
   conda create -n docking_env python=3.11
   conda activate docking_env
   ```

3. 在该环境中安装 Python 依赖：

   ```bash
   cd C:\\你的\\项目路径\\批量对接脚本
   pip install -r requirements.txt
   ```

4. 安装 OpenBabel：
   - 方式 A（推荐）：在当前环境中使用 conda-forge：

     ```bash
     conda install -n docking_env -c conda-forge openbabel
     ```

   - 方式 B：从 OpenBabel 官网下载 Windows 安装包并安装，确保 `obabel.exe` 所在目录加入 PATH。

5. 安装 AutoDock Vina：
   - 从官网 https://vina.scripps.edu/ 下载 Windows 版本压缩包
   - 解压后将 `vina.exe` 所在目录加入系统 PATH，或在运行脚本前将当前目录切换到包含 `vina.exe` 的路径。

6. （可选）在 Windows 上通过 Conda 安装 PDBFixer + OpenMM：

   ```bash
   conda install -n docking_env -c conda-forge pdbfixer openmm
   ```

完成以上步骤后，在 Anaconda Prompt 中执行：

```bash
conda activate docking_env
cd C:\\你的\\项目路径\\批量对接脚本
python batch_docking.py
```

### 3️⃣ Windows WSL（高级用户，可选）

如果你在 Windows 上已经使用 WSL（Windows Subsystem for Linux），可以直接按照“macOS / Linux”一节中的步骤，在 WSL 内创建 `docking_env` 环境并安装依赖，然后在 WSL 中运行脚本。这样可以获得更接近 Linux 服务器的环境行为。

##  功能介绍

一站式分子对接自动化解决方案，支持从原始数据到最终可视化的完整流程：

### ✨ 核心功能

1. **智能受体准备**
   - 自动从RCSB PDB数据库下载蛋白结构（CIF格式）
   - CIF → PDB → PDBQT全自动转换
   - 支持OpenBabel和Meeko两种处理方法
   - 原始文件与处理文件分目录管理

2. **高效配体准备**
   - SMILES直接转换为3D构象
   - 自动能量最小化优化
   - 智能处理复杂立体化学
   - 支持大分子和长链化合物

3. **批量对接引擎**
   - AutoDock Vina多线程并行计算
   - 每个对接对独立输出目录
   - 完整保存所有对接姿态
   - 自动生成受体-配体复合物PDB

4. **专业级可视化**
   - 对接结果热图（20种顶刊配色预设）
   - Top结果柱状图
   - PNG + 300dpi PDF双格式输出
   - 支持自定义配色方案

## 💻 环境要求

### 系统要求
- macOS / Linux / Windows（含 Windows WSL）
- Python 3.7+
- 8GB+ RAM（推荐16GB用于大规模对接）

### 必需软件

#### 1. Python依赖（一键安装）

项目根目录已经提供 `requirements.txt`，可直接通过以下命令安装本脚本所需的全部 Python 依赖（不包含 OpenBabel / AutoDock Vina 等外部软件）：

```bash
pip install -r requirements.txt
```

`requirements.txt` 中主要包含：
- pandas / numpy / matplotlib / seaborn / requests
- meeko（可选的受体/配体处理方案）
- rdkit-pypi（提供 RDKit 功能，供 Meeko 使用）

> 说明：PDBFixer + OpenMM 用于高级的“缺失残基/缺失侧链修补”，推荐在 Conda 环境中单独安装，因此没有写入 requirements.txt，详见下文“高级可选依赖”。

#### 2. OpenBabel

```bash
# macOS:
brew install open-babel

# Linux (Ubuntu/Debian):
sudo apt-get install openbabel
```

#### 3. AutoDock Vina

```bash
# macOS:
brew install vina

conda install -c conda-forge vina

# Linux: 从官网下载 https://vina.scripps.edu/
```

#### 4. 高级可选依赖：PDBFixer + OpenMM（用于补残基/补侧链）

如果希望启用脚本中“缺失残基/缺失侧链/缺失原子修补”功能，推荐使用 Conda 新建独立环境，并通过 conda-forge 安装：

```bash
conda create -n docking_env -c conda-forge python=3.11 pdbfixer openmm
conda activate docking_env
```

在该环境中运行脚本即可自动启用 PDBFixer 相关功能。

### 配体准备方法选择

#### 方法1：OpenBabel（推荐新手）
```bash
# 已在上面安装，无需额外操作
```

#### 方法2：Meeko（推荐进阶用户）
```bash
pip install rdkit meeko
```

## 📁 目录结构

```
批量对接脚本/
├── batch_docking.py              # 主程序脚本
├── docking_config.json           # 配置文件
├── 测试受体.csv                  # 受体输入文件
├── 测试配体.csv                  # 配体输入文件
└── docking_workspace/            # 工作目录（自动创建）
    ├── receptors/
    │   ├── raw/                  # 原始PDB下载文件（CIF格式）
    │   └── prepared/             # 处理后的受体文件（PDB、PDBQT）
    ├── ligands/
    │   ├── raw/                  # SMILES生成的3D结构（PDB）
    │   └── prepared/             # 处理后的配体文件（PDBQT）
    └── results/
        ├── 受体名_vs_配体名/        # 每个对接对一个目录
        │   ├── docked_poses.pdbqt    # 所有对接姿态
        │   ├── best_complex.pdb      # 最佳复合物结构
        │   ├── docking_log.txt       # 完整对接日志
        │   └── summary.txt           # 对接摘要
        ├── docking_results.csv       # 汇总结果表
        ├── docking_heatmap.png/pdf   # 热图
        └── top_docking_results.png/pdf  # Top结果图
```

## ⚙️ 配置文件详解

编辑 `docking_config.json` 自定义对接参数：

```json
{
  "receptor_csv": "测试受体.csv",
  "ligand_csv": "测试配体.csv",
  "work_dir": "./docking_workspace",
  "receptor_method": "openbabel",
  "box_center": [0, 0, 0],
  "box_size": [25, 25, 25],
  "exhaustiveness": 8,
  "num_modes": 9,
  "energy_range": 3,
  "cpu": 4,
  "heatmap_colormap": "RdYlGn_r",
  "receptor_cleanup": {
    "enable": true,
    "enable_altloc": true,
    "remove_waters": true,
    "remove_ions": true,
    "remove_buffers": true,
    "remove_organic_ligands": true,
    "keep_hetero_resnames": []
  },
  "pdbfixer": {
    "enable": false,
    "add_missing_residues": true,
    "add_missing_atoms": true,
    "add_hydrogens": true,
    "ph": 7.4
  }
}
```

### ✅ 输入模式：CSV vs 本地文件/压缩包（local）

脚本现在支持两种输入方式：

- **CSV 模式（默认）**：受体来自 PDB ID（脚本自动下载 CIF），配体来自 SMILES（脚本自动生成 3D 并转 PDBQT）
- **local 模式**：受体与配体都来自你已下载好的本地文件夹/单文件/压缩包（脚本自动扫描、必要时解压并转 PDBQT）

#### 1) CSV 模式（默认，不用改）

```json
{
  "receptor_input_mode": "csv",
  "ligand_input_mode": "csv",
  "receptor_csv": "测试受体.csv",
  "ligand_csv": "测试配体.csv"
}
```

#### 2) local 模式（本地文件/压缩包）

把 `receptor_input_mode` / `ligand_input_mode` 改成 `local`，并填写本地路径：

```json
{
  "receptor_input_mode": "local",
  "ligand_input_mode": "local",
  "receptor_local_path": "/abs/path/to/receptors.zip",
  "ligand_local_path": "/abs/path/to/ligands_folder",
  "local_recursive": true,
  "local_extract_archives": true,
  "receptor_local_extensions": [".pdb", ".cif", ".mmcif"],
  "ligand_local_extensions": [".pdb", ".mol", ".mol2", ".sdf", ".pdbqt"],
  "local_ligand_add_hydrogens": true,
  "local_ligand_gen3d": false,
  "local_ligand_split_multimol_sdf": false
}
```

local 模式要点：

- `receptor_local_path` / `ligand_local_path` 支持：**文件夹** / **单个文件** / **压缩包**
- 压缩包支持：`.zip` / `.tar` / `.tgz` / `.tar.gz` / `.tar.bz2` / `.tar.xz`
- `local_extract_archives=true` 时，压缩包会解压到 `docking_workspace/inputs/extracted/` 下再扫描
- 配体如果是 `.pdbqt` 会直接使用；否则会用 OpenBabel 转换为 `.pdbqt`

其中：

- `receptor_cleanup` 控制 ALTLOC 多构象处理、水/离子/缓冲剂/共晶配体的有意识删留：
  - `enable`: 总开关，控制是否启用受体清理流程
  - `enable_altloc`: 是否按占有率自动选择 ALTLOC 中最佳构象
  - `remove_waters`: 是否删除 HOH/WAT 等水分子
  - `remove_ions`: 是否删除 Na⁺/Cl⁻/Mg²⁺ 等常见无机离子
  - `remove_buffers`: 是否删除 GOL/PEG/MES/HEPES 等缓冲剂和添加剂
  - `remove_organic_ligands`: 是否删除其它有机小分子/共晶配体
  - `keep_hetero_resnames`: 三字母残基名白名单（如 `["HEM", "NAG"]`），始终保留

- `pdbfixer` 控制是否使用 PDBFixer 自动补残基/补侧链和加氢：
  - `enable`: 是否启用 PDBFixer（默认 false，需在 Conda 环境中安装依赖后再改为 true）
  - `add_missing_residues`: 是否尝试补齐缺失残基/loop
  - `add_missing_atoms`: 是否补全缺失原子
  - `add_hydrogens`: 是否根据 pH 自动加氢
  - `ph`: 加氢时使用的 pH 值

### 🎨 热图配色预设（20种顶刊级配色）

#### Nature系列
- `RdYlGn_r` - 红黄绿反转（推荐，负值=绿色=强结合）
- `RdYlBu_r` - 红黄蓝反转
- `viridis` - 感知均匀配色
- `plasma` - 高对比度配色

#### Science系列
- `coolwarm` - 冷暖对比
- `seismic` - 地震图配色
- `bwr` - 蓝白红
- `RdBu_r` - 红蓝反转

#### Cell系列
- `YlOrRd` - 黄橙红
- `YlGnBu` - 黄绿蓝
- `PuRd` - 紫红
- `BuPu` - 蓝紫

#### PNAS系列
- `Spectral_r` - 光谱反转
- `RdYlGn` - 红黄绿
- `PRGn_r` - 紫绿反转
- `PiYG_r` - 粉黄绿反转

#### 其他高级配色
- `turbo` - 高动态范围
- `cividis` - 色盲友好
- `inferno` - 火焰
- `magma` - 岩浆

💡 **使用方法**：将配色名称复制到配置文件的 `heatmap_colormap` 字段即可

### ⚠️ 重要：设置对接盒子

对接盒子决定了配体在受体的哪个区域进行对接。有几种设置方法：

#### 方法1：如果知道活性位点坐标
```json
{
  "box_center": [15.2, -3.5, 22.8],  # 活性位点坐标
  "box_size": [20, 20, 20]            # 覆盖活性位点的盒子
}
```

#### 方法2：全蛋白对接（不确定活性位点时）
```json
{
  "box_center": [0, 0, 0],
  "box_size": [40, 40, 40]            # 较大的盒子覆盖整个蛋白
}
```

#### 如何获取活性位点坐标？

1. **使用PyMOL**：
   ```
   - 打开受体PDB文件
   - 鼠标点击活性位点残基
   - 在底部查看坐标
   ```

2. **使用Chimera**：
   ```
   - Tools > Structure Analysis > Axes/Planes/Centroids
   - 选择活性位点残基
   - 查看质心坐标
   ```

3. **根据已知配体**：
   如果PDB文件包含配体，使用其坐标作为中心

## 🚀 快速开始

### 1️⃣ 准备输入文件

你可以选择 **CSV 模式** 或 **local 模式**：

#### 方案 A：CSV 模式（默认）

**测试受体.csv** 格式：
```csv
PRO,ID
Thioredoxin Reductase 1,3EAN
Thrombin,2ZO3
Dihydrofolate Reductase,1S3U
```

**测试配体.csv** 格式：
```csv
CHE,SMILES
Inermine,O1[C@@H]2[C@H](C3=CC4=C(OCO4)C=C13)COC1=C2C=CC(=C1)O
Medicarpin,COC1=CC2=C(C=C1)[C@@H]3COC4=C([C@@H]3O2)C=CC(=C4)O
```

#### 方案 B：local 模式（本地文件夹/压缩包）

准备两份输入（任意一种组合即可）：

- 受体：`pdb/cif/mmcif` 文件夹或压缩包
- 配体：`pdb/mol/mol2/sdf/pdbqt` 文件夹或压缩包

然后在 `docking_config.json` 中设置：

```json
{
  "receptor_input_mode": "local",
  "ligand_input_mode": "local",
  "receptor_local_path": "/abs/path/to/receptors.zip",
  "ligand_local_path": "/abs/path/to/ligands.zip"
}
```

### 2️⃣ 配置参数（可选）

根据需要编辑 `docking_config.json`：
- 调整对接盒子参数
- 选择喜欢的热图配色
- 设置CPU核心数

### 3️⃣ 运行脚本

```bash
cd /Users/dingtunan/生信分析/批量对接脚本
python batch_docking.py
```

脚本将自动完成：
1. ✓ 下载受体结构（4个受体）
2. ✓ 准备受体和配体文件（18个配体）
3. ✓ 执行批量对接（72个对接任务）
4. ✓ 生成可视化图表

### 4️⃣ 查看结果

#### 📊 汇总结果
- `docking_results.csv` - 所有对接数据（可用Excel打开）
- `docking_heatmap.png/pdf` - 对接亲和力热图
- `top_docking_results.png/pdf` - 最佳结果排名

#### 📂 详细结果
每个对接对的独立目录包含：
- `docked_poses.pdbqt` - 9个对接姿态
- `best_complex.pdb` - 最佳复合物（可用PyMOL/Chimera查看）
- `docking_log.txt` - 完整日志
- `summary.txt` - 对接摘要

## 📈 结果解读

### 对接亲和力（Binding Affinity）

- 单位：kcal/mol
- **数值越负，结合越强**
- 典型范围：
  - ≤ -10 kcal/mol：非常强的结合
  - -8 到 -10：强结合
  - -6 到 -8：中等结合
  - ≥ -6：弱结合

### 热图颜色（以RdYlGn_r为例）

- 🟢 **绿色**：强结合（≤ -8 kcal/mol）
- 🟡 **黄色**：中等结合（-6 到 -8 kcal/mol）
- 🔴 **红色**：弱结合（≥ -6 kcal/mol）

> 💡 不同配色方案颜色含义可能不同，但数值始终表示结合亲和力

## 常见问题

### 1. 下载受体失败

**原因**：网络连接或PDB ID错误

**解决**：
- 检查PDB ID是否正确
- 访问 https://www.rcsb.org/ 确认ID存在
- 检查网络连接

### 2. 配体准备失败

**OpenBabel方法失败**：
```bash
# 检查OpenBabel安装
obabel --version

# 测试SMILES转换
obabel -:"CCO" -opdb -O test.pdb --gen3d
```

**Meeko方法失败**：
```bash
# 检查RDKit和Meeko安装
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "from meeko import MoleculePreparation; print('Meeko OK')"

# 如果缺少，安装：
pip install rdkit meeko
```

### 3. Vina对接失败

```bash
# 检查Vina安装
vina --help

# macOS安装Vina
brew install vina

# 或手动下载：https://vina.scripps.edu/
```

### 4. 受体准备工具选择

**🔍 脚本智能检测机制**：

脚本会自动检测系统中可用的受体准备工具，并按以下优先级使用：

1. **MGLTools的prepare_receptor4**（最专业，如果已安装）
   ```bash
   # Windows/Linux安装
   # 下载：http://mgltools.scripps.edu/downloads
   
   # 验证安装
   prepare_receptor4 -h
   ```

2. **OpenBabel**（可靠的备选方案，适合所有平台）
   ```bash
   # 已安装则自动使用
   obabel --version
   ```

✅ **自动回退机制**：
- 已安装MGLTools → 自动使用`prepare_receptor4`
- 未安装MGLTools → 自动回退到OpenBabel
- 无需修改配置，脚本根据环境自动选择

⚠️ **注意**：macOS ARM架构不支持MGLTools，但OpenBabel效果同样专业

### 5. 内存不足

如果受体或配体数量很多，可以分批处理：
- 将CSV文件拆分为多个小文件
- 分别运行脚本
- 合并结果

### 6. 对接盒子设置不当

如果对接结果不理想，可能是盒子设置问题：
- 确保盒子包含活性位点
- 增大盒子尺寸（如从20x20x20改为30x30x30）
- 如果不确定，使用全蛋白对接（40x40x40）

## 📤 输出格式说明

### 图片格式
所有图表同时输出两种格式：
- **PNG格式** - 用于预览、PPT插入（300 DPI高清）
- **PDF格式** - 用于论文投稿、矢量编辑（300 DPI）

### 为什么需要PDF？
- ✓ 无损缩放，适合任意尺寸打印
- ✓ 期刊投稿首选格式
- ✓ 可用Adobe Illustrator/Inkscape编辑
- ✓ 文件更小，便于存档

## 🔧 进阶使用

### 自定义单个受体的对接盒子

如果不同受体需要不同的盒子参数，可以修改脚本中的 `run_vina_docking` 方法，为每个受体设置特定参数。

### 并行加速

对于大量对接任务，可以使用Python的 `multiprocessing` 模块并行处理：

```python
from multiprocessing import Pool

# 在batch_docking方法中
with Pool(processes=4) as pool:
    results = pool.starmap(run_docking_single, docking_tasks)
```

### 更精确的能量最小化

在配体准备中增加优化步数：

```python
# OpenBabel方法
--minimize --steps 1000 --ff MMFF94

# RDKit/Meeko方法
AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
```

## 引用

如果使用此脚本发表论文，请引用：

- **AutoDock Vina**: Trott, O., & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. *Journal of computational chemistry*, 31(2), 455-461.

- **OpenBabel**: O'Boyle, N. M., et al. (2011). Open Babel: An open chemical toolbox. *Journal of cheminformatics*, 3(1), 33.

- **RDKit**: RDKit: Open-source cheminformatics. http://www.rdkit.org

- **Meeko**: Forli, S., et al. (2016). Computational protein–ligand docking and virtual drug screening with the AutoDock suite. *Nature protocols*, 11(5), 905-919.

## 许可证

本脚本采用 MIT 许可证。

## 联系方式

如有问题或建议，请提交 Issue 或联系作者。

---

**祝您对接顺利！** 🧬💊
