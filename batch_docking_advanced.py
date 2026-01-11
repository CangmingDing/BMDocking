#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
批量分子对接脚本
功能：
1. 从PDB数据库下载受体蛋白结构
2. 准备受体和配体文件（支持OpenBabel和Meeko两种方法）
3. 使用AutoDock Vina进行批量对接
4. 生成对接结果热图
"""

import os
import sys
import pandas as pd
import requests
import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Tuple, Dict, Optional
import shutil
import zipfile
import tarfile
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class MolecularDocking:
    """分子对接主类"""
    
    def __init__(self, config_file: str = "docking_config_advanced.json"):
        """
        初始化对接系统
        
        Args:
            config_file: 配置文件路径
        """
        self.config = self.load_config(config_file)
        self.work_dir = Path(self.config.get('work_dir', './docking_workspace'))

        self._receptor_center_cache: Dict[str, List[float]] = {}
        self._receptor_size_cache: Dict[str, List[float]] = {}
        self._receptor_pdb_cache: Dict[str, str] = {}
        self._receptor_pockets_cache: Dict[str, List[Dict]] = {}
        
        # 创建分层目录结构
        self.receptor_raw_dir = self.work_dir / 'receptors' / 'raw'
        self.receptor_prepared_dir = self.work_dir / 'receptors' / 'prepared'
        self.ligand_raw_dir = self.work_dir / 'ligands' / 'raw'
        self.ligand_prepared_dir = self.work_dir / 'ligands' / 'prepared'
        self.results_dir = self.work_dir / 'results'
        
        # 创建所有工作目录
        for dir_path in [self.receptor_raw_dir, self.receptor_prepared_dir, 
                         self.ligand_raw_dir, self.ligand_prepared_dir, self.results_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def load_config(self, config_file: str) -> dict:
        """加载配置文件"""
        if os.path.exists(config_file):
            with open(config_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        else:
            # 默认配置（新增受体清理与PDBFixer相关参数）
            return {
                'work_dir': './docking_workspace',
                # 输入模式：csv 或 local
                'receptor_input_mode': 'csv',
                'ligand_input_mode': 'csv',
                # local 模式下：可以填“文件夹路径”或“压缩包路径”或“单个结构文件路径”
                'receptor_local_path': '',
                'ligand_local_path': '',
                'local_recursive': True,
                'local_extract_archives': True,
                # local 配体转换参数
                'local_ligand_add_hydrogens': True,
                'local_ligand_gen3d': False,
                'local_ligand_split_multimol_sdf': False,
                'receptor_method': 'openbabel',
                'receptor_cif_to_pdb_fallback_enable': True,
                'receptor_cif_to_pdb_preferred_tool': 'gemmi',
                'receptor_cif_to_pdb_openbabel_timeout_sec': 60,
                'ligand_method': 'meeko',  # 'openbabel' or 'meeko'
                'box_center': [0, 0, 0],
                'box_size': [20, 20, 20],
                'exhaustiveness': 8,
                'num_modes': 9,
                'energy_range': 3,
                # 受体结构清理相关配置
                'receptor_cleanup': {
                    'enable': True,
                    # ALTLOC 多构象处理开关
                    'enable_altloc': True,
                    # 是否删除水分子
                    'remove_waters': True,
                    # 是否删除无机离子（Na+/Cl-/Mg2+ 等）
                    'remove_ions': True,
                    # 是否删除常见缓冲剂小分子（GOL/PEG/MES 等）
                    'remove_buffers': True,
                    # 是否删除其它有机小分子/共晶配体
                    'remove_organic_ligands': True,
                    # 需要强制保留的异源分子三字母残基名列表，例如 ["HEM", "NAG"]
                    'keep_hetero_resnames': []
                },
                # PDBFixer 相关配置
                'pdbfixer': {
                    'enable': True,
                    # 是否尝试补全部分缺失残基（loop）信息
                    'add_missing_residues': True,
                    # 是否添加缺失侧链/主链原子
                    'add_missing_atoms': True,
                    # 是否根据 pH 添加氢原子
                    'add_hydrogens': True,
                    # 用于加氢的 pH 值
                    'ph': 7.4
                }
            }

    def _normalize_local_source(self, source_path: str) -> Path:
        """规范化本地输入路径，支持相对路径与用户目录

        Args:
            source_path: 文件/文件夹/压缩包路径

        Returns:
            规范化后的绝对路径
        """
        expanded = os.path.expanduser(str(source_path)).strip()
        return Path(expanded).resolve()

    def _extract_archive_to_dir(self, archive_path: Path, extract_root: Path) -> Path:
        """将压缩包解压到工作目录下的指定位置

        Args:
            archive_path: 压缩包路径（zip/tar/tar.gz/tgz/tar.bz2 等）
            extract_root: 解压输出根目录

        Returns:
            解压后的目录路径
        """
        extract_root.mkdir(parents=True, exist_ok=True)
        ts = datetime.now().strftime('%Y%m%d_%H%M%S')
        out_dir = extract_root / f"{archive_path.stem}_{ts}"
        out_dir.mkdir(parents=True, exist_ok=True)

        lower_name = archive_path.name.lower()
        if lower_name.endswith('.zip'):
            with zipfile.ZipFile(str(archive_path), 'r') as zf:
                zf.extractall(str(out_dir))
            return out_dir

        if (
            lower_name.endswith('.tar')
            or lower_name.endswith('.tgz')
            or lower_name.endswith('.tar.gz')
            or lower_name.endswith('.tar.bz2')
            or lower_name.endswith('.tar.xz')
        ):
            with tarfile.open(str(archive_path), 'r:*') as tf:
                tf.extractall(str(out_dir))
            return out_dir

        raise ValueError(f"不支持的压缩包格式: {archive_path}")

    def _resolve_local_input_root(self, source: str, kind: str) -> Path:
        """解析本地输入：支持文件夹/单文件/压缩包

        Args:
            source: 用户填写的本地路径
            kind: 'receptors' 或 'ligands'

        Returns:
            可遍历的目录路径（如果是单文件则返回其父目录，遍历时会仅选择该文件）
        """
        source_path = self._normalize_local_source(source)
        if not source_path.exists():
            raise FileNotFoundError(f"本地输入路径不存在: {source_path}")

        if source_path.is_dir():
            return source_path

        if source_path.is_file():
            if not self.config.get('local_extract_archives', True):
                return source_path.parent

            lower_name = source_path.name.lower()
            is_archive = any(
                lower_name.endswith(suffix)
                for suffix in ['.zip', '.tar', '.tgz', '.tar.gz', '.tar.bz2', '.tar.xz']
            )
            if is_archive:
                extract_root = self.work_dir / 'inputs' / 'extracted' / kind
                return self._extract_archive_to_dir(source_path, extract_root)
            return source_path.parent

        raise ValueError(f"无法识别的输入路径: {source_path}")

    def _collect_files_by_ext(self, root: Path, exts: List[str], recursive: bool) -> List[Path]:
        """在目录中按扩展名收集文件

        Args:
            root: 搜索根目录
            exts: 扩展名列表（如 ['.pdb', '.cif']）
            recursive: 是否递归搜索

        Returns:
            匹配到的文件路径列表
        """
        normalized_exts = set(e.lower() if e.startswith('.') else f".{e.lower()}" for e in exts)
        pattern = '**/*' if recursive else '*'
        files = []
        for p in root.glob(pattern):
            if p.is_file() and p.suffix.lower() in normalized_exts:
                files.append(p)
        files.sort(key=lambda x: x.name.lower())
        return files

    def _safe_stem(self, file_path: Path) -> str:
        """根据文件名生成安全的条目名称"""
        stem = file_path.stem
        safe = "".join(c for c in stem if c.isalnum() or c in (' ', '-', '_')).strip()
        return safe if safe else stem

    def _cache_receptor_box_from_pdb(self, receptor_pdb: Path, receptor_pdbqt: Path) -> None:
        center_cfg = self.config.get('box_center', 'auto')
        size_cfg = self.config.get('box_size', None)
        center_mode = str(center_cfg).strip().lower() if isinstance(center_cfg, str) else ''
        size_mode = str(size_cfg).strip().lower() if isinstance(size_cfg, str) else ''

        use_ligand_center = bool(self.config.get('box_auto_use_cocrystal_ligand', True)) and center_mode in (
            'auto',
            'auto_ligand',
            'ligand',
            'cocrystal',
        )
        use_ligand_size = size_mode in ('auto', 'auto_ligand', 'ligand', 'cocrystal')
        if not (use_ligand_center or use_ligand_size):
            return

        inferred = self._infer_box_from_cocrystal_ligand(receptor_pdb)
        if inferred is None:
            return
        center, size = inferred

        cache_key = os.path.abspath(str(receptor_pdbqt))
        if use_ligand_center:
            self._receptor_center_cache[cache_key] = center
        if use_ligand_size:
            self._receptor_size_cache[cache_key] = size

    def _read_receptor_protein_atoms_from_pdb(self, receptor_pdb: str) -> Optional[np.ndarray]:
        try:
            coords: List[Tuple[float, float, float]] = []
            with open(receptor_pdb, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if not line.startswith('ATOM'):
                        continue
                    if len(line) < 54:
                        continue
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except Exception:
                        continue
                    coords.append((x, y, z))
            if not coords:
                return None
            return np.array(coords, dtype=float)
        except Exception:
            return None

    def _cluster_pocket_centers(self, centers: List[List[float]], min_dist: float) -> List[List[float]]:
        kept: List[List[float]] = []
        for c in centers:
            ok = True
            for k in kept:
                dx = c[0] - k[0]
                dy = c[1] - k[1]
                dz = c[2] - k[2]
                if (dx * dx + dy * dy + dz * dz) <= (min_dist * min_dist):
                    ok = False
                    break
            if ok:
                kept.append(c)
        return kept

    def _predict_pockets_surface_random(self, receptor_pdb: str) -> List[Dict]:
        atoms = self._read_receptor_protein_atoms_from_pdb(receptor_pdb)
        if atoms is None or atoms.size == 0:
            return []

        n_points = int(self.config.get('pocket_random_points', 8000))
        top_keep = int(self.config.get('pocket_top_points', 400))
        dmin = float(self.config.get('pocket_dist_min', 2.5))
        dmax = float(self.config.get('pocket_dist_max', 6.5))
        chunk = int(self.config.get('pocket_distance_chunk', 200))
        cluster_radius = float(self.config.get('pocket_cluster_radius', 8.0))
        max_pockets = int(self.config.get('pocket_max_pockets', 5))
        seed = int(self.config.get('pocket_random_seed', 42))

        min_xyz = atoms.min(axis=0)
        max_xyz = atoms.max(axis=0)
        pad = float(self.config.get('pocket_bbox_pad', 4.0))
        min_xyz = min_xyz - pad
        max_xyz = max_xyz + pad

        rng = np.random.default_rng(seed)
        pts = rng.uniform(low=min_xyz, high=max_xyz, size=(n_points, 3)).astype(float)

        min_dists = np.full((n_points,), np.inf, dtype=float)
        for i in range(0, n_points, chunk):
            p = pts[i:i + chunk]
            diff = atoms[None, :, :] - p[:, None, :]
            d2 = np.sum(diff * diff, axis=2)
            min_dists[i:i + len(p)] = np.sqrt(np.min(d2, axis=1))

        mask = (min_dists >= dmin) & (min_dists <= dmax)
        if not np.any(mask):
            return []

        cand_pts = pts[mask]
        cand_scores = min_dists[mask]
        order = np.argsort(-cand_scores)
        if order.size > top_keep:
            order = order[:top_keep]
        cand_pts = cand_pts[order]
        cand_scores = cand_scores[order]

        centers: List[List[float]] = []
        for p in cand_pts:
            centers.append([float(p[0]), float(p[1]), float(p[2])])
        centers = self._cluster_pocket_centers(centers, min_dist=cluster_radius)
        centers = centers[:max_pockets]

        size_cfg = self.config.get('pocket_box_size', self.config.get('box_size_default', [35, 35, 35]))
        if isinstance(size_cfg, (list, tuple)) and len(size_cfg) == 3:
            try:
                size = [float(size_cfg[0]), float(size_cfg[1]), float(size_cfg[2])]
            except Exception:
                size = [35.0, 35.0, 35.0]
        else:
            size = [35.0, 35.0, 35.0]

        pockets: List[Dict] = []
        for idx, c in enumerate(centers, start=1):
            pockets.append({'name': f'surface_random_{idx}', 'center': c, 'size': size})
        return pockets

    def _get_pockets_for_receptor(self, receptor_pdbqt: str) -> List[Dict]:
        key = os.path.abspath(str(receptor_pdbqt))
        cached = self._receptor_pockets_cache.get(key)
        if cached is not None:
            return cached

        pocket_enable = bool(self.config.get('pocket_enable', False))
        if not pocket_enable:
            pockets = [{'name': 'single', 'center': self._get_box_center_for_receptor(receptor_pdbqt), 'size': self._get_box_size_for_receptor(receptor_pdbqt)}]
            self._receptor_pockets_cache[key] = pockets
            return pockets

        algos_cfg = self.config.get('pocket_algorithms', ['cocrystal', 'surface_random', 'receptor_center'])
        if isinstance(algos_cfg, (list, tuple)):
            algos = [str(x).strip().lower() for x in algos_cfg if str(x).strip()]
        else:
            algos = [x.strip().lower() for x in str(algos_cfg).split(',') if x.strip()]

        pockets: List[Dict] = []
        receptor_pdb = self._receptor_pdb_cache.get(key)

        if 'cocrystal' in algos and receptor_pdb:
            inferred = self._infer_box_from_cocrystal_ligand(Path(receptor_pdb))
            if inferred is not None:
                c, s = inferred
                pockets.append({'name': 'cocrystal', 'center': c, 'size': s})

        if 'surface_random' in algos and receptor_pdb:
            pockets.extend(self._predict_pockets_surface_random(receptor_pdb))

        if 'receptor_center' in algos:
            pockets.append({'name': 'receptor_center', 'center': self._get_box_center_for_receptor(receptor_pdbqt), 'size': self._get_box_size_for_receptor(receptor_pdbqt)})

        max_pockets = int(self.config.get('pocket_max_pockets', 5))
        uniq: List[Dict] = []
        seen = set()
        for p in pockets:
            c = p.get('center')
            s = p.get('size')
            if not c or not s:
                continue
            try:
                ck = (round(float(c[0]), 1), round(float(c[1]), 1), round(float(c[2]), 1), round(float(s[0]), 1), round(float(s[1]), 1), round(float(s[2]), 1))
            except Exception:
                continue
            if ck in seen:
                continue
            seen.add(ck)
            uniq.append(p)
            if len(uniq) >= max_pockets:
                break

        if not uniq:
            uniq = [{'name': 'single', 'center': self._get_box_center_for_receptor(receptor_pdbqt), 'size': self._get_box_size_for_receptor(receptor_pdbqt)}]

        self._receptor_pockets_cache[key] = uniq
        return uniq

    def _infer_box_from_cocrystal_ligand(self, receptor_pdb: Path) -> Optional[Tuple[List[float], List[float]]]:
        if (not receptor_pdb.exists()) or receptor_pdb.stat().st_size == 0:
            return None

        keep_hetero = set(str(x).strip().upper() for x in (self.config.get('receptor_cleanup', {}) or {}).get('keep_hetero_resnames', []) if str(x).strip())
        margin = float(self.config.get('box_auto_margin', 10.0))
        min_size = float(self.config.get('box_auto_min_size', 18.0))
        max_size = float(self.config.get('box_auto_max_size', 40.0))

        water = {'HOH', 'WAT', 'H2O', 'DOD'}
        ions = {
            'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'MN', 'FE', 'CU', 'CO', 'NI', 'CD', 'HG', 'BR', 'I', 'F', 'CS', 'SR', 'BA'
        }
        buffers = {
            'SO4', 'PO4', 'EDO', 'GOL', 'PEG', 'PG4', 'DMS', 'MPD', 'ACT', 'FMT', 'TRS', 'MES', 'HEP', 'BME', 'ACE', 'NAG'
        }

        groups: Dict[Tuple[str, str, str], List[Tuple[float, float, float]]] = {}
        try:
            with open(receptor_pdb, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if not line.startswith('HETATM'):
                        continue
                    if len(line) < 54:
                        continue

                    resname = line[17:20].strip().upper()
                    chain = line[21].strip()
                    resseq = line[22:26].strip()

                    if resname in water:
                        continue
                    if resname in ions and resname not in keep_hetero:
                        continue
                    if resname in buffers and resname not in keep_hetero:
                        continue

                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except Exception:
                        continue

                    key = (resname, chain, resseq)
                    groups.setdefault(key, []).append((x, y, z))
        except Exception:
            return None

        if not groups:
            return None

        best_key = max(groups.keys(), key=lambda k: len(groups[k]))
        coords = groups.get(best_key, [])
        if not coords:
            return None

        xs = [c[0] for c in coords]
        ys = [c[1] for c in coords]
        zs = [c[2] for c in coords]

        center = [float(sum(xs) / len(xs)), float(sum(ys) / len(ys)), float(sum(zs) / len(zs))]

        dx = (max(xs) - min(xs)) + 2.0 * margin
        dy = (max(ys) - min(ys)) + 2.0 * margin
        dz = (max(zs) - min(zs)) + 2.0 * margin
        size = [
            float(min(max(dx, min_size), max_size)),
            float(min(max(dy, min_size), max_size)),
            float(min(max(dz, min_size), max_size)),
        ]

        return center, size

    def _looks_like_pdb(self, pdb_path: Path) -> bool:
        """快速判断文件是否像有效 PDB（非严格校验）"""
        if (not pdb_path.exists()) or pdb_path.stat().st_size == 0:
            return False
        try:
            with pdb_path.open('r', encoding='utf-8', errors='ignore') as f:
                for _ in range(300):
                    line = f.readline()
                    if not line:
                        break
                    if line.startswith(('ATOM', 'HETATM')):
                        return True
        except Exception:
            return False
        return False

    def _convert_cif_to_pdb_openbabel(self, cif_file: Path, pdb_file: Path, timeout_sec: int) -> Tuple[bool, str]:
        """使用 OpenBabel 将 CIF/mmCIF 转为 PDB（带超时与最小有效性检查）"""
        if pdb_file.exists():
            pdb_file.unlink()
        cmd = ['obabel', '-icif', str(cif_file), '-opdb', '-O', str(pdb_file), '-d']
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
        except subprocess.TimeoutExpired:
            return False, f'OpenBabel CIF→PDB 超时（>{timeout_sec}s）'
        if result.returncode != 0:
            err = (result.stderr or result.stdout or '')[:200]
            return False, f'OpenBabel CIF→PDB 失败: {err}'
        if not self._looks_like_pdb(pdb_file):
            return False, 'OpenBabel CIF→PDB 产物为空/不含 ATOM'
        return True, ''

    def _convert_cif_to_pdb_gemmi(self, cif_file: Path, pdb_file: Path) -> Tuple[bool, str]:
        """使用 gemmi 将 CIF/mmCIF 转为 PDB（最小有效性检查）"""
        try:
            import gemmi
        except Exception as e:
            return False, f'gemmi 不可用: {e}'
        if pdb_file.exists():
            pdb_file.unlink()
        try:
            st = gemmi.read_structure(str(cif_file))
            st.write_pdb(str(pdb_file))
        except Exception as e:
            return False, f'gemmi CIF→PDB 失败: {e}'
        if not self._looks_like_pdb(pdb_file):
            return False, 'gemmi CIF→PDB 产物为空/不含 ATOM'
        return True, ''

    def _convert_cif_to_pdb_openmm(self, cif_file: Path, pdb_file: Path) -> Tuple[bool, str]:
        """使用 OpenMM 将 CIF/mmCIF 转为 PDB（最小有效性检查）"""
        try:
            from openmm.app import PDBFile, PDBxFile
        except Exception as e:
            return False, f'OpenMM 不可用: {e}'
        if pdb_file.exists():
            pdb_file.unlink()
        try:
            px = PDBxFile(str(cif_file))
            with open(pdb_file, 'w', encoding='utf-8') as f:
                PDBFile.writeFile(px.topology, px.positions, f, keepIds=True)
        except Exception as e:
            return False, f'OpenMM CIF→PDB 失败: {e}'
        if not self._looks_like_pdb(pdb_file):
            return False, 'OpenMM CIF→PDB 产物为空/不含 ATOM'
        return True, ''

    def _convert_cif_to_pdb_with_fallback(self, cif_file: Path, pdb_file: Path) -> Tuple[bool, str]:
        """CIF/mmCIF 转 PDB：按配置执行 OpenBabel(<=timeout)→gemmi→OpenMM 降级"""
        fallback_enable = bool(self.config.get('receptor_cif_to_pdb_fallback_enable', True))
        preferred = str(self.config.get('receptor_cif_to_pdb_preferred_tool', 'openbabel')).strip().lower()
        ob_timeout = int(self.config.get('receptor_cif_to_pdb_openbabel_timeout_sec', 60))

        supported = ['openbabel', 'gemmi', 'openmm']
        if preferred not in supported:
            preferred = 'openbabel'

        if fallback_enable:
            base_order = ['openbabel', 'gemmi', 'openmm']
            order = [preferred] + [t for t in base_order if t != preferred]
        else:
            order = [preferred]

        last_msg = ''
        for tool in order:
            if tool == 'openbabel':
                ok, msg = self._convert_cif_to_pdb_openbabel(cif_file, pdb_file, timeout_sec=ob_timeout)
            elif tool == 'gemmi':
                ok, msg = self._convert_cif_to_pdb_gemmi(cif_file, pdb_file)
            else:
                ok, msg = self._convert_cif_to_pdb_openmm(cif_file, pdb_file)

            if ok:
                return True, ''
            last_msg = msg

        return False, last_msg

    def prepare_receptor_from_structure_file(self, structure_file: str, method: Optional[str] = None) -> Optional[str]:
        """从本地结构文件准备受体（支持 pdb / cif / mmcif）

        Args:
            structure_file: 受体结构文件路径
            method: 'openbabel' 或 'meeko'

        Returns:
            PDBQT 文件路径（失败返回 None）
        """
        if method is None:
            method = self.config.get('receptor_method', 'openbabel')

        in_path = Path(structure_file)
        suffix = in_path.suffix.lower()
        base_name = self._safe_stem(in_path)

        pdb_file = self.receptor_prepared_dir / f"{base_name}.pdb"
        pdbqt_file = self.receptor_prepared_dir / f"{base_name}.pdbqt"

        print(f"正在准备受体: {base_name} (本地文件, 使用{method})...")

        try:
            if suffix in ['.cif', '.mmcif']:
                ok, msg = self._convert_cif_to_pdb_with_fallback(in_path, pdb_file)
                if not ok:
                    print(f"  ✗ CIF转PDB失败: {msg[:200]}")
                    return None
            elif suffix == '.pdb':
                shutil.copyfile(str(in_path), str(pdb_file))
            else:
                print(f"  ✗ 不支持的受体结构格式: {in_path.name}")
                return None

            self._cache_receptor_box_from_pdb(pdb_file, pdbqt_file)

            pdb_file = Path(self.preprocess_receptor_pdb(str(pdb_file)))

            if method.lower() == 'meeko':
                success = self.prepare_receptor_meeko(str(pdb_file), str(pdbqt_file))
                if not success:
                    print("  Meeko失败，回退到OpenBabel...")
                    method = 'openbabel'

            if method.lower() == 'openbabel':
                prepare_receptor4_bin = shutil.which('prepare_receptor4')
                if prepare_receptor4_bin:
                    try:
                        cmd = [prepare_receptor4_bin, '-r', str(pdb_file), '-o', str(pdbqt_file), '-A', 'hydrogens']
                        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                        if result.returncode != 0:
                            raise RuntimeError('prepare_receptor4 failed')
                    except Exception:
                        cmd = f'obabel -ipdb "{pdb_file}" -opdbqt -O "{pdbqt_file}" -p 7.4 -xr'
                        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                        if result.returncode != 0:
                            print(f"  ✗ PDB转PDBQT失败: {result.stderr[:200]}")
                            return None
                else:
                    cmd = f'obabel -ipdb "{pdb_file}" -opdbqt -O "{pdbqt_file}" -p 7.4 -xr'
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    if result.returncode != 0:
                        print(f"  ✗ PDB转PDBQT失败: {result.stderr[:200]}")
                        return None

            cache_key = os.path.abspath(str(pdbqt_file))
            self._receptor_pdb_cache[cache_key] = str(pdb_file)
            self._receptor_pockets_cache.pop(cache_key, None)
            print(f"  ✓ 受体准备完成: {pdbqt_file}")
            return str(pdbqt_file)
        except Exception as e:
            print(f"  ✗ 受体准备失败: {e}")
            return None

    def prepare_ligand_from_file(self, ligand_file: str) -> Tuple[Optional[str], Optional[str]]:
        """从本地配体文件准备配体（支持 pdb/mol/mol2/sdf/pdbqt）

        Args:
            ligand_file: 配体文件路径

        Returns:
            (raw_file, pdbqt_file) 或 (None, None)
        """
        in_path = Path(ligand_file)
        suffix = in_path.suffix.lower()
        base_name = self._safe_stem(in_path)

        raw_copy = self.ligand_raw_dir / f"{base_name}{suffix}"
        pdbqt_out = self.ligand_prepared_dir / f"{base_name}.pdbqt"

        try:
            if suffix == '.pdbqt':
                if in_path.resolve() != pdbqt_out.resolve():
                    shutil.copyfile(str(in_path), str(pdbqt_out))
                return None, str(pdbqt_out)

            if in_path.resolve() != raw_copy.resolve():
                shutil.copyfile(str(in_path), str(raw_copy))
            else:
                raw_copy = in_path

            fmt_map = {
                '.pdb': 'pdb',
                '.mol': 'mol',
                '.mol2': 'mol2',
                '.sdf': 'sdf'
            }
            if suffix not in fmt_map:
                print(f"    ✗ 不支持的配体格式: {in_path.name}")
                return None, None

            add_h = self.config.get('local_ligand_add_hydrogens', True)
            gen3d = self.config.get('local_ligand_gen3d', False)
            split_sdf = self.config.get('local_ligand_split_multimol_sdf', False)

            add_h_flag = '--addh' if add_h else ''
            gen3d_flag = '--gen3d' if gen3d else ''

            if suffix == '.sdf' and split_sdf:
                split_dir = self.ligand_prepared_dir / f"{base_name}_split"
                split_dir.mkdir(parents=True, exist_ok=True)
                split_prefix = split_dir / f"{base_name}_"
                cmd = (
                    f'obabel -i{fmt_map[suffix]} "{raw_copy}" -opdbqt -O "{split_prefix}.pdbqt" '
                    f'-m {add_h_flag} {gen3d_flag} -xn'
                )
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
                if result.returncode != 0:
                    print(f"    OpenBabel转换PDBQT失败: {result.stderr[:200]}")
                    return str(raw_copy), None

                pdbqt_files = sorted(split_dir.glob('*.pdbqt'))
                if not pdbqt_files:
                    return str(raw_copy), None
                shutil.copyfile(str(pdbqt_files[0]), str(pdbqt_out))
                return str(raw_copy), str(pdbqt_out)

            cmd = (
                f'obabel -i{fmt_map[suffix]} "{raw_copy}" -opdbqt -O "{pdbqt_out}" '
                f'{add_h_flag} {gen3d_flag} -xn'
            )
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
            if result.returncode != 0 or not pdbqt_out.exists() or pdbqt_out.stat().st_size == 0:
                if result.stderr:
                    print(f"    OpenBabel转换PDBQT失败: {result.stderr[:200]}")
                return str(raw_copy), None

            return str(raw_copy), str(pdbqt_out)
        except subprocess.TimeoutExpired:
            print("    处理超时（配体可能太大或格式异常）")
            return None, None
        except Exception as e:
            print(f"    本地配体处理失败: {e}")
            return None, None
    
    def download_pdb_structure(self, pdb_id: str, pro_name: str) -> str:
        """
        从RCSB PDB数据库下载蛋白质结构
        
        Args:
            pdb_id: PDB ID
            pro_name: 蛋白质名称
            
        Returns:
            下载的CIF文件路径
        """
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        # 清理文件名中的特殊字符
        safe_name = "".join(c for c in pro_name if c.isalnum() or c in (' ', '-', '_')).strip()
        if not safe_name:
            safe_name = pdb_id
        
        output_file = self.receptor_raw_dir / f"{safe_name}_{pdb_id}.cif"
        
        print(f"正在下载 {pdb_id} ({pro_name})...")
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            
            with open(output_file, 'wb') as f:
                f.write(response.content)
            
            print(f"  ✓ 下载成功: {output_file}")
            return str(output_file)
        except Exception as e:
            print(f"  ✗ 下载失败: {e}")
            return None
    
    def prepare_receptor(self, cif_file: str, method: str = None) -> str:
        """
        准备受体文件（CIF转PDB，添加氢原子，转换为PDBQT）
        
        Args:
            cif_file: CIF文件路径
            method: 'openbabel' 或 'meeko'，默认使用配置文件中的方法
            
        Returns:
            PDBQT文件路径
        """
        if method is None:
            method = self.config.get('receptor_method', 'openbabel')
        
        base_name = Path(cif_file).stem
        pdb_file = self.receptor_prepared_dir / f"{base_name}.pdb"
        pdbqt_file = self.receptor_prepared_dir / f"{base_name}.pdbqt"
        
        print(f"正在准备受体: {base_name} (使用{method})...")
        
        try:
            # 1. CIF转PDB（OpenBabel 超时后自动降级到 gemmi / OpenMM）
            ok, msg = self._convert_cif_to_pdb_with_fallback(Path(cif_file), pdb_file)
            if not ok:
                print(f"  ✗ CIF转PDB失败: {msg[:200]}")
                return None

            self._cache_receptor_box_from_pdb(pdb_file, pdbqt_file)

            # 2. 对受体 PDB 进行预处理：ALTLOC、多构象、水/小分子清理、缺失残基/原子修补
            pdb_file = Path(self.preprocess_receptor_pdb(str(pdb_file)))

            # 3. 根据方法选择处理方式
            if method.lower() == 'meeko':
                # 使用Meeko准备受体
                success = self.prepare_receptor_meeko(str(pdb_file), str(pdbqt_file))
                if not success:
                    print("  Meeko失败，回退到OpenBabel...")
                    method = 'openbabel'
            
            if method.lower() == 'openbabel':
                prepare_receptor4_bin = shutil.which('prepare_receptor4')
                if prepare_receptor4_bin:
                    try:
                        cmd = [prepare_receptor4_bin, '-r', str(pdb_file), '-o', str(pdbqt_file), '-A', 'hydrogens']
                        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                        if result.returncode != 0:
                            raise RuntimeError('prepare_receptor4 failed')
                    except Exception:
                        cmd = f'obabel -ipdb "{pdb_file}" -opdbqt -O "{pdbqt_file}" -p 7.4 -xr'
                        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                        if result.returncode != 0:
                            print(f"  ✗ PDB转PDBQT失败: {result.stderr[:200]}")
                            return None
                else:
                    cmd = f'obabel -ipdb "{pdb_file}" -opdbqt -O "{pdbqt_file}" -p 7.4 -xr'
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    if result.returncode != 0:
                        print(f"  ✗ PDB转PDBQT失败: {result.stderr[:200]}")
                        return None
            
            cache_key = os.path.abspath(str(pdbqt_file))
            self._receptor_pdb_cache[cache_key] = str(pdb_file)
            self._receptor_pockets_cache.pop(cache_key, None)
            print(f"  ✓ 受体准备完成: {pdbqt_file}")
            return str(pdbqt_file)
            
        except Exception as e:
            print(f"  ✗ 受体准备失败: {e}")
            return None

    def preprocess_receptor_pdb(self, pdb_file: str) -> str:
        """受体 PDB 预处理：ALTLOC 处理、溶剂/小分子清理、缺失残基与缺失原子修补

        Args:
            pdb_file: 初始受体 PDB 文件路径

        Returns:
            预处理后的 PDB 文件路径（若预处理失败则退回原始路径）
        """
        cleanup_cfg = self.config.get('receptor_cleanup', {}) or {}
        pdbfixer_cfg = self.config.get('pdbfixer', {}) or {}

        current_pdb = Path(pdb_file)

        try:
            # 1. ALTLOC 多构象处理：按占有率选择最佳构象
            if cleanup_cfg.get('enable', True) and cleanup_cfg.get('enable_altloc', True):
                altloc_pdb = current_pdb.with_name(current_pdb.stem + '_altloc.pdb')
                self._select_altloc_highest_occupancy(str(current_pdb), str(altloc_pdb))
                current_pdb = altloc_pdb

            # 2. 使用 PDBFixer 修补缺失残基/原子与加氢
            if pdbfixer_cfg.get('enable', True):
                fixed_pdb = current_pdb.with_name(current_pdb.stem + '_fixed.pdb')
                fix_success = self._fix_missing_residues_and_atoms(
                    input_pdb=str(current_pdb),
                    output_pdb=str(fixed_pdb),
                    add_missing_residues=pdbfixer_cfg.get('add_missing_residues', True),
                    add_missing_atoms=pdbfixer_cfg.get('add_missing_atoms', True),
                    add_hydrogens=pdbfixer_cfg.get('add_hydrogens', True),
                    ph=float(pdbfixer_cfg.get('ph', 7.4))
                )
                if fix_success:
                    current_pdb = fixed_pdb

            # 3. 水、离子、缓冲剂、共晶配体的有意识删留
            if cleanup_cfg.get('enable', True):
                cleaned_pdb = current_pdb.with_name(current_pdb.stem + '_clean.pdb')
                self._cleanup_solvents_and_ligands(
                    input_pdb=str(current_pdb),
                    output_pdb=str(cleaned_pdb),
                    cleanup_cfg=cleanup_cfg
                )
                current_pdb = cleaned_pdb

            return str(current_pdb)
        except Exception as e:
            print(f"  ⚠ 受体预处理出现问题，将使用原始 PDB：{e}")
            return pdb_file

    def _select_altloc_highest_occupancy(self, input_pdb: str, output_pdb: str) -> None:
        """在 PDB 中按占有率选择 ALTLOC 多构象中占优的构象

        处理逻辑：
        1. 针对每个原子（同一链/残基号/原子名）统计所有 ALTLOC 构象
        2. 若存在多个 ALTLOC（A/B/...），选择占有率最高的那个保留
        3. 写出新的 PDB 文件，并将保留的 ALTLOC 标记统一为空格

        Args:
            input_pdb: 输入 PDB 文件路径
            output_pdb: 输出 PDB 文件路径
        """
        from collections import defaultdict

        with open(input_pdb, 'r') as f:
            lines = f.readlines()

        # 记录每个原子（不含 ALTLOC）的所有候选构象
        groups = defaultdict(list)

        for idx, line in enumerate(lines):
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            # PDB 列字段参考：
            #  1-6  record name
            # 13-16 atom name
            # 17    altLoc
            # 18-20 resName
            # 22    chainID
            # 23-26 resSeq
            # 27    iCode
            # 55-60 occupancy
            if len(line) < 60:
                continue

            altloc = line[16]
            atom_name = line[12:16]
            res_name = line[17:20]
            chain_id = line[21]
            res_seq = line[22:26]
            i_code = line[26]

            key = (chain_id, res_seq, i_code, res_name, atom_name)

            occ_str = line[54:60].strip()
            try:
                occupancy = float(occ_str) if occ_str else 0.0
            except ValueError:
                occupancy = 0.0

            groups[key].append((altloc, occupancy, idx))

        # 为存在 ALTLOC 的原子选择占有率最高的构象
        best_altloc_by_key = {}
        for key, entries in groups.items():
            # 判断是否真的存在多构象（即 altloc 字段中出现非空字符）
            has_altloc = any(e[0].strip() for e in entries)
            if not has_altloc:
                continue

            best_char = None
            best_occ = -1.0
            for altloc, occ, _ in entries:
                alt_char = altloc.strip() or ' '
                if occ > best_occ:
                    best_occ = occ
                    best_char = alt_char

            best_altloc_by_key[key] = best_char

        # 第二次遍历：只写出被选中的 ALTLOC，并将 ALTLOC 字段统一为空格
        output_lines = []
        for idx, line in enumerate(lines):
            if not line.startswith(('ATOM', 'HETATM')):
                output_lines.append(line)
                continue
            if len(line) < 60:
                output_lines.append(line)
                continue

            atom_name = line[12:16]
            res_name = line[17:20]
            chain_id = line[21]
            res_seq = line[22:26]
            i_code = line[26]
            altloc = line[16]

            key = (chain_id, res_seq, i_code, res_name, atom_name)
            chosen_alt = best_altloc_by_key.get(key)

            # 如果该原子没有 ALTLOC 组信息，则原样保留
            if chosen_alt is None:
                output_lines.append(line)
                continue

            this_alt = altloc.strip() or ' '
            if this_alt != chosen_alt:
                # 非最佳占有率构象，跳过
                continue

            # 统一将 ALTLOC 字段改为空格，避免后续程序误识别
            if altloc != ' ':
                line = line[:16] + ' ' + line[17:]

            output_lines.append(line)

        with open(output_pdb, 'w') as f:
            f.writelines(output_lines)

    def _cleanup_solvents_and_ligands(self, input_pdb: str, output_pdb: str, cleanup_cfg: dict) -> None:
        """根据配置对水、离子、缓冲剂和共晶配体进行有意识删留

        Args:
            input_pdb: 输入 PDB 文件路径
            output_pdb: 输出 PDB 文件路径
            cleanup_cfg: 清理相关配置字典
        """
        # 需要强制保留的 HETATM 残基名（例如 HEM/NAG 等重要辅基）
        keep_hetero = set(
            res.strip().upper()
            for res in cleanup_cfg.get('keep_hetero_resnames', [])
            if isinstance(res, str) and res.strip()
        )

        # 常见水分子、无机离子和缓冲剂残基名集合
        water_resnames = {"HOH", "WAT", "H2O"}
        ion_resnames = {
            "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO",
            "NI", "SR", "CS", "RB", "CD"
        }
        buffer_resnames = {
            "HEP", "HEPES", "TRS", "TRIS", "MES", "GOL", "PEG", "MPD",
            "BME", "DMS", "DMSO", "ACE", "SO4", "PO4", "BOG", "NAG", "MAN", "BMA"
        }

        def classify_hetero(resname: str) -> str:
            name = resname.upper()
            if name in water_resnames:
                return 'water'
            if name in ion_resnames:
                return 'ion'
            if name in buffer_resnames:
                return 'buffer'
            return 'ligand'

        remove_waters = cleanup_cfg.get('remove_waters', True)
        remove_ions = cleanup_cfg.get('remove_ions', True)
        remove_buffers = cleanup_cfg.get('remove_buffers', True)
        remove_organics = cleanup_cfg.get('remove_organic_ligands', True)
        enable_cleanup = cleanup_cfg.get('enable', True)

        output_lines = []
        with open(input_pdb, 'r') as f:
            for line in f:
                if not line.startswith(('ATOM', 'HETATM')):
                    # 其它记录（标题、CONECT、TER 等）直接保留
                    output_lines.append(line)
                    continue

                record_name = line[:6].strip().upper()
                resname = line[17:20].strip().upper()

                if record_name == 'ATOM':
                    # 标准氨基酸/核酸原子，全部保留
                    output_lines.append(line)
                    continue

                # 下面只处理 HETATM
                if resname in keep_hetero:
                    output_lines.append(line)
                    continue

                if not enable_cleanup:
                    output_lines.append(line)
                    continue

                category = classify_hetero(resname)
                remove_flag = False

                if category == 'water':
                    remove_flag = remove_waters
                elif category == 'ion':
                    remove_flag = remove_ions
                elif category == 'buffer':
                    remove_flag = remove_buffers
                else:
                    remove_flag = remove_organics

                if remove_flag:
                    # 跳过需要删除的溶剂/小分子
                    continue

                output_lines.append(line)

        with open(output_pdb, 'w') as f:
            f.writelines(output_lines)

    def _fix_missing_residues_and_atoms(
        self,
        input_pdb: str,
        output_pdb: str,
        add_missing_residues: bool = True,
        add_missing_atoms: bool = True,
        add_hydrogens: bool = True,
        ph: float = 7.4
    ) -> bool:
        """使用 PDBFixer 修补缺失残基/缺失原子并根据 pH 加氢

        Args:
            input_pdb: 输入 PDB 文件路径
            output_pdb: 输出 PDB 文件路径
            add_missing_residues: 是否补全缺失残基/loop
            add_missing_atoms: 是否补全缺失原子
            add_hydrogens: 是否添加氢原子
            ph: 加氢时使用的 pH 值

        Returns:
            是否修补成功（失败时调用方可以回退到原始结构）
        """
        try:
            from pdbfixer import PDBFixer
            try:
                from openmm.app import PDBFile
            except ImportError:
                from simtk.openmm.app import PDBFile
        except ImportError as e:
            print(f"  ⚠ 未找到 PDBFixer/OpenMM，跳过缺失残基与缺失原子修补: {e}")
            return False

        try:
            self._sanitize_pdb_numeric_fields_for_pdbfixer(input_pdb)
            fixer = PDBFixer(filename=input_pdb)

            if add_missing_residues:
                fixer.findMissingResidues()

            if add_missing_atoms:
                fixer.findMissingAtoms()
                fixer.addMissingAtoms()

            if add_hydrogens:
                fixer.addMissingHydrogens(pH=ph)

            with open(output_pdb, 'w') as out:
                PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)

            return True
        except Exception as e:
            print(f"  ⚠ PDBFixer 修补失败，将使用原始结构: {e}")
            return False

    def _sanitize_pdb_numeric_fields_for_pdbfixer(self, pdb_file: str) -> None:
        """将 PDB 中占有率/温度因子字段的异常值修正为数值，避免 PDBFixer 解析失败"""
        pdb_path = Path(pdb_file)
        if (not pdb_path.exists()) or pdb_path.stat().st_size == 0:
            return

        def _fmt_float_field(raw: str, default: float) -> str:
            try:
                val = float(str(raw).strip())
            except Exception:
                val = float(default)
            return f"{val:6.2f}"

        changed = False
        out_lines: List[str] = []
        with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    base = line.rstrip('\n')
                    if len(base) < 66:
                        base = base.ljust(66)
                    occ = base[54:60]
                    bfac = base[60:66]
                    new_occ = _fmt_float_field(occ, 1.0)
                    new_bfac = _fmt_float_field(bfac, 0.0)
                    if new_occ != occ or new_bfac != bfac:
                        changed = True
                        base = base[:54] + new_occ + new_bfac + base[66:]
                    out_lines.append(base + '\n')
                else:
                    out_lines.append(line)

        if changed:
            with open(pdb_path, 'w', encoding='utf-8') as out:
                out.writelines(out_lines)
    
    def prepare_receptor_meeko(self, pdb_file: str, pdbqt_file: str) -> bool:
        """
        使用Meeko准备受体
        
        Args:
            pdb_file: PDB文件路径
            pdbqt_file: 输出PDBQT文件路径
            
        Returns:
            是否成功
        """
        try:
            from meeko import PDBQTWriterLegacy, MoleculePreparation
            from rdkit import Chem
            from rdkit.Chem import AllChem
            
            # 读取PDB文件
            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                return False
            
            # 添加氢原子
            mol = Chem.AddHs(mol, addCoords=True)
            
            # 使用Meeko准备
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            
            # 写入PDBQT
            for setup in mol_setups:
                pdbqt_string = PDBQTWriterLegacy.write_string(setup)
                with open(pdbqt_file, 'w') as f:
                    f.write(pdbqt_string[0] if isinstance(pdbqt_string, tuple) else pdbqt_string)
                break
            
            return os.path.exists(pdbqt_file)
            
        except Exception as e:
            print(f"    Meeko处理受体失败: {e}")
            return False
    
    def smiles_to_3d_openbabel(self, smiles: str, ligand_name: str) -> tuple:
        """
        使用OpenBabel将SMILES转换为3D结构并准备配体
        
        Args:
            smiles: SMILES字符串
            ligand_name: 配体名称（安全文件名）
            
        Returns:
            (pdb_file, pdbqt_file) 或 (None, None)
        """
        try:
            smiles = str(smiles).strip().replace('\r', '').replace('\n', '')

            # 分别保存到raw和prepared目录
            pdb_file = self.ligand_raw_dir / f"{ligand_name}.pdb"
            pdbqt_file = self.ligand_prepared_dir / f"{ligand_name}.pdbqt"
            
            pdb_file = str(pdb_file)
            pdbqt_file = str(pdbqt_file)
            
            # 删除旧文件避免误判
            for f in [pdb_file, pdbqt_file]:
                if os.path.exists(f):
                    os.remove(f)

            large_smiles = len(smiles) >= 120
            gen3d_timeout = int(self.config.get('ligand_openbabel_gen3d_timeout_sec', 300 if large_smiles else 120))
            pdbqt_timeout = int(self.config.get('ligand_openbabel_pdbqt_timeout_sec', 60))

            gen3d_attempts = [
                ['obabel', f'-:{smiles}', '-opdbqt', '-O', pdbqt_file, '--gen3d', '--addh', '-xn'],
                ['obabel', f'-:{smiles}', '-opdbqt', '-O', pdbqt_file, '--gen3d', '--addh', '--minimize', '--steps', '200', '--ff', 'MMFF94', '-xn'],
                ['obabel', f'-:{smiles}', '-opdbqt', '-O', pdbqt_file, '--gen3d', '--addh', '--minimize', '--steps', '200', '--ff', 'UFF', '-xn'],
            ]

            last_err = ''
            gen3d_ok = False
            for cmd in gen3d_attempts:
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=gen3d_timeout)
                except subprocess.TimeoutExpired:
                    last_err = f"gen3d timeout after {gen3d_timeout}s"
                    continue

                if result.returncode == 0 and os.path.exists(pdbqt_file) and os.path.getsize(pdbqt_file) > 0:
                    gen3d_ok = True
                    break
                last_err = (result.stderr or result.stdout or '')[:400]

            if not gen3d_ok:
                print("    OpenBabel生成3D结构失败")
                if last_err:
                    print(f"    错误信息: {last_err}")
                return None, None
            
            cmd = ['obabel', '-ipdbqt', pdbqt_file, '-opdb', '-O', pdb_file]
            subprocess.run(cmd, capture_output=True, text=True, timeout=pdbqt_timeout)

            return pdb_file, pdbqt_file
            
        except subprocess.TimeoutExpired:
            print(f"    处理超时（分子可能太复杂）")
            return None, None

    def smiles_to_3d_rdkit(self, smiles: str, ligand_name: str) -> tuple:
        """使用RDKit生成3D构象，再用OpenBabel转为PDBQT"""
        try:
            smiles = str(smiles).strip().replace('\r', '').replace('\n', '')
            if not smiles:
                print("    SMILES为空")
                return None, None

            from rdkit import Chem
            from rdkit.Chem import AllChem

            pdb_file = self.ligand_raw_dir / f"{ligand_name}.pdb"
            sdf_file = self.ligand_raw_dir / f"{ligand_name}.sdf"
            pdbqt_file = self.ligand_prepared_dir / f"{ligand_name}.pdbqt"

            pdb_file = str(pdb_file)
            sdf_file = str(sdf_file)
            pdbqt_file = str(pdbqt_file)

            for f in [pdb_file, sdf_file, pdbqt_file]:
                if os.path.exists(f):
                    os.remove(f)

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print("    无法解析SMILES")
                return None, None

            mol = Chem.AddHs(mol)

            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            embed_res = AllChem.EmbedMolecule(mol, params)
            if embed_res != 0:
                print("    生成3D坐标失败")
                return None, None

            if AllChem.MMFFHasAllMoleculeParams(mol):
                AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)
            else:
                AllChem.UFFOptimizeMolecule(mol, maxIters=2000)

            Chem.MolToPDBFile(mol, pdb_file)
            w = Chem.SDWriter(sdf_file)
            w.write(mol)
            w.close()

            timeout_sec = int(self.config.get('ligand_rdkit_to_pdbqt_timeout_sec', 120))
            cmd = ['obabel', '-isdf', sdf_file, '-opdbqt', '-O', pdbqt_file, '--addh', '-xn']
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
            if result.returncode != 0 or (not os.path.exists(pdbqt_file)) or os.path.getsize(pdbqt_file) == 0:
                err = (result.stderr or result.stdout or '')[:400]
                if err:
                    print(f"    RDKit->PDBQT失败: {err}")
                return pdb_file if os.path.exists(pdb_file) else None, None

            return pdb_file if os.path.exists(pdb_file) else None, pdbqt_file

        except ImportError as e:
            print(f"    缺少RDKit相关依赖: {e}")
            return None, None
        except subprocess.TimeoutExpired:
            print("    RDKit->PDBQT超时（分子可能太复杂）")
            return None, None
        except Exception as e:
            print(f"    RDKit处理失败: {e}")
            return None, None
        except Exception as e:
            print(f"    OpenBabel处理失败: {e}")
            return None, None
    
    def smiles_to_3d_meeko(self, smiles: str, output_file: str) -> bool:
        """
        使用RDKit + Meeko将SMILES转换为3D结构并准备配体
        
        Args:
            smiles: SMILES字符串
            output_file: 输出文件路径（不含扩展名）
            
        Returns:
            是否成功
        """
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from meeko import MoleculePreparation
            from meeko import PDBQTWriterLegacy
            
            # 1. SMILES -> RDKit分子
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"    无法解析SMILES")
                return False
            
            # 2. 添加氢原子
            mol = Chem.AddHs(mol)
            
            # 3. 生成3D坐标
            result = AllChem.EmbedMolecule(mol, randomSeed=42)
            if result != 0:
                print(f"    生成3D坐标失败")
                return False
            
            # 4. 能量最小化
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            
            # 5. 使用Meeko准备配体
            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)
            
            # 6. 写入PDBQT文件
            pdbqt_file = f"{output_file}.pdbqt"
            for setup in mol_setups:
                pdbqt_string = PDBQTWriterLegacy.write_string(setup)
                with open(pdbqt_file, 'w') as f:
                    f.write(pdbqt_string)
                break  # 只使用第一个构象
            
            return os.path.exists(pdbqt_file)
            
        except ImportError as e:
            print(f"    缺少必要的库: {e}")
            print("    请安装: pip install rdkit meeko")
            return False

        except Exception as e:
            print(f"    Meeko处理失败: {e}")
            return False

    def prepare_ligand_entries(self, smiles: str, ligand_name: str) -> List[Dict]:
        safe_name = "".join(c for c in str(ligand_name) if c.isalnum() or c in (' ', '-', '_')).strip()
        if not safe_name:
            safe_name = "ligand"

        enum_enable = bool(self.config.get('ligand_enumeration_enable', False))
        if not enum_enable:
            pdb_file, pdbqt_file = self.prepare_ligand(smiles, ligand_name)
            if pdbqt_file:
                return [{
                    'name': safe_name,
                    'parent_name': safe_name,
                    'smiles': str(smiles).strip(),
                    'pdb': pdb_file,
                    'pdbqt': pdbqt_file
                }]
            return []

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
        except Exception:
            pdb_file, pdbqt_file = self.prepare_ligand(smiles, ligand_name)
            if pdbqt_file:
                return [{
                    'name': safe_name,
                    'parent_name': safe_name,
                    'smiles': str(smiles).strip(),
                    'pdb': pdb_file,
                    'pdbqt': pdbqt_file
                }]
            return []

        smiles_clean = str(smiles).strip().replace('\r', '').replace('\n', '')
        if not smiles_clean:
            return []

        mol0 = Chem.MolFromSmiles(smiles_clean)
        if mol0 is None:
            return []

        max_tautomers = int(self.config.get('ligand_enum_max_tautomers', 8))
        max_conformers = int(self.config.get('ligand_enum_max_conformers', 10))
        max_variants = int(self.config.get('ligand_enum_max_variants', 24))
        prune_rms = float(self.config.get('ligand_enum_prune_rms_thresh', 0.5))
        energy_window = float(self.config.get('ligand_enum_energy_window_kcal', 5.0))
        random_seed = int(self.config.get('ligand_enum_random_seed', 42))
        embed_confs = int(self.config.get('ligand_enum_embed_conformers', max(max_conformers * 2, 20)))
        timeout_sec = int(self.config.get('ligand_rdkit_to_pdbqt_timeout_sec', 120))

        lig_ph = self.config.get('ligand_ph', None)
        ph_flag: List[str] = []
        try:
            if lig_ph is not None and str(lig_ph).strip() != '':
                ph_flag = ['-p', str(float(lig_ph))]
        except Exception:
            ph_flag = []

        tautomers: List['Chem.Mol'] = []
        try:
            from rdkit.Chem.MolStandardize import rdMolStandardize
            enumerator = rdMolStandardize.TautomerEnumerator()
            enumerated = enumerator.Enumerate(mol0)
            tautomers = list(enumerated)
        except Exception:
            tautomers = [mol0]

        seen = set()
        uniq_tautomers: List['Chem.Mol'] = []
        for m in tautomers:
            try:
                smi = Chem.MolToSmiles(m, isomericSmiles=True)
            except Exception:
                continue
            if smi in seen:
                continue
            seen.add(smi)
            uniq_tautomers.append(m)
            if len(uniq_tautomers) >= max_tautomers:
                break

        prepared: List[Dict] = []
        variant_no = 0
        for t_idx, t_mol in enumerate(uniq_tautomers, start=1):
            if variant_no >= max_variants:
                break

            m = Chem.AddHs(Chem.Mol(t_mol))

            params = AllChem.ETKDGv3()
            params.randomSeed = random_seed
            params.pruneRmsThresh = prune_rms
            try:
                conf_ids = list(AllChem.EmbedMultipleConfs(m, numConfs=embed_confs, params=params))
            except Exception:
                conf_ids = []

            if not conf_ids:
                try:
                    ok = AllChem.EmbedMolecule(m, params)
                    conf_ids = [0] if ok == 0 else []
                except Exception:
                    conf_ids = []

            if not conf_ids:
                continue

            energies: List[Tuple[int, float]] = []
            try:
                res = AllChem.MMFFOptimizeMoleculeConfs(m, numThreads=0, maxIters=500)
                for cid, (status, e) in zip(conf_ids, res):
                    if status == 0 and e is not None:
                        energies.append((int(cid), float(e)))
            except Exception:
                energies = []

            if not energies:
                try:
                    res = AllChem.UFFOptimizeMoleculeConfs(m, numThreads=0, maxIters=1000)
                    for cid, (status, e) in zip(conf_ids, res):
                        if status == 0 and e is not None:
                            energies.append((int(cid), float(e)))
                except Exception:
                    energies = []

            if energies:
                energies.sort(key=lambda x: x[1])
                best_e = energies[0][1]
                selected = [x for x in energies if x[1] <= best_e + energy_window]
                selected = selected[:max_conformers]
                selected_conf_ids = [cid for cid, _ in selected]
            else:
                selected_conf_ids = [int(x) for x in conf_ids[:max_conformers]]

            for c_idx, conf_id in enumerate(selected_conf_ids, start=1):
                if variant_no >= max_variants:
                    break

                variant_no += 1
                variant_suffix = f"__t{t_idx:02d}_c{c_idx:02d}"
                variant_name = f"{safe_name}{variant_suffix}"

                sdf_file = self.ligand_raw_dir / f"{variant_name}.sdf"
                pdb_file = self.ligand_raw_dir / f"{variant_name}.pdb"
                pdbqt_file = self.ligand_prepared_dir / f"{variant_name}.pdbqt"

                for fp in [sdf_file, pdb_file, pdbqt_file]:
                    try:
                        if fp.exists():
                            fp.unlink()
                    except Exception:
                        pass

                try:
                    writer = Chem.SDWriter(str(sdf_file))
                    writer.write(m, confId=int(conf_id))
                    writer.close()
                except Exception:
                    continue

                try:
                    Chem.MolToPDBFile(m, str(pdb_file), confId=int(conf_id))
                except Exception:
                    pdb_file = None

                cmd = ['obabel', '-isdf', str(sdf_file), '-opdbqt', '-O', str(pdbqt_file), '--addh', '-xn']
                if ph_flag:
                    cmd = ['obabel', '-isdf', str(sdf_file), '-opdbqt', '-O', str(pdbqt_file)] + ph_flag + ['--addh', '-xn']
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_sec)
                    if result.returncode != 0 or (not pdbqt_file.exists()) or pdbqt_file.stat().st_size == 0:
                        continue
                except Exception:
                    continue

                prepared.append({
                    'name': variant_name,
                    'parent_name': safe_name,
                    'smiles': smiles_clean,
                    'pdb': str(pdb_file) if pdb_file and Path(pdb_file).exists() else None,
                    'pdbqt': str(pdbqt_file)
                })

        if prepared:
            print(f"正在准备配体: {ligand_name} (枚举模式)...")
            print(f"  ✓ 生成 {len(prepared)} 个配体变体用于对接")
            return prepared

        pdb_file, pdbqt_file = self.prepare_ligand(smiles, ligand_name)
        if pdbqt_file:
            return [{
                'name': safe_name,
                'parent_name': safe_name,
                'smiles': smiles_clean,
                'pdb': pdb_file,
                'pdbqt': pdbqt_file
            }]
        return []
    
    def prepare_ligand(self, smiles: str, ligand_name: str) -> tuple:
        """
        准备配体文件（支持OpenBabel / RDKit / Meeko）
        
        Args:
            smiles: SMILES字符串
            ligand_name: 配体名称
            
        Returns:
            (pdb_file, pdbqt_file) 或 (None, None)
        """
        # 清理文件名
        safe_name = "".join(c for c in ligand_name if c.isalnum() or c in (' ', '-', '_')).strip()
        if not safe_name:
            safe_name = "ligand"
        
        ligand_method = str(self.config.get('ligand_method', 'openbabel')).strip().lower()

        if ligand_method == 'meeko':
            print(f"正在准备配体: {ligand_name} (使用Meeko)...")
            output_base = str(self.ligand_prepared_dir / safe_name)
            ok = self.smiles_to_3d_meeko(smiles, output_base)
            pdb_file = str(self.ligand_raw_dir / f"{safe_name}.pdb")
            pdbqt_file = f"{output_base}.pdbqt"
            if not ok or not os.path.exists(pdbqt_file):
                print(f"  ✗ 配体准备失败")
                return None, None
            print(f"  ✓ 配体准备完成")
            print(f"    Prepared: {pdbqt_file}")
            return pdb_file if os.path.exists(pdb_file) else None, pdbqt_file

        if ligand_method == 'rdkit':
            print(f"正在准备配体: {ligand_name} (使用RDKit)...")
            pdb_file, pdbqt_file = self.smiles_to_3d_rdkit(smiles, safe_name)
            if pdbqt_file and os.path.exists(pdbqt_file):
                print(f"  ✓ 配体准备完成")
                if pdb_file:
                    print(f"    Raw: {pdb_file}")
                print(f"    Prepared: {pdbqt_file}")
                return pdb_file, pdbqt_file
            print(f"  ✗ 配体准备失败")
            return None, None

        print(f"正在准备配体: {ligand_name} (使用OpenBabel)...")

        pdb_file, pdbqt_file = self.smiles_to_3d_openbabel(smiles, safe_name)

        if not pdbqt_file:
            fallback_cfg = self.config.get('ligand_method_fallback', 'meeko')
            if isinstance(fallback_cfg, (list, tuple)):
                fallback_methods = [str(x).strip().lower() for x in fallback_cfg if str(x).strip()]
            else:
                fallback_methods = [m.strip().lower() for m in str(fallback_cfg).split(',') if m.strip()]

            for method in fallback_methods:
                if method == 'rdkit':
                    print(f"  ⚠ OpenBabel失败，尝试使用RDKit回退...")
                    pdb_file2, pdbqt2 = self.smiles_to_3d_rdkit(smiles, safe_name)
                    if pdbqt2 and os.path.exists(pdbqt2):
                        print(f"  ✓ 配体准备完成")
                        print(f"    Prepared: {pdbqt2}")
                        return pdb_file2, pdbqt2
                elif method == 'meeko':
                    print(f"  ⚠ OpenBabel失败，尝试使用Meeko回退...")
                    output_base = str(self.ligand_prepared_dir / safe_name)
                    ok = self.smiles_to_3d_meeko(smiles, output_base)
                    pdbqt2 = f"{output_base}.pdbqt"
                    if ok and os.path.exists(pdbqt2):
                        print(f"  ✓ 配体准备完成")
                        print(f"    Prepared: {pdbqt2}")
                        pdb_file2 = str(self.ligand_raw_dir / f"{safe_name}.pdb")
                        return pdb_file2 if os.path.exists(pdb_file2) else None, pdbqt2

        if pdbqt_file and os.path.exists(pdbqt_file):
            print(f"  ✓ 配体准备完成")
            print(f"    Raw: {pdb_file}")
            print(f"    Prepared: {pdbqt_file}")
            return pdb_file, pdbqt_file
        else:
            print(f"  ✗ 配体准备失败")
            return None, None

    def _compute_pdbqt_bounds(self, pdbqt_file: str) -> Optional[Tuple[float, float, float, float, float, float]]:
        """计算PDBQT中所有原子坐标的包围盒"""
        try:
            min_x = min_y = min_z = float('inf')
            max_x = max_y = max_z = float('-inf')
            found = False

            with open(pdbqt_file, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    if not line.startswith(('ATOM', 'HETATM')):
                        continue

                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except Exception:
                        parts = line.split()
                        if len(parts) < 6:
                            continue
                        try:
                            x = float(parts[5])
                            y = float(parts[6])
                            z = float(parts[7])
                        except Exception:
                            continue

                    found = True
                    min_x = min(min_x, x)
                    min_y = min(min_y, y)
                    min_z = min(min_z, z)
                    max_x = max(max_x, x)
                    max_y = max(max_y, y)
                    max_z = max(max_z, z)

            if not found:
                return None
            return min_x, max_x, min_y, max_y, min_z, max_z
        except Exception:
            return None

    def _get_box_center_for_receptor(self, receptor_pdbqt: str) -> List[float]:
        """获取指定受体的对接盒中心坐标（支持配置为auto自动计算）"""
        center_cfg = self.config.get('box_center', [0, 0, 0])
        if isinstance(center_cfg, (list, tuple)) and len(center_cfg) == 3:
            try:
                return [float(center_cfg[0]), float(center_cfg[1]), float(center_cfg[2])]
            except Exception:
                pass

        center_mode = str(center_cfg).strip().lower() if isinstance(center_cfg, str) else ''

        if isinstance(center_cfg, str) and center_mode not in ('auto', 'auto_ligand', 'ligand', 'cocrystal'):
            try:
                parts = [p.strip() for p in center_cfg.split(',')]
                if len(parts) == 3:
                    return [float(parts[0]), float(parts[1]), float(parts[2])]
            except Exception:
                pass

        cache_key = os.path.abspath(receptor_pdbqt)
        cached = self._receptor_center_cache.get(cache_key)
        if cached is not None:
            return cached

        bounds = self._compute_pdbqt_bounds(receptor_pdbqt)
        if bounds is None:
            center = [0.0, 0.0, 0.0]
        else:
            min_x, max_x, min_y, max_y, min_z, max_z = bounds
            center = [
                (min_x + max_x) / 2.0,
                (min_y + max_y) / 2.0,
                (min_z + max_z) / 2.0,
            ]

        self._receptor_center_cache[cache_key] = center
        return center

    def _get_box_size_for_receptor(self, receptor_pdbqt: str) -> List[float]:
        size_cfg = self.config.get('box_size', [20, 20, 20])
        if isinstance(size_cfg, (list, tuple)) and len(size_cfg) == 3:
            try:
                return [float(size_cfg[0]), float(size_cfg[1]), float(size_cfg[2])]
            except Exception:
                pass

        size_mode = str(size_cfg).strip().lower() if isinstance(size_cfg, str) else ''
        if isinstance(size_cfg, str) and size_mode in ('auto', 'auto_ligand', 'ligand', 'cocrystal'):
            cache_key = os.path.abspath(receptor_pdbqt)
            cached = self._receptor_size_cache.get(cache_key)
            if cached is not None:
                return cached
            default_size = self.config.get('box_size_default', [35, 35, 35])
            if isinstance(default_size, (list, tuple)) and len(default_size) == 3:
                try:
                    return [float(default_size[0]), float(default_size[1]), float(default_size[2])]
                except Exception:
                    pass
            return [35.0, 35.0, 35.0]

        if isinstance(size_cfg, str):
            try:
                parts = [p.strip() for p in size_cfg.split(',')]
                if len(parts) == 3:
                    return [float(parts[0]), float(parts[1]), float(parts[2])]
            except Exception:
                pass

        return [35.0, 35.0, 35.0]
    
    def run_vina_docking(
        self,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        receptor_name: str,
        ligand_name: str,
        receptor_id: str = '',
        center: Optional[List[float]] = None,
        size: Optional[List[float]] = None,
        pocket_name: str = '',
        pocket_index: Optional[int] = None,
    ) -> Tuple[str, float]:
        """
        运行AutoDock Vina对接
        
        Args:
            receptor_pdbqt: 受体PDBQT文件
            ligand_pdbqt: 配体PDBQT文件
            receptor_name: 受体名称
            ligand_name: 配体名称
            
        Returns:
            (输出目录路径, 最佳对接得分)
        """
        receptor_storage = receptor_name
        if str(receptor_id).strip() != '':
            receptor_storage = f"{receptor_name}_{receptor_id}"

        pair_dir = self.results_dir / f"{receptor_storage}_vs_{ligand_name}"
        if pocket_name or (pocket_index is not None):
            safe_pocket = "".join(c for c in str(pocket_name) if c.isalnum() or c in ('-', '_')).strip('_')
            if not safe_pocket:
                safe_pocket = 'pocket'
            idx_tag = ''
            if pocket_index is not None:
                try:
                    idx_tag = f"{int(pocket_index) + 1:02d}_"
                except Exception:
                    idx_tag = ''
            docking_dir = pair_dir / f"pocket_{idx_tag}{safe_pocket}"
        else:
            docking_dir = pair_dir

        docking_dir.mkdir(parents=True, exist_ok=True)
        
        # 定义输出文件
        output_poses = docking_dir / "docked_poses.pdbqt"
        output_complex = docking_dir / "best_complex.pdb"
        log_file = docking_dir / "docking_log.txt"
        
        if center is None:
            center = self._get_box_center_for_receptor(receptor_pdbqt)
        if size is None:
            size = self._get_box_size_for_receptor(receptor_pdbqt)
        
        cmd = (
            f'vina --receptor "{receptor_pdbqt}" '
            f'--ligand "{ligand_pdbqt}" '
            f'--out "{output_poses}" '
            f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} "
            f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} "
            f"--exhaustiveness {self.config['exhaustiveness']} "
            f"--num_modes {self.config['num_modes']} "
            f"--energy_range {self.config['energy_range']} "
            f"--cpu {self.config.get('cpu', 4)}"
        )
        
        try:
            # 重定向输出到日志文件
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                   text=True, timeout=300)
            
            # 保存完整日志
            with open(log_file, 'w', encoding='utf-8') as f:
                f.write("=" * 60 + "\n")
                f.write("AutoDock Vina Docking Log\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"Receptor: {receptor_pdbqt}\n")
                f.write(f"Ligand: {ligand_pdbqt}\n")
                f.write(f"Output: {output_poses}\n\n")
                f.write("=" * 60 + "\n")
                f.write("Vina Output:\n")
                f.write("=" * 60 + "\n")
                f.write(result.stdout)
                if result.stderr:
                    f.write("\n" + "=" * 60 + "\n")
                    f.write("Stderr:\n")
                    f.write("=" * 60 + "\n")
                    f.write(result.stderr)
            
            # 检查是否成功
            if result.returncode != 0:
                print(f"  ✗ Vina执行失败: {result.stderr[:200]}")
                return None, None
            
            # 解析对接得分
            score = self.parse_vina_output(result.stdout)
            
            # 生成最佳复合物PDB文件（合并受体和最佳配体姿态）
            try:
                self.create_complex_pdb(receptor_pdbqt, output_poses, output_complex)
            except Exception as e:
                print(f"  ⚠ 生成复合物文件失败: {e}")
            
            # 保存对接摘要
            summary_file = docking_dir / "summary.txt"
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write(f"Receptor: {receptor_name}\n")
                if str(receptor_id).strip() != '':
                    f.write(f"Receptor_ID: {receptor_id}\n")
                f.write(f"Ligand: {ligand_name}\n")
                if pocket_name or (pocket_index is not None):
                    f.write(f"Pocket: {pocket_name}\n")
                    if pocket_index is not None:
                        f.write(f"Pocket_Index: {pocket_index + 1}\n")
                    f.write(f"Box_Center: {center[0]}, {center[1]}, {center[2]}\n")
                    f.write(f"Box_Size: {size[0]}, {size[1]}, {size[2]}\n")
                f.write(f"Best Affinity: {score} kcal/mol\n")
                f.write(f"\nFiles Generated:\n")
                f.write(f"  - docked_poses.pdbqt: All docking poses\n")
                f.write(f"  - best_complex.pdb: Best pose with receptor\n")
                f.write(f"  - docking_log.txt: Complete docking log\n")
            
            return str(docking_dir), score
            
        except subprocess.TimeoutExpired:
            print(f"  ✗ 对接超时")
            return None, None
        except Exception as e:
            print(f"  ✗ 对接失败: {e}")
            return None, None
    
    def create_complex_pdb(self, receptor_pdbqt: str, ligand_poses_pdbqt: str, output_pdb: str):
        """
        创建受体-配体复合物PDB文件（最佳姿态）
        
        Args:
            receptor_pdbqt: 受体PDBQT文件
            ligand_poses_pdbqt: 配体姿态PDBQT文件
            output_pdb: 输出复合物PDB文件
        """
        try:
            # 提取最佳配体姿态（第一个MODEL）
            best_ligand_pdbqt = str(output_pdb).replace('.pdb', '_ligand_best.pdbqt')
            with open(ligand_poses_pdbqt, 'r') as f:
                lines = f.readlines()
            
            # 提取第一个模型
            in_model = False
            model_lines = []
            for line in lines:
                if line.startswith('MODEL 1'):
                    in_model = True
                    continue
                if line.startswith('ENDMDL'):
                    break
                if in_model:
                    model_lines.append(line)
            
            with open(best_ligand_pdbqt, 'w') as f:
                f.writelines(model_lines)
            
            # 转换为PDB
            receptor_pdb = str(output_pdb).replace('.pdb', '_receptor.pdb')
            ligand_pdb = str(output_pdb).replace('.pdb', '_ligand_best.pdb')
            
            # 受体PDBQT转PDB
            cmd = f'obabel -ipdbqt "{receptor_pdbqt}" -opdb -O "{receptor_pdb}"'
            subprocess.run(cmd, shell=True, capture_output=True, timeout=30)
            
            # 配体PDBQT转PDB
            cmd = f'obabel -ipdbqt "{best_ligand_pdbqt}" -opdb -O "{ligand_pdb}"'
            subprocess.run(cmd, shell=True, capture_output=True, timeout=30)
            
            # 合并受体和配体
            if os.path.exists(receptor_pdb) and os.path.exists(ligand_pdb):
                with open(output_pdb, 'w') as out:
                    # 写入受体
                    with open(receptor_pdb, 'r') as f:
                        for line in f:
                            if not line.startswith('END'):
                                out.write(line)
                    
                    # 写入配体
                    out.write('TER\n')
                    with open(ligand_pdb, 'r') as f:
                        for line in f:
                            if line.startswith(('ATOM', 'HETATM')):
                                out.write(line)
                    
                    out.write('END\n')
        
        except Exception as e:
            raise Exception(f"创建复合物失败: {e}")
    
    def parse_vina_output(self, output_text: str) -> float:
        """解析Vina输出获取最佳对接得分"""
        try:
            lines = output_text.split('\n')
            for i, line in enumerate(lines):
                # 查找结果表格的标题行（包含mode和affinity）
                if 'mode' in line.lower() and 'affinity' in line.lower():
                    # 跳过标题行和分隔线，查找第一行数据
                    for j in range(i + 1, min(i + 10, len(lines))):
                        score_line = lines[j].strip()
                        # 跳过空行、分隔线和表头续行
                        if not score_line or '----' in score_line or score_line.startswith('|'):
                            continue
                        # 找到数据行
                        parts = score_line.split()
                        # 第一列是mode编号，第二列是affinity得分
                        if len(parts) >= 2:
                            try:
                                # 验证第一列是数字（mode编号）
                                int(parts[0])
                                # 返回第二列作为得分
                                return float(parts[1])
                            except ValueError:
                                continue
            return None
        except Exception as e:
            print(f"  解析得分失败: {e}")
            return None

    def _run_docking_matrix(self, receptors: List[Dict], ligands: List[Dict]) -> pd.DataFrame:
        """对已准备好的受体与配体进行批量对接并保存结果

        Args:
            receptors: 受体列表，每个元素至少包含 name/id/pdbqt
            ligands: 配体列表，每个元素至少包含 name/pdbqt

        Returns:
            对接结果 DataFrame
        """
        print("步骤 3/3: 执行批量对接")
        print("-" * 60)
        results = []
        total_dockings = len(receptors) * len(ligands)
        current = 0

        for receptor in receptors:
            pockets = self._get_pockets_for_receptor(receptor['pdbqt'])
            for ligand in ligands:
                current += 1
                parent_name = ligand.get('parent_name', ligand['name'])
                ligand_run_name = ligand['name']
                if parent_name != ligand_run_name:
                    print(f"\n[{current}/{total_dockings}] 对接: {receptor['name']} vs {parent_name} ({ligand_run_name})")
                else:
                    print(f"\n[{current}/{total_dockings}] 对接: {receptor['name']} vs {ligand_run_name}")

                best_score = None
                best_output_dir = None
                best_pocket = None
                best_pocket_index = None

                if pockets and len(pockets) > 1:
                    print(f"  口袋数量: {len(pockets)}")

                for pi, pocket in enumerate(pockets or []):
                    pocket_name = str(pocket.get('name', 'pocket'))
                    center = pocket.get('center')
                    size = pocket.get('size')

                    if pockets and len(pockets) > 1:
                        print(f"  - Pocket {pi + 1}: {pocket_name}")

                    output_dir, score = self.run_vina_docking(
                        receptor['pdbqt'],
                        ligand['pdbqt'],
                        receptor['name'],
                        ligand_run_name,
                        receptor_id=receptor.get('id', ''),
                        center=center,
                        size=size,
                        pocket_name=pocket_name,
                        pocket_index=pi,
                    )

                    if score is None:
                        if pockets and len(pockets) > 1:
                            print("    ✗ 该口袋对接失败或无法解析得分")
                        continue

                    if pockets and len(pockets) > 1:
                        print(f"    ✓ 得分: {score:.2f} kcal/mol")

                    if (best_score is None) or (score < best_score):
                        best_score = score
                        best_output_dir = output_dir
                        best_pocket = pocket_name
                        best_pocket_index = pi + 1

                if best_score is not None:
                    print(f"  ✓ 最佳得分: {best_score:.2f} kcal/mol")
                    print(f"    最佳结果目录: {best_output_dir}")
                    results.append({
                        'Receptor': receptor['name'],
                        'Receptor_ID': receptor.get('id', ''),
                        'Ligand': parent_name,
                        'Ligand_Variant': ligand_run_name if parent_name != ligand_run_name else '',
                        'Ligand_SMILES': ligand.get('smiles', ''),
                        'Pocket': best_pocket if best_pocket is not None else '',
                        'Pocket_Index': best_pocket_index if best_pocket_index is not None else '',
                        'Affinity': best_score,
                        'Output_Dir': best_output_dir
                    })
                else:
                    print("  ✗ 所有口袋对接均失败或无法解析得分")

        results_df = pd.DataFrame(results)
        results_file = self.results_dir / 'docking_results.csv'
        results_df.to_csv(results_file, index=False, encoding='utf-8-sig')

        enum_enable = bool(self.config.get('ligand_enumeration_enable', False))
        if enum_enable and (not results_df.empty) and ('Ligand_Variant' in results_df.columns):
            has_variants = (results_df['Ligand_Variant'].astype(str).str.strip() != '').any()
            if has_variants:
                try:
                    grouped = results_df.dropna(subset=['Affinity']).copy()
                    idx = grouped.groupby(['Receptor', 'Receptor_ID', 'Ligand'])['Affinity'].idxmin()
                    best_df = grouped.loc[idx].reset_index(drop=True)
                    best_file = self.results_dir / 'docking_results_best.csv'
                    best_df.to_csv(best_file, index=False, encoding='utf-8-sig')
                    print(f"最佳变体结果已保存至: {best_file}")
                    results_df = best_df
                except Exception:
                    pass

        print("\n" + "=" * 60)
        print(f"批量对接完成！成功对接: {len(results)}/{total_dockings}")
        print(f"结果已保存至: {results_file}")
        print("=" * 60 + "\n")

        return results_df

    def _prepare_receptors_from_csv(self, receptor_csv: str) -> List[Dict]:
        """从 CSV 下载并准备受体"""
        receptor_df = pd.read_csv(receptor_csv)
        receptor_df = receptor_df.dropna(subset=['ID'])

        print(f"受体数量: {len(receptor_df)}")
        print("\n步骤 1/3: 下载并准备受体（CSV）")
        print("-" * 60)

        receptors = []
        for _, row in receptor_df.iterrows():
            pdb_id = str(row['ID']).strip()
            pro_name = str(row.get('PRO', pdb_id)).strip()

            cif_file = self.download_pdb_structure(pdb_id, pro_name)
            if cif_file is None:
                continue

            pdbqt_file = self.prepare_receptor(cif_file)
            if pdbqt_file:
                receptors.append({
                    'name': pro_name,
                    'id': pdb_id,
                    'pdbqt': pdbqt_file
                })

        print(f"\n成功准备 {len(receptors)}/{len(receptor_df)} 个受体\n")
        return receptors

    def _prepare_ligands_from_csv(self, ligand_csv: str) -> List[Dict]:
        """从 CSV（SMILES）准备配体"""
        ligand_df = pd.read_csv(ligand_csv)
        ligand_df = ligand_df.dropna(subset=['SMILES'])

        print(f"配体数量: {len(ligand_df)}\n")
        print("步骤 2/3: 准备配体（CSV/SMILES）")
        print("-" * 60)

        ligands = []
        for _, row in ligand_df.iterrows():
            che_name = str(row.get('CHE', 'ligand')).strip()
            smiles = str(row['SMILES']).strip()

            entries = self.prepare_ligand_entries(smiles, che_name)
            for ent in entries:
                ligands.append(ent)

        print(f"\n成功准备 {len(ligands)}/{len(ligand_df)} 个配体\n")
        return ligands

    def _prepare_receptors_from_local(self, receptor_local_path: str) -> List[Dict]:
        """从本地文件夹/压缩包准备受体"""
        src = self._normalize_local_source(receptor_local_path)
        root = self._resolve_local_input_root(receptor_local_path, kind='receptors')

        exts = self.config.get('receptor_local_extensions', ['.pdb', '.cif', '.mmcif'])
        recursive = self.config.get('local_recursive', True)

        if src.is_file() and src.suffix.lower() in {'.pdb', '.cif', '.mmcif'}:
            receptor_files = [src]
        else:
            receptor_files = self._collect_files_by_ext(root, exts=exts, recursive=recursive)

        print(f"受体文件数量: {len(receptor_files)}")
        print("\n步骤 1/3: 准备受体（本地文件/压缩包）")
        print("-" * 60)

        receptors = []
        for fp in receptor_files:
            pdbqt_file = self.prepare_receptor_from_structure_file(str(fp))
            if pdbqt_file:
                name = self._safe_stem(fp)
                receptors.append({
                    'name': name,
                    'id': name,
                    'pdbqt': pdbqt_file
                })

        print(f"\n成功准备 {len(receptors)}/{len(receptor_files)} 个受体\n")
        return receptors

    def _prepare_ligands_from_local(self, ligand_local_path: str) -> List[Dict]:
        """从本地文件夹/压缩包准备配体"""
        src = self._normalize_local_source(ligand_local_path)
        root = self._resolve_local_input_root(ligand_local_path, kind='ligands')

        exts = self.config.get('ligand_local_extensions', ['.pdb', '.mol', '.mol2', '.sdf', '.pdbqt'])
        recursive = self.config.get('local_recursive', True)

        if src.is_file() and src.suffix.lower() in {'.pdb', '.mol', '.mol2', '.sdf', '.pdbqt'}:
            ligand_files = [src]
        else:
            ligand_files = self._collect_files_by_ext(root, exts=exts, recursive=recursive)

        print(f"配体文件数量: {len(ligand_files)}\n")
        print("步骤 2/3: 准备配体（本地文件/压缩包）")
        print("-" * 60)

        ligands = []
        for fp in ligand_files:
            print(f"正在准备配体: {fp.name} (本地文件)...")
            raw_file, pdbqt_file = self.prepare_ligand_from_file(str(fp))
            if pdbqt_file:
                ligands.append({
                    'name': self._safe_stem(fp),
                    'parent_name': self._safe_stem(fp),
                    'pdb': raw_file,
                    'pdbqt': pdbqt_file
                })
            else:
                print("  ✗ 配体准备失败")

        print(f"\n成功准备 {len(ligands)}/{len(ligand_files)} 个配体\n")
        return ligands

    def batch_docking_from_config(self) -> pd.DataFrame:
        """根据配置自动选择输入模式并执行批量对接"""
        print("\n" + "=" * 60)
        print("开始批量分子对接")
        print("=" * 60 + "\n")

        receptor_mode = str(self.config.get('receptor_input_mode', 'csv')).strip().lower()
        ligand_mode = str(self.config.get('ligand_input_mode', 'csv')).strip().lower()

        receptors: List[Dict]
        ligands: List[Dict]

        if receptor_mode == 'local':
            receptor_local_path = self.config.get('receptor_local_path', '')
            if not receptor_local_path:
                raise ValueError('receptor_input_mode=local 时必须设置 receptor_local_path')
            receptors = self._prepare_receptors_from_local(receptor_local_path)
        else:
            receptor_csv = self.config.get('receptor_csv', '测试受体.csv')
            if not os.path.isabs(receptor_csv):
                receptor_csv = os.path.abspath(receptor_csv)
            receptors = self._prepare_receptors_from_csv(receptor_csv)

        if ligand_mode == 'local':
            ligand_local_path = self.config.get('ligand_local_path', '')
            if not ligand_local_path:
                raise ValueError('ligand_input_mode=local 时必须设置 ligand_local_path')
            ligands = self._prepare_ligands_from_local(ligand_local_path)
        else:
            ligand_csv = self.config.get('ligand_csv', '测试配体.csv')
            if not os.path.isabs(ligand_csv):
                ligand_csv = os.path.abspath(ligand_csv)
            ligands = self._prepare_ligands_from_csv(ligand_csv)

        if not receptors:
            print("✗ 未成功准备任何受体，终止")
            return pd.DataFrame([])
        if not ligands:
            print("✗ 未成功准备任何配体，终止")
            return pd.DataFrame([])

        return self._run_docking_matrix(receptors, ligands)
    
    def batch_docking(self, receptor_csv: str, ligand_csv: str) -> pd.DataFrame:
        """
        批量对接
        
        Args:
            receptor_csv: 受体CSV文件路径
            ligand_csv: 配体CSV文件路径
            
        Returns:
            对接结果DataFrame
        """
        print("\n" + "="*60)
        print("开始批量分子对接")
        print("="*60 + "\n")
        
        # 读取CSV文件
        receptor_df = pd.read_csv(receptor_csv)
        ligand_df = pd.read_csv(ligand_csv)
        
        # 清理空行
        receptor_df = receptor_df.dropna(subset=['ID'])
        ligand_df = ligand_df.dropna(subset=['SMILES'])
        
        print(f"受体数量: {len(receptor_df)}")
        print(f"配体数量: {len(ligand_df)}\n")
        
        # 步骤1：下载并准备受体
        print("步骤 1/3: 下载并准备受体")
        print("-" * 60)
        receptors = []
        for idx, row in receptor_df.iterrows():
            pdb_id = row['ID'].strip()
            pro_name = row['PRO'].strip()
            
            # 下载CIF文件
            cif_file = self.download_pdb_structure(pdb_id, pro_name)
            if cif_file is None:
                continue
            
            # 准备受体
            pdbqt_file = self.prepare_receptor(cif_file)
            if pdbqt_file:
                receptors.append({
                    'name': pro_name,
                    'id': pdb_id,
                    'pdbqt': pdbqt_file
                })
        
        print(f"\n成功准备 {len(receptors)}/{len(receptor_df)} 个受体\n")
        
        # 步骤2：准备配体
        print("步骤 2/3: 准备配体")
        print("-" * 60)
        ligands = []
        for idx, row in ligand_df.iterrows():
            che_name = row['CHE'].strip()
            smiles = row['SMILES'].strip()

            entries = self.prepare_ligand_entries(smiles, che_name)
            for ent in entries:
                ligands.append(ent)
        
        print(f"\n成功准备 {len(ligands)}/{len(ligand_df)} 个配体\n")
        
        return self._run_docking_matrix(receptors, ligands)
    
    def plot_heatmap(self, results_df: pd.DataFrame, output_base: str = None):
        """
        绘制对接结果热图
        
        Args:
            results_df: 对接结果DataFrame
            output_base: 输出文件路径（不含扩展名）
        """
        if output_base is None:
            output_base = str(self.results_dir / 'docking_heatmap')
        else:
            output_base = str(output_base).replace('.png', '').replace('.pdf', '')
        
        print("正在绘制对接结果热图...")

        plot_df = results_df.copy()
        plot_df['Receptor_Label'] = plot_df['Receptor']

        # 创建数据透视表（遇到重复键时取最优亲和力）
        pivot_data = plot_df.pivot_table(
            index='Ligand',
            columns='Receptor_Label',
            values='Affinity',
            aggfunc='min'
        )
        
        # 设置图形大小
        fig_width = max(10, len(pivot_data.columns) * 1.2)
        fig_height = max(8, len(pivot_data.index) * 0.5)
        plt.figure(figsize=(fig_width, fig_height))
        
        # 获取配色方案
        colormap = self.config.get('heatmap_colormap', 'RdYlGn_r')
        
        # 绘制热图
        sns.heatmap(
            pivot_data,
            annot=True,
            fmt='.2f',
            cmap=colormap,
            cbar_kws={'label': 'Binding Affinity (kcal/mol)'},
            linewidths=0.5,
            linecolor='gray'
        )
        
        plt.title('Molecular Docking Results Heatmap', fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Receptor', fontsize=12, fontweight='bold')
        plt.ylabel('Ligand', fontsize=12, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        plt.yticks(rotation=0)
        plt.tight_layout()
        
        # 保存PNG格式
        png_file = f"{output_base}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', format='png')
        print(f"  ✓ 热图已保存 (PNG): {png_file}")
        
        # 保存PDF格式
        pdf_file = f"{output_base}.pdf"
        plt.savefig(pdf_file, dpi=300, bbox_inches='tight', format='pdf')
        print(f"  ✓ 热图已保存 (PDF): {pdf_file}")
        
        plt.close()
        
        # 额外绘制最佳对接结果柱状图
        self.plot_top_results(results_df)
    
    def plot_top_results(self, results_df: pd.DataFrame, top_n: int = 10):
        """绘制最佳对接结果柱状图"""
        output_base = str(self.results_dir / 'top_docking_results')
        
        # 按亲和力排序，选择top N
        top_results = results_df.nsmallest(top_n, 'Affinity')
        
        plt.figure(figsize=(12, 6))
        
        # 创建标签
        labels = [f"{row['Ligand']}\nvs\n{row['Receptor']}" for _, row in top_results.iterrows()]
        
        colors = plt.cm.RdYlGn_r(np.linspace(0.3, 0.9, len(top_results)))
        
        bars = plt.bar(range(len(top_results)), top_results['Affinity'], color=colors)
        
        plt.xlabel('Receptor-Ligand Pair', fontsize=12, fontweight='bold')
        plt.ylabel('Binding Affinity (kcal/mol)', fontsize=12, fontweight='bold')
        plt.title(f'Top {top_n} Docking Results', fontsize=14, fontweight='bold', pad=20)
        plt.xticks(range(len(top_results)), labels, rotation=45, ha='right', fontsize=9)
        plt.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
        plt.grid(axis='y', alpha=0.3)
        
        # 在柱子上标注数值
        for i, (bar, val) in enumerate(zip(bars, top_results['Affinity'])):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height(), 
                    f'{val:.2f}', ha='center', va='bottom', fontsize=9)
        
        plt.tight_layout()
        
        # 保存PNG格式
        png_file = f"{output_base}.png"
        plt.savefig(png_file, dpi=300, bbox_inches='tight', format='png')
        print(f"  ✓ Top结果图已保存 (PNG): {png_file}")
        
        # 保存PDF格式
        pdf_file = f"{output_base}.pdf"
        plt.savefig(pdf_file, dpi=300, bbox_inches='tight', format='pdf')
        print(f"  ✓ Top结果图已保存 (PDF): {pdf_file}")
        
        plt.close()


def main():
    """主函数"""
    # 配置文件路径
    config_file = "docking_config_advanced.json" if os.path.exists("docking_config_advanced.json") else "docking_config.json"
    
    # 创建对接系统实例
    docking_system = MolecularDocking(config_file)
    
    # 根据配置文件自动选择输入模式并执行批量对接
    results_df = docking_system.batch_docking_from_config()
    
    # 绘制热图
    if len(results_df) > 0:
        docking_system.plot_heatmap(results_df)
        print("\n✓ 所有任务完成！")
    else:
        print("\n✗ 没有成功的对接结果，无法绘制热图")


if __name__ == "__main__":
    main()
