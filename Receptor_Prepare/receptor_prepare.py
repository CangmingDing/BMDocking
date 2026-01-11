#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import shutil
import tarfile
import zipfile
import subprocess
from pathlib import Path
from datetime import datetime
from typing import List, Tuple, Dict, Optional

import pandas as pd
import requests


class ReceptorPrepare:
    def __init__(self, config_file: str = "docking_config.json", work_dir: Optional[str] = None):
        self.config = self.load_config(config_file)
        if work_dir is not None and str(work_dir).strip():
            self.config['work_dir'] = str(work_dir).strip()

        self.work_dir = Path(self.config.get('work_dir', './docking_workspace'))
        self.receptor_raw_dir = self.work_dir / 'receptors' / 'raw'
        self.receptor_prepared_dir = self.work_dir / 'receptors' / 'prepared'
        self.inputs_extracted_dir = self.work_dir / 'inputs' / 'extracted' / 'receptors'
        for dir_path in [self.receptor_raw_dir, self.receptor_prepared_dir, self.inputs_extracted_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)

    def load_config(self, config_file: str) -> dict:
        if os.path.exists(config_file):
            with open(config_file, 'r', encoding='utf-8') as f:
                return json.load(f)

        return {
            'work_dir': './docking_workspace',
            'receptor_method': 'openbabel',
            'receptor_cif_to_pdb_fallback_enable': True,
            'receptor_cif_to_pdb_preferred_tool': 'openbabel',
            'receptor_cif_to_pdb_openbabel_timeout_sec': 60,
            'local_recursive': True,
            'local_extract_archives': True,
            'receptor_cleanup': {
                'enable': True,
                'enable_altloc': True,
                'remove_waters': True,
                'remove_ions': True,
                'remove_buffers': True,
                'remove_organic_ligands': True,
                'keep_hetero_resnames': []
            },
            'pdbfixer': {
                'enable': True,
                'add_missing_residues': True,
                'add_missing_atoms': True,
                'add_hydrogens': True,
                'ph': 7.4
            }
        }

    def _safe_stem(self, file_path: Path) -> str:
        stem = file_path.stem
        safe = "".join(c for c in stem if c.isalnum() or c in (' ', '-', '_')).strip()
        return safe if safe else stem

    def _looks_like_pdb(self, pdb_path: Path) -> bool:
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

    def _normalize_local_source(self, source_path: str) -> Path:
        expanded = os.path.expanduser(str(source_path)).strip()
        return Path(expanded).resolve()

    def _extract_archive_to_dir(self, archive_path: Path, extract_root: Path) -> Path:
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

    def _resolve_local_input_root(self, source: str) -> Tuple[Path, Optional[Path]]:
        source_path = self._normalize_local_source(source)
        if not source_path.exists():
            raise FileNotFoundError(f"本地输入路径不存在: {source_path}")

        if source_path.is_dir():
            return source_path, None

        if source_path.is_file():
            lower_name = source_path.name.lower()
            is_archive = any(
                lower_name.endswith(suffix)
                for suffix in ['.zip', '.tar', '.tgz', '.tar.gz', '.tar.bz2', '.tar.xz']
            )
            if is_archive and self.config.get('local_extract_archives', True):
                out_dir = self._extract_archive_to_dir(source_path, self.inputs_extracted_dir)
                return out_dir, None
            return source_path.parent, source_path

        raise ValueError(f"无法识别的输入路径: {source_path}")

    def _collect_files_by_ext(self, root: Path, exts: List[str], recursive: bool) -> List[Path]:
        normalized_exts = set(e.lower() if e.startswith('.') else f".{e.lower()}" for e in exts)
        pattern = '**/*' if recursive else '*'
        files = []
        for p in root.glob(pattern):
            if p.is_file() and p.suffix.lower() in normalized_exts:
                files.append(p)
        files.sort(key=lambda x: x.name.lower())
        return files

    def _convert_cif_to_pdb_openbabel(self, cif_file: Path, pdb_file: Path, timeout_sec: int) -> Tuple[bool, str]:
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

    def preprocess_receptor_pdb(self, pdb_file: str) -> str:
        cleanup_cfg = self.config.get('receptor_cleanup', {}) or {}
        pdbfixer_cfg = self.config.get('pdbfixer', {}) or {}

        current_pdb = Path(pdb_file)

        try:
            if cleanup_cfg.get('enable', True) and cleanup_cfg.get('enable_altloc', True):
                altloc_pdb = current_pdb.with_name(current_pdb.stem + '_altloc.pdb')
                self._select_altloc_highest_occupancy(str(current_pdb), str(altloc_pdb))
                current_pdb = altloc_pdb

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
        from collections import defaultdict

        with open(input_pdb, 'r') as f:
            lines = f.readlines()

        groups = defaultdict(list)

        for idx, line in enumerate(lines):
            if not line.startswith(('ATOM', 'HETATM')):
                continue
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

        best_altloc_by_key = {}
        for key, entries in groups.items():
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

        output_lines = []
        for _, line in enumerate(lines):
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

            if chosen_alt is None:
                output_lines.append(line)
                continue

            this_alt = altloc.strip() or ' '
            if this_alt != chosen_alt:
                continue

            if altloc != ' ':
                line = line[:16] + ' ' + line[17:]

            output_lines.append(line)

        with open(output_pdb, 'w') as f:
            f.writelines(output_lines)

    def _cleanup_solvents_and_ligands(self, input_pdb: str, output_pdb: str, cleanup_cfg: dict) -> None:
        keep_hetero = set(
            res.strip().upper()
            for res in cleanup_cfg.get('keep_hetero_resnames', [])
            if isinstance(res, str) and res.strip()
        )

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
                    output_lines.append(line)
                    continue

                record_name = line[:6].strip().upper()
                resname = line[17:20].strip().upper()

                if record_name == 'ATOM':
                    output_lines.append(line)
                    continue

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
                    continue

                output_lines.append(line)

        with open(output_pdb, 'w') as f:
            f.writelines(output_lines)

    def _sanitize_pdb_numeric_fields_for_pdbfixer(self, pdb_file: str) -> None:
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

    def _fix_missing_residues_and_atoms(
        self,
        input_pdb: str,
        output_pdb: str,
        add_missing_residues: bool = True,
        add_missing_atoms: bool = True,
        add_hydrogens: bool = True,
        ph: float = 7.4
    ) -> bool:
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

    def prepare_receptor_meeko(self, pdb_file: str, pdbqt_file: str) -> bool:
        try:
            from meeko import PDBQTWriterLegacy, MoleculePreparation
            from rdkit import Chem

            mol = Chem.MolFromPDBFile(pdb_file, removeHs=False)
            if mol is None:
                return False
            mol = Chem.AddHs(mol, addCoords=True)

            preparator = MoleculePreparation()
            mol_setups = preparator.prepare(mol)

            for setup in mol_setups:
                pdbqt_string = PDBQTWriterLegacy.write_string(setup)
                with open(pdbqt_file, 'w') as f:
                    f.write(pdbqt_string[0] if isinstance(pdbqt_string, tuple) else pdbqt_string)
                break

            return os.path.exists(pdbqt_file)
        except Exception as e:
            print(f"    Meeko处理受体失败: {e}")
            return False

    def prepare_receptor_from_structure_file(self, structure_file: str, method: Optional[str] = None) -> Optional[str]:
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

            pdb_file = Path(self.preprocess_receptor_pdb(str(pdb_file)))

            if str(method).lower() == 'meeko':
                success = self.prepare_receptor_meeko(str(pdb_file), str(pdbqt_file))
                if not success:
                    print("  Meeko失败，回退到OpenBabel...")
                    method = 'openbabel'

            if str(method).lower() == 'openbabel':
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

            print(f"  ✓ 受体准备完成: {pdbqt_file}")
            return str(pdbqt_file)
        except Exception as e:
            print(f"  ✗ 受体准备失败: {e}")
            return None

    def download_pdb_structure(self, pdb_id: str, pro_name: str) -> Optional[str]:
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
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

    def prepare_receptor(self, cif_file: str, method: Optional[str] = None) -> Optional[str]:
        if method is None:
            method = self.config.get('receptor_method', 'openbabel')

        base_name = Path(cif_file).stem
        pdb_file = self.receptor_prepared_dir / f"{base_name}.pdb"
        pdbqt_file = self.receptor_prepared_dir / f"{base_name}.pdbqt"

        print(f"正在准备受体: {base_name} (使用{method})...")

        try:
            ok, msg = self._convert_cif_to_pdb_with_fallback(Path(cif_file), pdb_file)
            if not ok:
                print(f"  ✗ CIF转PDB失败: {msg[:200]}")
                return None

            pdb_file = Path(self.preprocess_receptor_pdb(str(pdb_file)))

            if str(method).lower() == 'meeko':
                success = self.prepare_receptor_meeko(str(pdb_file), str(pdbqt_file))
                if not success:
                    print("  Meeko失败，回退到OpenBabel...")
                    method = 'openbabel'

            if str(method).lower() == 'openbabel':
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

            print(f"  ✓ 受体准备完成: {pdbqt_file}")
            return str(pdbqt_file)
        except Exception as e:
            print(f"  ✗ 受体准备失败: {e}")
            return None

    def prepare_receptors_from_csv(self, receptor_csv: str, method: Optional[str] = None) -> List[Dict]:
        receptor_df = pd.read_csv(receptor_csv)
        receptor_df = receptor_df.dropna(subset=['ID'])

        receptors = []
        for _, row in receptor_df.iterrows():
            pdb_id = str(row['ID']).strip()
            pro_name = str(row.get('PRO', pdb_id)).strip()

            cif_file = self.download_pdb_structure(pdb_id, pro_name)
            if cif_file is None:
                continue

            pdbqt_file = self.prepare_receptor(cif_file, method=method)
            if pdbqt_file:
                receptors.append({'name': pro_name, 'id': pdb_id, 'pdbqt': pdbqt_file})

        return receptors

    def prepare_receptors_from_local(self, local_path: str, method: Optional[str] = None) -> List[Dict]:
        root, only_file = self._resolve_local_input_root(local_path)
        recursive = bool(self.config.get('local_recursive', True))
        candidates = self._collect_files_by_ext(root, exts=['.pdb', '.cif', '.mmcif'], recursive=recursive)

        if only_file is not None:
            only_file_abs = only_file.resolve()
            candidates = [p for p in candidates if p.resolve() == only_file_abs]

        receptors = []
        for p in candidates:
            pdbqt_file = self.prepare_receptor_from_structure_file(str(p), method=method)
            if pdbqt_file:
                receptors.append({'name': self._safe_stem(p), 'id': '', 'pdbqt': pdbqt_file})
        return receptors


def _parse_args(argv: List[str]):
    import argparse

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('--config', default='docking_config.json')
    parser.add_argument('--work_dir', default=None)
    parser.add_argument('--mode', choices=['csv', 'local'], required=True)
    parser.add_argument('--receptor_csv', default=None)
    parser.add_argument('--receptor_local_path', default=None)
    parser.add_argument('--method', choices=['openbabel', 'meeko'], default=None)
    return parser.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]

    args = _parse_args(argv)
    system = ReceptorPrepare(config_file=args.config, work_dir=args.work_dir)

    mode = args.mode or str(system.config.get('mode', '')).strip().lower()
    receptor_csv = args.receptor_csv or system.config.get('receptor_csv')
    receptor_local_path = args.receptor_local_path or system.config.get('receptor_local_path')
    method = args.method or system.config.get('receptor_method')

    if mode == 'csv':
        if not receptor_csv:
            raise ValueError('mode=csv 时必须提供 --receptor_csv 或在配置中提供 receptor_csv')
        receptors = system.prepare_receptors_from_csv(str(receptor_csv), method=method)
    elif mode == 'local':
        if not receptor_local_path:
            raise ValueError('mode=local 时必须提供 --receptor_local_path 或在配置中提供 receptor_local_path')
        receptors = system.prepare_receptors_from_local(str(receptor_local_path), method=method)
    else:
        raise ValueError('必须通过 --mode 指定 csv 或 local，或在配置中提供 mode')

    print('\n' + '=' * 60)
    print(f"受体准备完成：{len(receptors)} 个")
    print(f"输出目录：{system.receptor_prepared_dir}")
    print('=' * 60)

    return 0


if __name__ == '__main__':
    raise SystemExit(main())
