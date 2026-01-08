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
from typing import List, Tuple, Dict
import warnings
warnings.filterwarnings('ignore')

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS', 'SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False


class MolecularDocking:
    """分子对接主类"""
    
    def __init__(self, config_file: str = "docking_config.json"):
        """
        初始化对接系统
        
        Args:
            config_file: 配置文件路径
        """
        self.config = self.load_config(config_file)
        self.work_dir = Path(self.config.get('work_dir', './docking_workspace'))
        
        # 受体目录结构
        self.receptor_dir = self.work_dir / 'receptors'
        self.receptor_raw_dir = self.receptor_dir / 'raw'
        self.receptor_prepared_dir = self.receptor_dir / 'prepared'
        
        # 配体目录结构
        self.ligand_dir = self.work_dir / 'ligands'
        self.ligand_raw_dir = self.ligand_dir / 'raw'
        self.ligand_prepared_dir = self.ligand_dir / 'prepared'
        
        # 结果目录
        self.results_dir = self.work_dir / 'results'
        
        # 创建所有目录
        for dir_path in [
            self.receptor_raw_dir, self.receptor_prepared_dir,
            self.ligand_raw_dir, self.ligand_prepared_dir,
            self.results_dir
        ]:
            dir_path.mkdir(parents=True, exist_ok=True)
    
    def load_config(self, config_file: str) -> dict:
        """加载配置文件"""
        if os.path.exists(config_file):
            with open(config_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        else:
            # 默认配置
            return {
                'work_dir': './docking_workspace',
                'ligand_method': 'meeko',  # 'openbabel' or 'meeko'
                'box_center': [0, 0, 0],
                'box_size': [20, 20, 20],
                'exhaustiveness': 8,
                'num_modes': 9,
                'energy_range': 3
            }
    
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
            # 1. CIF转PDB（路径加引号处理空格）
            cmd = f'obabel -icif "{cif_file}" -opdb -O "{pdb_file}" -d'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            if result.returncode != 0:
                print(f"  ✗ CIF转PDB失败: {result.stderr[:200]}")
                return None
            
            # 2. 根据方法选择处理方式
            if method.lower() == 'meeko':
                # 使用Meeko准备受体
                success = self.prepare_receptor_meeko(str(pdb_file), str(pdbqt_file))
                if not success:
                    print("  Meeko失败，回退到OpenBabel...")
                    method = 'openbabel'
            
            if method.lower() == 'openbabel':
                # 使用OpenBabel准备受体
                # 首先尝试ADT的prepare_receptor4
                try:
                    cmd = f'prepare_receptor4 -r "{pdb_file}" -o "{pdbqt_file}" -A hydrogens'
                    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
                    if result.returncode != 0:
                        raise Exception("prepare_receptor4 failed")
                except:
                    # 备选方案：使用obabel
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
    
    def smiles_to_3d_openbabel(self, smiles: str, output_file: str) -> bool:
        """
        使用OpenBabel将SMILES转换为3D结构并准备配体
        
        Args:
            smiles: SMILES字符串
            output_file: 输出文件路径（不含扩展名）
            
        Returns:
            是否成功
        """
        try:
            # raw目录中的PDB文件
            pdb_file_raw = f"{output_file}_raw.pdb"
            # prepared目录中的最终PDBQT文件
            pdbqt_file = f"{output_file}.pdbqt"
            
            # 删除旧文件避免误判
            for f in [pdb_file_raw, pdbqt_file]:
                if os.path.exists(f):
                    os.remove(f)
            
            # SMILES -> 3D PDB (使用更安全的引号处理)
            # 将SMILES用双引号包裹，避免shell解析问题
            smiles_escaped = smiles.replace('"', '\\"')
            cmd = f'obabel -:"{smiles_escaped}" -opdb -O "{pdb_file_raw}" --gen3d --addh --minimize --steps 200 --ff MMFF94'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
            
            if result.returncode != 0 or not os.path.exists(pdb_file_raw) or os.path.getsize(pdb_file_raw) == 0:
                print(f"    OpenBabel生成3D结构失败")
                if result.stderr:
                    print(f"    错误信息: {result.stderr[:200]}")
                return False
            
            # PDB -> PDBQT
            cmd = f'obabel -ipdb "{pdb_file_raw}" -opdbqt -O "{pdbqt_file}" -xn'
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=30)
            
            if result.returncode != 0 or not os.path.exists(pdbqt_file) or os.path.getsize(pdbqt_file) == 0:
                print(f"    OpenBabel转换PDBQT失败")
                if result.stderr:
                    print(f"    错误信息: {result.stderr[:200]}")
                return False
            
            return True
            
        except subprocess.TimeoutExpired:
            print(f"    处理超时（分子可能太复杂）")
            return False
        except Exception as e:
            print(f"    OpenBabel处理失败: {e}")
            return False
    
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
    
    def prepare_ligand(self, smiles: str, ligand_name: str) -> str:
        """
        准备配体文件（固定使用OpenBabel）
        
        Args:
            smiles: SMILES字符串
            ligand_name: 配体名称
            
        Returns:
            PDBQT文件路径
        """
        # 清理文件名
        safe_name = "".join(c for c in ligand_name if c.isalnum() or c in (' ', '-', '_')).strip()
        if not safe_name:
            safe_name = "ligand"
        
        # raw文件保存到raw目录，prepared文件保存到prepared目录
        output_base_raw = self.ligand_raw_dir / safe_name
        output_base_prepared = self.ligand_prepared_dir / safe_name
        pdbqt_file = f"{output_base_prepared}.pdbqt"
        
        print(f"正在准备配体: {ligand_name} (使用OpenBabel)...")
        
        # 固定使用OpenBabel
        # 传入raw和prepared的基础路径
        success = self.smiles_to_3d_openbabel(smiles, str(output_base_raw), str(output_base_prepared))
        
        if success and os.path.exists(pdbqt_file):
            print(f"  ✓ 配体准备完成: {pdbqt_file}")
            return pdbqt_file
        else:
            print(f"  ✗ 配体准备失败")
            return None
    
    def run_vina_docking(self, receptor_pdbqt: str, ligand_pdbqt: str, 
                        receptor_name: str, ligand_name: str) -> Tuple[str, float]:
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
        # 为每个对接对创建独立目录
        docking_pair_name = f"{receptor_name}_vs_{ligand_name}"
        pair_dir = self.results_dir / docking_pair_name
        pair_dir.mkdir(parents=True, exist_ok=True)
        
        # 输出文件
        output_poses = pair_dir / "docked_poses.pdbqt"  # 所有对接姿态
        log_file = pair_dir / "vina_log.txt"  # Vina日志
        best_pose = pair_dir / "best_pose.pdbqt"  # 最佳姿态
        
        # 构建Vina命令
        center = self.config['box_center']
        size = self.config['box_size']
        
        cmd = (
            f'vina --receptor "{receptor_pdbqt}" '
            f'--ligand "{ligand_pdbqt}" '
            f'--out "{output_poses}" '
            f"--center_x {center[0]} --center_y {center[1]} --center_z {center[2]} "
            f"--size_x {size[0]} --size_y {size[1]} --size_z {size[2]} "
            f"--exhaustiveness {self.config['exhaustiveness']} "
            f"--num_modes {self.config['num_modes']} "
            f"--energy_range {self.config['energy_range']}"
        )
        
        try:
            # 运行Vina
            result = subprocess.run(cmd, shell=True, capture_output=True, 
                                   text=True, timeout=300)
            
            # 保存完整日志
            with open(log_file, 'w') as f:
                f.write("=== VINA OUTPUT ===\n")
                f.write(result.stdout)
                if result.stderr:
                    f.write("\n=== STDERR ===\n")
                    f.write(result.stderr)
            
            # 检查是否成功
            if result.returncode != 0:
                print(f"  ✗ Vina执行失败: {result.stderr[:200]}")
                return None, None
            
            # 解析对接得分
            score = self.parse_vina_output(result.stdout)
            
            # 提取最佳姿态（第一个模型）
            if os.path.exists(output_poses):
                self.extract_best_pose(output_poses, best_pose)
                
                # 生成复合物PDB（受体+最佳配体姿态）
                self.create_complex(receptor_pdbqt, best_pose, pair_dir / "complex.pdb")
            
            return str(pair_dir), score
            
        except subprocess.TimeoutExpired:
            print(f"  ✗ 对接超时")
            return None, None
        except Exception as e:
            print(f"  ✗ 对接失败: {e}")
            return None, None
    
    def parse_vina_output(self, output_text: str) -> float:
        """解析Vina输出获取最佳对接得分"""
        try:
            lines = output_text.split('\n')
            for i, line in enumerate(lines):
                # 查找结果表格的标题行
                if 'affinity' in line.lower() and 'kcal/mol' in line.lower():
                    # 跳过标题行和分隔线，查找第一行数据
                    for j in range(i + 1, min(i + 10, len(lines))):
                        score_line = lines[j].strip()
                        # 跳过空行和分隔线
                        if not score_line or score_line.startswith('-') or score_line.startswith('|'):
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
            
            pdbqt_file = self.prepare_ligand(smiles, che_name)
            if pdbqt_file:
                ligands.append({
                    'name': che_name,
                    'smiles': smiles,
                    'pdbqt': pdbqt_file
                })
        
        print(f"\n成功准备 {len(ligands)}/{len(ligand_df)} 个配体\n")
        
        # 步骤3：批量对接
        print("步骤 3/3: 执行批量对接")
        print("-" * 60)
        results = []
        total_dockings = len(receptors) * len(ligands)
        current = 0
        
        for receptor in receptors:
            for ligand in ligands:
                current += 1
                print(f"\n[{current}/{total_dockings}] 对接: {receptor['name']} vs {ligand['name']}")
                
                output_prefix = f"{receptor['id']}_{ligand['name']}"
                output_file, score = self.run_vina_docking(
                    receptor['pdbqt'],
                    ligand['pdbqt'],
                    output_prefix
                )
                
                if score is not None:
                    print(f"  ✓ 对接完成，得分: {score:.2f} kcal/mol")
                    results.append({
                        'Receptor': receptor['name'],
                        'Receptor_ID': receptor['id'],
                        'Ligand': ligand['name'],
                        'Affinity': score,
                        'Output': output_file
                    })
                else:
                    print(f"  ✗ 对接失败或无法解析得分")
        
        # 保存结果
        results_df = pd.DataFrame(results)
        results_file = self.results_dir / 'docking_results.csv'
        results_df.to_csv(results_file, index=False, encoding='utf-8-sig')
        
        print("\n" + "="*60)
        print(f"批量对接完成！成功对接: {len(results)}/{total_dockings}")
        print(f"结果已保存至: {results_file}")
        print("="*60 + "\n")
        
        return results_df
    
    def plot_heatmap(self, results_df: pd.DataFrame, output_file: str = None):
        """
        绘制对接结果热图
        
        Args:
            results_df: 对接结果DataFrame
            output_file: 输出图片文件路径
        """
        if output_file is None:
            output_file = self.results_dir / 'docking_heatmap.png'
        
        print("正在绘制对接结果热图...")
        
        # 创建数据透视表
        pivot_data = results_df.pivot(
            index='Ligand',
            columns='Receptor',
            values='Affinity'
        )
        
        # 设置图形大小
        fig_width = max(10, len(pivot_data.columns) * 1.2)
        fig_height = max(8, len(pivot_data.index) * 0.5)
        plt.figure(figsize=(fig_width, fig_height))
        
        # 绘制热图
        sns.heatmap(
            pivot_data,
            annot=True,
            fmt='.2f',
            cmap='RdYlGn_r',  # 红黄绿反转（更负=更好=更绿）
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
        
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ 热图已保存: {output_file}")
        plt.close()
        
        # 额外绘制最佳对接结果柱状图
        self.plot_top_results(results_df)
    
    def plot_top_results(self, results_df: pd.DataFrame, top_n: int = 10):
        """绘制最佳对接结果柱状图"""
        output_file = self.results_dir / 'top_docking_results.png'
        
        # 按亲和力排序，选择top N
        top_results = results_df.nsmallest(top_n, 'Affinity')
        
        plt.figure(figsize=(12, 6))
        
        # 创建标签
        labels = [f"{row['Ligand']}\nvs\n{row['Receptor']}" 
                 for _, row in top_results.iterrows()]
        
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
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"  ✓ Top结果图已保存: {output_file}")
        plt.close()


def main():
    """主函数"""
    # 配置文件路径
    config_file = "docking_config.json"
    
    # 创建对接系统实例
    docking_system = MolecularDocking(config_file)
    
    # 从配置文件读取CSV文件路径
    receptor_csv = docking_system.config.get('receptor_csv', '测试受体.csv')
    ligand_csv = docking_system.config.get('ligand_csv', '测试配体.csv')
    
    # 如果是相对路径，转换为绝对路径
    if not os.path.isabs(receptor_csv):
        receptor_csv = os.path.abspath(receptor_csv)
    if not os.path.isabs(ligand_csv):
        ligand_csv = os.path.abspath(ligand_csv)
    
    # 执行批量对接
    results_df = docking_system.batch_docking(receptor_csv, ligand_csv)
    
    # 绘制热图
    if len(results_df) > 0:
        docking_system.plot_heatmap(results_df)
        print("\n✓ 所有任务完成！")
    else:
        print("\n✗ 没有成功的对接结果，无法绘制热图")


if __name__ == "__main__":
    main()
