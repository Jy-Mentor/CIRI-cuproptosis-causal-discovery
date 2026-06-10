#!/usr/bin/env python3
"""
下载PDB结构文件
支持从RCSB PDB数据库下载蛋白质结构
"""

import os
import sys
import urllib.request
import urllib.error
from pathlib import Path


def download_pdb(pdb_id: str, output_dir: str = ".", format: str = "pdb") -> str:
    """
    从RCSB PDB下载结构文件
    
    Parameters:
    -----------
    pdb_id : str
        PDB ID (如: 3RNM)
    output_dir : str
        输出目录路径
    format : str
        文件格式: "pdb" (经典格式) 或 "cif" (mmCIF格式)
    
    Returns:
    --------
    str : 下载文件的完整路径
    """
    pdb_id = pdb_id.upper().strip()
    
    # 确保输出目录存在
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # 构建下载URL
    if format.lower() == "pdb":
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        filename = f"{pdb_id}.pdb"
    elif format.lower() == "cif":
        url = f"https://files.rcsb.org/download/{pdb_id}.cif"
        filename = f"{pdb_id}.cif"
    else:
        raise ValueError(f"不支持的格式: {format}。请使用 'pdb' 或 'cif'")
    
    output_file = output_path / filename
    
    # 下载文件
    try:
        print(f"正在下载 {pdb_id} 的PDB文件...")
        print(f"URL: {url}")
        
        # 设置User-Agent以避免被阻止
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.0'
        }
        
        request = urllib.request.Request(url, headers=headers)
        
        with urllib.request.urlopen(request, timeout=30) as response:
            if response.status == 200:
                with open(output_file, 'wb') as f:
                    f.write(response.read())
                print(f"✓ 成功下载: {output_file}")
                print(f"  文件大小: {output_file.stat().st_size / 1024:.2f} KB")
                return str(output_file)
            else:
                raise urllib.error.HTTPError(
                    url, response.status, 
                    f"HTTP错误: {response.status}", None, None
                )
                
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"✗ 错误: PDB ID '{pdb_id}' 不存在或文件不可用")
        else:
            print(f"✗ HTTP错误 {e.code}: {e.reason}")
        raise
    except Exception as e:
        print(f"✗ 下载失败: {e}")
        raise


def download_pdb_batch(pdb_ids: list, output_dir: str = ".", format: str = "pdb"):
    """
    批量下载多个PDB文件
    
    Parameters:
    -----------
    pdb_ids : list
        PDB ID列表
    output_dir : str
        输出目录
    format : str
        文件格式
    """
    success_count = 0
    failed_ids = []
    
    print(f"开始批量下载 {len(pdb_ids)} 个PDB文件...\n")
    
    for i, pdb_id in enumerate(pdb_ids, 1):
        print(f"[{i}/{len(pdb_ids)}] ", end="")
        try:
            download_pdb(pdb_id, output_dir, format)
            success_count += 1
        except Exception as e:
            failed_ids.append(pdb_id)
        print()
    
    print(f"\n{'='*50}")
    print(f"下载完成: 成功 {success_count}/{len(pdb_ids)}")
    if failed_ids:
        print(f"失败列表: {', '.join(failed_ids)}")


def main():
    """主函数"""
    # 要下载的PDB ID列表 (基因对应的真实PDB结构)
    # DBT -> 2COO (二氢硫辛酰胺支链转酰基酶E2, E3结合域)
    # GCSH -> 2EDG (甘氨酸裂解系统H蛋白)
    # TYR -> 4P6R (酪氨酸酶)
    pdb_ids = ["2COO", "2EDG", "4P6R"]
    output_dir = "."
    
    print("="*50)
    print("PDB文件批量下载工具")
    print("="*50)
    print(f"目标PDB: {', '.join(pdb_ids)}")
    print()
    
    # 批量下载PDB格式
    download_pdb_batch(pdb_ids, output_dir, format="pdb")
    
    print("\n" + "="*50)
    print("所有下载任务完成!")
    print("="*50)


if __name__ == "__main__":
    main()
