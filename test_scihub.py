#!/usr/bin/env python3
"""
使用scihub包下载论文
"""

from scihub import SciHub

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
OUTPUT_FILE = "vitamin_c_ferroaging.pdf"

def main():
    print("="*70)
    print("使用 scihub Python包下载")
    print("="*70)
    print(f"论文: Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
    print(f"DOI: {DOI}")
    print("="*70)
    
    try:
        # 初始化SciHub
        sh = SciHub()
        
        # 使用DOI下载
        print(f"\n尝试下载: {DOI}")
        result = sh.download(DOI, OUTPUT_FILE)
        
        if result:
            print(f"✓ 下载成功!")
            print(f"✓ 保存为: {OUTPUT_FILE}")
        else:
            print("✗ 下载失败")
            
    except Exception as e:
        print(f"✗ 错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
