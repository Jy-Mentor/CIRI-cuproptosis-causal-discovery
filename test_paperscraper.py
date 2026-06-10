#!/usr/bin/env python3
"""
使用paperscraper下载论文
"""

import asyncio
from paperscraper import PaperScraper

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
OUTPUT_FILE = "vitamin_c_ferroaging.pdf"

async def main():
    print("="*70)
    print("使用 paperscraper 下载")
    print("="*70)
    print(f"论文: Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
    print(f"DOI: {DOI}")
    print("="*70)
    
    try:
        # 初始化PaperScraper
        scraper = PaperScraper()
        
        # 使用DOI下载
        print(f"\n尝试下载: {DOI}")
        
        # 搜索论文
        papers = await scraper.search_papers(DOI, limit=1)
        
        if papers:
            print(f"找到论文: {papers[0]}")
            
            # 尝试下载PDF
            pdf_path = await scraper.download_pdf(papers[0], OUTPUT_FILE)
            
            if pdf_path:
                print(f"✓ 下载成功!")
                print(f"✓ 保存为: {pdf_path}")
            else:
                print("✗ 下载失败")
        else:
            print("✗ 未找到论文")
            
    except Exception as e:
        print(f"✗ 错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    asyncio.run(main())
