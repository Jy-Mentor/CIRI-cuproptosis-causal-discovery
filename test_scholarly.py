#!/usr/bin/env python3
"""
使用scholarly搜索论文
"""

from scholarly import scholarly

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
TITLE = "Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates"

def main():
    print("="*70)
    print("使用 scholarly 搜索")
    print("="*70)
    print(f"论文: {TITLE}")
    print(f"DOI: {DOI}")
    print("="*70)
    
    try:
        # 搜索论文
        print(f"\n搜索论文...")
        search_query = scholarly.search_pubs(TITLE)
        
        # 获取第一个结果
        paper = next(search_query)
        
        print(f"\n找到论文:")
        print(f"标题: {paper.get('bib', {}).get('title', 'N/A')}")
        print(f"作者: {', '.join(paper.get('bib', {}).get('author', []))}")
        print(f"年份: {paper.get('bib', {}).get('pub_year', 'N/A')}")
        print(f"DOI: {paper.get('bib', {}).get('doi', 'N/A')}")
        
        # 尝试获取PDF链接
        print(f"\n论文信息:")
        for key, value in paper.items():
            if value and key != 'bib':
                print(f"  {key}: {value}")
        
        # 检查是否有epub或PDF
        if 'epub' in paper:
            print(f"\n✓ 找到epub链接: {paper['epub']}")
        if 'pdf' in paper:
            print(f"\n✓ 找到PDF链接: {paper['pdf']}")
            
    except StopIteration:
        print("✗ 未找到论文")
    except Exception as e:
        print(f"✗ 错误: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
