#!/usr/bin/env python3
"""
文献下载脚本 - 尝试多个Sci-Hub镜像下载论文
"""

import requests
import time
import os

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
OUTPUT_FILE = "vitamin_c_ferroaging.pdf"

# Sci-Hub镜像列表
SCIHUB_MIRRORS = [
    "https://sci-hub.se",
    "https://sci-hub.st",
    "https://sci-hub.ru",
    "https://sci-hub.wf",
    "https://sci-hub.ren",
]

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
}

def download_from_scihub():
    """尝试从Sci-Hub下载论文"""
    
    for mirror in SCIHUB_MIRRORS:
        try:
            print(f"\n尝试从 {mirror} 下载...")
            
            # 先尝试DOI
            url = f"{mirror}/{DOI}"
            print(f"访问: {url}")
            
            response = requests.get(url, headers=HEADERS, timeout=30, allow_redirects=True)
            print(f"状态码: {response.status_code}")
            
            if response.status_code == 200:
                # 查找PDF链接
                content = response.text
                
                # 查找嵌入的PDF
                if '.pdf' in content:
                    print("找到PDF链接!")
                    
                    # 尝试提取PDF URL
                    import re
                    
                    # 查找各种可能的PDF链接格式
                    patterns = [
                        r'location\.href=\'([^\']+\.pdf)\'',
                        r'href=\'([^\']+\.pdf)\'',
                        r'src=\'([^\']+\.pdf)\'',
                        r'"([^"]+\.pdf)"',
                        r'https://[^\s"<>]+\.pdf',
                    ]
                    
                    pdf_url = None
                    for pattern in patterns:
                        matches = re.findall(pattern, content)
                        if matches:
                            pdf_url = matches[0]
                            print(f"找到PDF URL: {pdf_url}")
                            break
                    
                    if pdf_url:
                        # 处理相对URL
                        if pdf_url.startswith('//'):
                            pdf_url = 'https:' + pdf_url
                        elif pdf_url.startswith('/'):
                            pdf_url = mirror + pdf_url
                        elif not pdf_url.startswith('http'):
                            pdf_url = mirror + '/' + pdf_url
                        
                        print(f"下载PDF: {pdf_url}")
                        
                        # 下载PDF
                        pdf_response = requests.get(pdf_url, headers=HEADERS, timeout=60)
                        
                        if pdf_response.status_code == 200 and len(pdf_response.content) > 1000:
                            with open(OUTPUT_FILE, 'wb') as f:
                                f.write(pdf_response.content)
                            print(f"✓ 成功下载! 文件大小: {len(pdf_response.content)} bytes")
                            print(f"✓ 保存为: {OUTPUT_FILE}")
                            return True
                        else:
                            print(f"✗ PDF下载失败或文件太小")
                
                # 如果没有找到PDF链接，保存HTML查看
                else:
                    print("未找到PDF链接，保存HTML查看...")
                    with open('debug_page.html', 'w', encoding='utf-8') as f:
                        f.write(content)
                    print("已保存 debug_page.html")
            else:
                print(f"✗ 访问失败: {response.status_code}")
                
        except Exception as e:
            print(f"✗ 错误: {e}")
        
        time.sleep(2)  # 避免请求过快
    
    return False

def try_unpaywall():
    """尝试使用Unpaywall API获取开放获取版本"""
    try:
        print("\n尝试Unpaywall API...")
        url = f"https://api.unpaywall.org/v2/{DOI}?email=user@example.com"
        response = requests.get(url, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            print(f"Unpaywall响应: {data}")
            
            if data.get('is_oa') and data.get('best_oa_location'):
                pdf_url = data['best_oa_location'].get('url_for_pdf')
                if pdf_url:
                    print(f"找到开放获取PDF: {pdf_url}")
                    pdf_response = requests.get(pdf_url, headers=HEADERS, timeout=60)
                    if pdf_response.status_code == 200:
                        with open(OUTPUT_FILE, 'wb') as f:
                            f.write(pdf_response.content)
                        print(f"✓ 成功下载开放获取版本!")
                        return True
    except Exception as e:
        print(f"✗ Unpaywall错误: {e}")
    
    return False

def main():
    print("="*60)
    print("文献下载工具")
    print("="*60)
    print(f"DOI: {DOI}")
    print(f"PMID: {PMID}")
    print(f"论文: Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
    print("="*60)
    
    # 先尝试Unpaywall（合法途径）
    if try_unpaywall():
        return
    
    # 再尝试Sci-Hub
    if download_from_scihub():
        return
    
    print("\n" + "="*60)
    print("所有下载尝试均失败")
    print("="*60)
    print("\n建议:")
    print("1. 通过机构图书馆访问Cell Metabolism官网")
    print("2. 联系作者请求全文: liuguanghui@ioz.ac.cn")
    print("3. 使用ResearchGate向作者请求")
    print("4. 访问: https://pubmed.ncbi.nlm.nih.gov/41819088/")

if __name__ == "__main__":
    main()
