#!/usr/bin/env python3
"""
文献下载脚本 V2 - 使用多种方法尝试下载
"""

import requests
import time
import re

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
OUTPUT_FILE = "vitamin_c_ferroaging.pdf"

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.5',
}

def try_direct_sci_hub():
    """直接尝试Sci-Hub的PDF链接格式"""
    
    # Sci-Hub常用的PDF存储域名
    pdf_domains = [
        "https://sci-hub.se",
        "https://sci-hub.st",
        "https://sci-hub.ru", 
        "https://sci-hub.wf",
        "https://sci-hub.mksa.top",
        "https://sci-hub.yncjkj.com",
    ]
    
    for domain in pdf_domains:
        try:
            print(f"\n尝试 {domain}...")
            
            # 访问论文页面
            url = f"{domain}/{DOI}"
            session = requests.Session()
            response = session.get(url, headers=HEADERS, timeout=30)
            
            print(f"页面状态码: {response.status_code}")
            
            if response.status_code == 200:
                content = response.text
                
                # 查找PDF链接 - 改进匹配模式
                pdf_patterns = [
                    r'<iframe[^>]+src=["\']([^"\']+\.pdf)["\']',
                    r'<embed[^>]+src=["\']([^"\']+\.pdf)["\']',
                    r'location\.href=["\']([^"\']+\.pdf)["\']',
                    r'location=["\']([^"\']+\.pdf)["\']',
                    r'href=["\']([^"\']+\.pdf)["\']',
                    r'src=["\']([^"\']+\.pdf)["\']',
                    r'(https?://[^\s"<>]+\.pdf)',
                    r'(//[^\s"<>]+\.pdf)',
                ]
                
                pdf_urls = []
                for pattern in pdf_patterns:
                    matches = re.findall(pattern, content, re.IGNORECASE)
                    pdf_urls.extend(matches)
                
                # 去重
                pdf_urls = list(set(pdf_urls))
                print(f"找到 {len(pdf_urls)} 个可能的PDF链接")
                
                for pdf_url in pdf_urls:
                    print(f"  尝试: {pdf_url[:80]}...")
                    
                    # 处理URL
                    if pdf_url.startswith('//'):
                        pdf_url = 'https:' + pdf_url
                    elif pdf_url.startswith('/'):
                        pdf_url = domain + pdf_url
                    elif not pdf_url.startswith('http'):
                        pdf_url = domain + '/' + pdf_url
                    
                    try:
                        pdf_response = session.get(pdf_url, headers=HEADERS, timeout=60, stream=True)
                        
                        if pdf_response.status_code == 200:
                            content_type = pdf_response.headers.get('content-type', '')
                            print(f"  Content-Type: {content_type}")
                            
                            # 检查是否是PDF
                            if 'pdf' in content_type.lower() or pdf_url.endswith('.pdf'):
                                # 下载内容
                                content_bytes = b''
                                for chunk in pdf_response.iter_content(chunk_size=8192):
                                    if chunk:
                                        content_bytes += chunk
                                
                                if len(content_bytes) > 10000:  # 至少10KB
                                    with open(OUTPUT_FILE, 'wb') as f:
                                        f.write(content_bytes)
                                    print(f"✓✓✓ 成功下载! 文件大小: {len(content_bytes)} bytes")
                                    return True
                    except Exception as e:
                        print(f"  错误: {e}")
                        continue
                        
        except Exception as e:
            print(f"✗ 错误: {e}")
        
        time.sleep(1)
    
    return False

def try_alternative_sources():
    """尝试其他来源"""
    
    # 尝试LibGen
    try:
        print("\n尝试 Library Genesis...")
        # LibGen通常需要搜索，这里简化处理
        pass
    except:
        pass
    
    # 尝试通过PubMed Central (PMC)
    try:
        print("\n尝试 PubMed Central...")
        pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMID{PMID}/"
        response = requests.get(pmc_url, headers=HEADERS, timeout=30)
        
        if response.status_code == 200 and 'pdf' in response.text.lower():
            print("PMC可能有全文，请手动访问:")
            print(pmc_url)
    except:
        pass
    
    return False

def main():
    print("="*70)
    print("文献下载工具 V2")
    print("="*70)
    print(f"论文: Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
    print(f"DOI: {DOI}")
    print(f"PMID: {PMID}")
    print("="*70)
    
    # 尝试Sci-Hub
    if try_direct_sci_hub():
        print("\n" + "="*70)
        print("✓ 下载成功!")
        print(f"✓ 文件保存为: {OUTPUT_FILE}")
        print("="*70)
        return
    
    # 尝试其他来源
    try_alternative_sources()
    
    print("\n" + "="*70)
    print("自动下载失败")
    print("="*70)
    print("\n替代方案:")
    print("1. 手动访问Sci-Hub:")
    print(f"   https://sci-hub.se/{DOI}")
    print(f"   https://sci-hub.st/{DOI}")
    print(f"   https://sci-hub.wf/{DOI}")
    print("\n2. 联系作者:")
    print("   刘光慧研究员: liuguanghui@ioz.ac.cn")
    print("   中国科学院动物研究所")
    print("\n3. 机构访问:")
    print("   https://www.cell.com/cell-metabolism/fulltext/S1550-4131(26)00053-7")
    print("\n4. PubMed:")
    print("   https://pubmed.ncbi.nlm.nih.gov/41819088/")

if __name__ == "__main__":
    main()
