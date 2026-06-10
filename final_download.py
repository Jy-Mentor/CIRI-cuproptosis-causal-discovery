#!/usr/bin/env python3
"""
最终尝试 - 使用多种方法下载论文
"""

import requests
from bs4 import BeautifulSoup
import time
import re

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
OUTPUT_FILE = "vitamin_c_ferroaging.pdf"

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.9,zh-CN;q=0.8,zh;q=0.7',
}

def try_sci_hub_mirror(mirror_url):
    """尝试从单个Sci-Hub镜像下载"""
    try:
        print(f"\n{'='*60}")
        print(f"尝试: {mirror_url}")
        print('='*60)
        
        session = requests.Session()
        
        # 访问论文页面
        url = f"{mirror_url}/{DOI}"
        print(f"访问: {url}")
        
        response = session.get(url, headers=HEADERS, timeout=30)
        print(f"状态码: {response.status_code}")
        
        if response.status_code != 200:
            print(f"✗ 页面访问失败")
            return False
        
        # 解析HTML
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # 查找PDF链接
        pdf_url = None
        
        # 方法1: 查找iframe
        iframe = soup.find('iframe')
        if iframe and iframe.get('src'):
            pdf_url = iframe['src']
            print(f"✓ 从iframe找到PDF: {pdf_url}")
        
        # 方法2: 查找embed
        if not pdf_url:
            embed = soup.find('embed')
            if embed and embed.get('src'):
                pdf_url = embed['src']
                print(f"✓ 从embed找到PDF: {pdf_url}")
        
        # 方法3: 查找button的onclick
        if not pdf_url:
            button = soup.find('button', onclick=True)
            if button:
                onclick = button['onclick']
                match = re.search(r'location\.href=[\'"]([^\'"]+)[\'"]', onclick)
                if match:
                    pdf_url = match.group(1)
                    print(f"✓ 从button找到PDF: {pdf_url}")
        
        # 方法4: 查找所有链接
        if not pdf_url:
            for link in soup.find_all('a', href=True):
                href = link['href']
                if '.pdf' in href.lower():
                    pdf_url = href
                    print(f"✓ 从链接找到PDF: {pdf_url}")
                    break
        
        # 方法5: 查找div中的文本
        if not pdf_url:
            divs = soup.find_all('div')
            for div in divs:
                text = div.get_text()
                if '.pdf' in text:
                    match = re.search(r'https?://[^\s<>"]+\.pdf', text)
                    if match:
                        pdf_url = match.group(0)
                        print(f"✓ 从div文本找到PDF: {pdf_url}")
                        break
        
        if not pdf_url:
            print("✗ 未找到PDF链接")
            # 保存HTML用于调试
            with open(f'debug_{mirror_url.split("/")[-1]}.html', 'w', encoding='utf-8') as f:
                f.write(soup.prettify())
            print(f"  已保存HTML用于调试")
            return False
        
        # 处理相对URL
        if pdf_url.startswith('//'):
            pdf_url = 'https:' + pdf_url
        elif pdf_url.startswith('/'):
            pdf_url = mirror_url + pdf_url
        elif not pdf_url.startswith('http'):
            pdf_url = mirror_url + '/' + pdf_url
        
        print(f"\n下载PDF: {pdf_url}")
        
        # 下载PDF
        pdf_response = session.get(pdf_url, headers=HEADERS, timeout=60, stream=True)
        
        if pdf_response.status_code == 200:
            content_type = pdf_response.headers.get('content-type', '')
            content_length = pdf_response.headers.get('content-length', 'unknown')
            print(f"Content-Type: {content_type}")
            print(f"Content-Length: {content_length}")
            
            # 检查是否是PDF
            is_pdf = 'pdf' in content_type.lower() or pdf_url.endswith('.pdf')
            
            if is_pdf or content_length == 'unknown' or int(content_length) > 10000:
                # 下载内容
                content = b''
                for chunk in pdf_response.iter_content(chunk_size=8192):
                    if chunk:
                        content += chunk
                
                if len(content) > 10000:  # 至少10KB
                    with open(OUTPUT_FILE, 'wb') as f:
                        f.write(content)
                    print(f"\n✓✓✓ 成功下载!")
                    print(f"✓ 文件大小: {len(content)} bytes ({len(content)/1024:.1f} KB)")
                    print(f"✓ 保存为: {OUTPUT_FILE}")
                    return True
                else:
                    print(f"✗ 文件太小: {len(content)} bytes")
            else:
                print(f"✗ 不是PDF文件")
        else:
            print(f"✗ PDF下载失败: {pdf_response.status_code}")
            
    except Exception as e:
        print(f"✗ 错误: {e}")
    
    return False

def main():
    print("="*70)
    print("文献下载工具 - 最终版")
    print("="*70)
    print(f"论文: Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
    print(f"DOI: {DOI}")
    print(f"PMID: {PMID}")
    print("="*70)
    
    # Sci-Hub镜像列表
    mirrors = [
        "https://sci-hub.se",
        "https://sci-hub.st",
        "https://sci-hub.ru",
        "https://sci-hub.wf",
        "https://sci-hub.mksa.top",
        "https://sci-hub.yncjkj.com",
        "https://sci-hub.hkvisa.net",
        "https://sci-hub.ren",
    ]
    
    for mirror in mirrors:
        if try_sci_hub_mirror(mirror):
            print("\n" + "="*70)
            print("🎉 下载成功!")
            print("="*70)
            return
        time.sleep(2)  # 避免请求过快
    
    print("\n" + "="*70)
    print("❌ 所有Sci-Hub镜像都失败了")
    print("="*70)
    print("\n原因分析:")
    print("1. 这篇论文太新了(2026年3月发表)")
    print("2. Sci-Hub通常需要3-6个月才能收录新论文")
    print("3. Cell Metabolism是付费期刊，Sci-Hub收录需要时间")
    print("\n获取全文的合法途径:")
    print("1. 📧 联系作者: liuguanghui@ioz.ac.cn (刘光慧研究员)")
    print("2. 🏫 通过机构图书馆访问Cell Metabolism官网")
    print("3. 🔬 在ResearchGate上向作者请求全文")
    print("4. ⏰ 等待几个月后Sci-Hub会收录")

if __name__ == "__main__":
    main()
