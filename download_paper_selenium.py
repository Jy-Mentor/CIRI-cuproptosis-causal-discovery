#!/usr/bin/env python3
"""
文献下载脚本 V3 - 使用Selenium模拟浏览器
"""

import time
import os

# 论文信息
DOI = "10.1016/j.cmet.2026.02.010"
PMID = "41819088"
OUTPUT_FILE = "vitamin_c_ferroaging.pdf"

def try_selenium_download():
    """使用Selenium下载"""
    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait
        from selenium.webdriver.support import expected_conditions as EC
        from selenium.webdriver.chrome.options import Options
        
        print("启动Chrome浏览器...")
        
        # 配置Chrome选项
        chrome_options = Options()
        chrome_options.add_argument('--headless')  # 无头模式
        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument('--disable-gpu')
        chrome_options.add_argument('--window-size=1920,1080')
        
        # 设置下载选项
        prefs = {
            "download.default_directory": os.getcwd(),
            "download.prompt_for_download": False,
            "plugins.always_open_pdf_externally": True
        }
        chrome_options.add_experimental_option("prefs", prefs)
        
        # 启动浏览器
        driver = webdriver.Chrome(options=chrome_options)
        
        # Sci-Hub镜像
        sci_hub_urls = [
            f"https://sci-hub.ru/{DOI}",
            f"https://sci-hub.wf/{DOI}",
            f"https://sci-hub.mksa.top/{DOI}",
        ]
        
        for url in sci_hub_urls:
            try:
                print(f"\n访问: {url}")
                driver.get(url)
                
                # 等待页面加载
                time.sleep(5)
                
                # 查找PDF嵌入元素
                print("查找PDF元素...")
                
                # 尝试找到iframe
                try:
                    iframe = driver.find_element(By.TAG_NAME, 'iframe')
                    pdf_url = iframe.get_attribute('src')
                    print(f"找到iframe PDF: {pdf_url}")
                    
                    # 直接访问PDF URL
                    driver.get(pdf_url)
                    time.sleep(5)
                    
                    # 保存页面内容
                    pdf_content = driver.page_source.encode('utf-8')
                    if len(pdf_content) > 10000:
                        with open(OUTPUT_FILE, 'wb') as f:
                            f.write(pdf_content)
                        print(f"✓ 下载成功! 大小: {len(pdf_content)} bytes")
                        driver.quit()
                        return True
                        
                except Exception as e:
                    print(f"iframe方法失败: {e}")
                
                # 尝试找到embed
                try:
                    embed = driver.find_element(By.TAG_NAME, 'embed')
                    pdf_url = embed.get_attribute('src')
                    print(f"找到embed PDF: {pdf_url}")
                    
                    driver.get(pdf_url)
                    time.sleep(5)
                    
                    pdf_content = driver.page_source.encode('utf-8')
                    if len(pdf_content) > 10000:
                        with open(OUTPUT_FILE, 'wb') as f:
                            f.write(pdf_content)
                        print(f"✓ 下载成功! 大小: {len(pdf_content)} bytes")
                        driver.quit()
                        return True
                        
                except Exception as e:
                    print(f"embed方法失败: {e}")
                
                # 尝试找到直接链接
                try:
                    links = driver.find_elements(By.TAG_NAME, 'a')
                    for link in links:
                        href = link.get_attribute('href')
                        if href and '.pdf' in href:
                            print(f"找到PDF链接: {href}")
                            driver.get(href)
                            time.sleep(5)
                            
                            pdf_content = driver.page_source.encode('utf-8')
                            if len(pdf_content) > 10000:
                                with open(OUTPUT_FILE, 'wb') as f:
                                    f.write(pdf_content)
                                print(f"✓ 下载成功! 大小: {len(pdf_content)} bytes")
                                driver.quit()
                                return True
                except Exception as e:
                    print(f"链接方法失败: {e}")
                    
            except Exception as e:
                print(f"访问失败: {e}")
        
        driver.quit()
        
    except ImportError:
        print("未安装selenium，请运行: pip install selenium")
        return False
    except Exception as e:
        print(f"Selenium错误: {e}")
        return False
    
    return False

def main():
    print("="*70)
    print("文献下载工具 V3 - Selenium版本")
    print("="*70)
    print(f"论文: Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
    print(f"DOI: {DOI}")
    print(f"PMID: {PMID}")
    print("="*70)
    
    if try_selenium_download():
        print("\n" + "="*70)
        print("✓ 下载成功!")
        print(f"✓ 文件保存为: {OUTPUT_FILE}")
        print("="*70)
    else:
        print("\n" + "="*70)
        print("Selenium下载也失败了")
        print("="*70)
        print("\n最终建议:")
        print("这篇论文太新了(2026年3月)，Sci-Hub可能还没有收录")
        print("\n获取全文的合法途径:")
        print("1. 通过大学/机构图书馆访问Cell Metabolism")
        print("2. 联系通讯作者刘光慧研究员: liuguanghui@ioz.ac.cn")
        print("3. 在ResearchGate上向作者请求全文")
        print("4. 等待Sci-Hub收录(通常需要几个月)")

if __name__ == "__main__":
    main()
