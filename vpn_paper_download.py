#!/usr/bin/env python3
"""
通过 VPN 下载 Cell 论文
使用 Selenium 保持 VPN 会话
"""

import time
import os
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.edge.options import Options as EdgeOptions
from selenium.webdriver.edge.service import Service as EdgeService

class VPNDownloader:
    def __init__(self):
        # 设置 Edge 选项
        edge_options = EdgeOptions()
        # 设置下载目录
        self.download_dir = os.path.join(os.getcwd(), "downloads")
        os.makedirs(self.download_dir, exist_ok=True)
        
        prefs = {
            "download.default_directory": self.download_dir,
            "download.prompt_for_download": False,
            "download.directory_upgrade": True,
            "plugins.always_open_pdf_externally": True,  # 直接下载PDF而不是打开
        }
        edge_options.add_experimental_option("prefs", prefs)
        
        # 其他选项
        edge_options.add_argument('--no-sandbox')
        edge_options.add_argument('--disable-dev-shm-usage')
        edge_options.add_argument('--disable-gpu')
        edge_options.add_argument('--window-size=1920,1080')
        
        self.driver = webdriver.Edge(options=edge_options)
        self.wait = WebDriverWait(self.driver, 20)
        
    def login_vpn(self):
        """登录 VPN"""
        print("正在登录 VPN...")
        self.driver.get("https://webvpn.gxtcmu.edu.cn/")
        
        # 等待用户手动登录
        print("请在浏览器中登录 VPN")
        print("登录成功后，按 Enter 继续...")
        input()
        
    def download_via_webofscience(self):
        """通过 Web of Science 下载"""
        try:
            # 访问 Web of Science
            print("访问 Web of Science...")
            self.driver.get("https://www.webofscience.com/")
            time.sleep(3)
            
            # 搜索论文
            search_box = self.wait.until(
                EC.presence_of_element_located((By.CSS_SELECTOR, "input[data-ta='search-box-input']"))
            )
            search_box.send_keys("Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
            
            # 点击搜索
            search_btn = self.driver.find_element(By.CSS_SELECTOR, "button[data-ta='search-button']")
            search_btn.click()
            time.sleep(3)
            
            # 找到论文并点击
            paper_link = self.wait.until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, "a[data-ta='title-link']"))
            )
            paper_link.click()
            time.sleep(2)
            
            # 查找全文链接
            full_text_btn = self.driver.find_element(By.CSS_SELECTOR, "button[data-ta='full-text-link']")
            full_text_btn.click()
            time.sleep(5)
            
            print("✓ 已跳转到出版商网站，请手动下载 PDF")
            print("下载完成后，按 Enter 继续...")
            input()
            
        except Exception as e:
            print(f"错误: {e}")
            
    def download_direct(self):
        """直接访问 Cell 官网下载"""
        try:
            # 方法1: 通过 VPN 代理访问 Cell 官网
            # 先访问 Cell 主站
            print("访问 Cell 官网...")
            self.driver.get("https://www.cell.com")
            time.sleep(3)
            
            # 搜索论文
            print("搜索论文...")
            search_box = self.wait.until(
                EC.presence_of_element_located((By.CSS_SELECTOR, "input[type='search'], #search-input, .search-input"))
            )
            search_box.send_keys("Vitamin C inhibits ACSL4 to alleviate ferro-aging in primates")
            
            # 点击搜索按钮
            search_btn = self.driver.find_element(By.CSS_SELECTOR, "button[type='submit'], .search-button")
            search_btn.click()
            time.sleep(3)
            
            # 点击第一篇结果
            print("点击搜索结果...")
            first_result = self.wait.until(
                EC.element_to_be_clickable((By.CSS_SELECTOR, ".search-result-item a, .article-title a"))
            )
            first_result.click()
            time.sleep(5)
            
            # 查找下载按钮
            print("查找下载按钮...")
            try:
                # 尝试多种选择器
                selectors = [
                    "a[href*='pdf']",
                    ".pdf-download",
                    "button[data-aa-button='download-pdf']",
                    ".download-pdf",
                    "a[title*='PDF']",
                    "a[title*='Download']"
                ]
                
                for selector in selectors:
                    try:
                        download_btn = self.driver.find_element(By.CSS_SELECTOR, selector)
                        print(f"✓ 找到下载按钮: {selector}")
                        download_btn.click()
                        print("✓ 已点击下载")
                        time.sleep(10)
                        return
                    except:
                        continue
                
                print("✗ 未找到下载按钮")
                print("当前页面标题:", self.driver.title)
                print("当前URL:", self.driver.current_url)
                print("\n请手动在浏览器中查找下载按钮")
                print("按 Enter 继续...")
                input()
                
            except Exception as e:
                print(f"查找下载按钮出错: {e}")
                print("请手动下载")
                input()
                
        except Exception as e:
            print(f"错误: {e}")
            import traceback
            traceback.print_exc()
            
    def check_download(self):
        """检查下载的文件"""
        files = os.listdir(self.download_dir)
        pdf_files = [f for f in files if f.endswith('.pdf')]
        
        if pdf_files:
            print(f"\n✓ 下载成功!")
            for f in pdf_files:
                filepath = os.path.join(self.download_dir, f)
                size = os.path.getsize(filepath)
                print(f"  {f} ({size/1024:.1f} KB)")
        else:
            print("\n✗ 未找到下载的 PDF 文件")
            
    def close(self):
        """关闭浏览器"""
        self.driver.quit()
        print("✓ 浏览器已关闭")

def main():
    print("="*70)
    print("VPN 论文下载工具")
    print("="*70)
    print("\n这个工具会打开 Edge 浏览器")
    print("你需要手动登录 VPN，然后脚本会尝试下载论文")
    print("="*70)
    
    downloader = VPNDownloader()
    
    try:
        # 登录 VPN
        downloader.login_vpn()
        
        # 尝试直接下载
        print("\n尝试直接访问 Cell 官网...")
        downloader.download_direct()
        
        # 检查下载
        downloader.check_download()
        
    except Exception as e:
        print(f"错误: {e}")
    finally:
        print("\n按 Enter 关闭浏览器...")
        input()
        downloader.close()

if __name__ == "__main__":
    main()
