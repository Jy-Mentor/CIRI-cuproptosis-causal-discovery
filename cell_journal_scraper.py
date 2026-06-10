#!/usr/bin/env python3
"""
Cell Journal Scraper - 爬取Cell期刊论文
支持: Cell, Cell Metabolism, Cell Stem Cell等Cell Press期刊
"""

import time
import json
import csv
from datetime import datetime
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import requests

class CellJournalScraper:
    """Cell期刊爬虫类"""
    
    def __init__(self, headless=True):
        """初始化爬虫"""
        self.base_url = "https://www.cell.com"
        self.journals = {
            'cell': '/cell/home',
            'cell-metabolism': '/cell-metabolism/home',
            'cell-stem-cell': '/cell-stem-cell/home',
            'cell-research': '/cell-research/home',
            'molecular-cell': '/molecular-cell/home',
        }
        
        # 设置Chrome选项
        chrome_options = Options()
        if headless:
            chrome_options.add_argument('--headless')
        chrome_options.add_argument('--no-sandbox')
        chrome_options.add_argument('--disable-dev-shm-usage')
        chrome_options.add_argument('--disable-gpu')
        chrome_options.add_argument('--window-size=1920,1080')
        chrome_options.add_argument('--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36')
        
        self.driver = webdriver.Chrome(options=chrome_options)
        self.wait = WebDriverWait(self.driver, 10)
        
    def search_papers(self, journal='cell-metabolism', query='', max_results=10):
        """
        搜索论文
        
        Args:
            journal: 期刊名称 (cell, cell-metabolism, cell-stem-cell等)
            query: 搜索关键词
            max_results: 最大结果数
        """
        papers = []
        
        try:
            # 构建搜索URL
            if journal in self.journals:
                search_url = f"{self.base_url}/action/doSearch?journalCode={journal.replace('-', '')}&searchText={query}"
            else:
                search_url = f"{self.base_url}/action/doSearch?searchText={query}"
            
            print(f"搜索URL: {search_url}")
            self.driver.get(search_url)
            
            # 等待页面加载
            time.sleep(3)
            
            # 获取论文列表
            articles = self.driver.find_elements(By.CSS_SELECTOR, '.search-result-item, .article-item')
            
            print(f"找到 {len(articles)} 篇文章")
            
            for i, article in enumerate(articles[:max_results]):
                try:
                    paper = self._extract_paper_info(article)
                    if paper:
                        papers.append(paper)
                        print(f"[{i+1}] {paper.get('title', 'N/A')[:80]}...")
                except Exception as e:
                    print(f"提取第{i+1}篇文章时出错: {e}")
                    continue
                    
        except Exception as e:
            print(f"搜索出错: {e}")
        
        return papers
    
    def _extract_paper_info(self, article_element):
        """提取单篇文章信息"""
        paper = {}
        
        try:
            # 标题
            title_elem = article_element.find_element(By.CSS_SELECTOR, '.title, .article-title, h3 a')
            paper['title'] = title_elem.text.strip()
            paper['url'] = title_elem.get_attribute('href')
        except:
            paper['title'] = 'N/A'
            paper['url'] = ''
        
        try:
            # 作者
            authors_elem = article_element.find_element(By.CSS_SELECTOR, '.authors, .author-list')
            paper['authors'] = authors_elem.text.strip()
        except:
            paper['authors'] = 'N/A'
        
        try:
            # 日期
            date_elem = article_element.find_element(By.CSS_SELECTOR, '.date, .published-date')
            paper['date'] = date_elem.text.strip()
        except:
            paper['date'] = 'N/A'
        
        try:
            # DOI
            doi_elem = article_element.find_element(By.CSS_SELECTOR, '[data-doi], .doi')
            paper['doi'] = doi_elem.get_attribute('data-doi') or doi_elem.text.strip()
        except:
            paper['doi'] = 'N/A'
        
        try:
            # 摘要
            abstract_elem = article_element.find_element(By.CSS_SELECTOR, '.abstract, .summary')
            paper['abstract'] = abstract_elem.text.strip()
        except:
            paper['abstract'] = 'N/A'
        
        return paper
    
    def get_paper_details(self, url):
        """获取论文详细信息"""
        try:
            self.driver.get(url)
            time.sleep(2)
            
            details = {}
            
            # 获取完整摘要
            try:
                abstract = self.driver.find_element(By.CSS_SELECTOR, '.abstract-content, #abstract')
                details['full_abstract'] = abstract.text.strip()
            except:
                details['full_abstract'] = 'N/A'
            
            # 获取关键词
            try:
                keywords = self.driver.find_element(By.CSS_SELECTOR, '.keywords')
                details['keywords'] = keywords.text.strip()
            except:
                details['keywords'] = 'N/A'
            
            # 获取PDF链接
            try:
                pdf_link = self.driver.find_element(By.CSS_SELECTOR, 'a[href*=".pdf"], .pdf-link')
                details['pdf_url'] = pdf_link.get_attribute('href')
            except:
                details['pdf_url'] = 'N/A'
            
            return details
            
        except Exception as e:
            print(f"获取详情出错: {e}")
            return {}
    
    def download_pdf(self, pdf_url, filename=None):
        """
        下载PDF（需要机构权限）
        """
        if not filename:
            filename = f"paper_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
        
        try:
            print(f"尝试下载PDF: {pdf_url}")
            
            # 使用selenium下载
            self.driver.get(pdf_url)
            time.sleep(5)  # 等待下载
            
            print(f"✓ PDF下载完成: {filename}")
            return True
            
        except Exception as e:
            print(f"✗ PDF下载失败: {e}")
            return False
    
    def save_to_csv(self, papers, filename='cell_papers.csv'):
        """保存到CSV"""
        if not papers:
            print("没有论文可保存")
            return
        
        keys = papers[0].keys()
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=keys)
            writer.writeheader()
            writer.writerows(papers)
        
        print(f"✓ 已保存 {len(papers)} 篇论文到 {filename}")
    
    def save_to_json(self, papers, filename='cell_papers.json'):
        """保存到JSON"""
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(papers, f, ensure_ascii=False, indent=2)
        
        print(f"✓ 已保存 {len(papers)} 篇论文到 {filename}")
    
    def close(self):
        """关闭浏览器"""
        self.driver.quit()
        print("✓ 浏览器已关闭")


def main():
    """主函数"""
    print("="*70)
    print("Cell Journal Scraper")
    print("Cell期刊论文爬虫工具")
    print("="*70)
    
    # 创建爬虫实例
    scraper = CellJournalScraper(headless=False)  # 设置为True可隐藏浏览器
    
    try:
        # 搜索Vitamin C论文
        print("\n搜索论文: Vitamin ACSL4 ferro-aging")
        papers = scraper.search_papers(
            journal='cell-metabolism',
            query='Vitamin ACSL4 ferro-aging',
            max_results=5
        )
        
        if papers:
            print(f"\n✓ 找到 {len(papers)} 篇论文")
            
            # 保存结果
            scraper.save_to_csv(papers, 'cell_papers.csv')
            scraper.save_to_json(papers, 'cell_papers.json')
            
            # 获取第一篇论文的详细信息
            if papers[0].get('url'):
                print("\n获取第一篇论文详情...")
                details = scraper.get_paper_details(papers[0]['url'])
                print(f"完整摘要: {details.get('full_abstract', 'N/A')[:200]}...")
                print(f"PDF链接: {details.get('pdf_url', 'N/A')}")
        else:
            print("✗ 未找到论文")
    
    except Exception as e:
        print(f"错误: {e}")
    
    finally:
        scraper.close()


if __name__ == "__main__":
    main()
