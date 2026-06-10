import requests, re, sys

session = requests.Session()
session.headers.update({
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.9',
})

request_kwargs = {
    'timeout': 30,
    'verify': False,
}

import urllib3
urllib3.disable_warnings()

url = 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/'
try:
    resp = session.get(url, **request_kwargs)
    print(f'Status: {resp.status_code}')
    print(f'Len: {len(resp.text)}')
    
    links = re.findall(r'href=[\"\\\']([^\"\\\']*.(?:xlsx?|csv|txt|zip))[\"\\\']', resp.text, re.I)
    print(f'Download links: {links}')
    
    # Look for download section
    if 'download' in resp.text.lower():
        print('Page contains download section')
    
    # Try to find the download URL
    for line in resp.text.split('\n'):
        if '.xls' in line.lower() or 'download' in line.lower():
            print(f'  -> {line.strip()[:200]}')
except Exception as e:
    print(f'Error: {type(e).__name__}: {e}')
    print(f'Content: {getattr(e, "response", None)}')