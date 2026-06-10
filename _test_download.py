import requests, sys

# Test miRTarBase URLs
urls = [
    ('miRTarBase 下载页', 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/'),
    ('hsa_MTI.xlsx', 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/downloads/hsa_MTI.xlsx'),
    ('MTI.xls', 'https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/MTI.xls'),
]

session = requests.Session()
session.headers.update({'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'})

for name, url in urls:
    try:
        resp = session.head(url, timeout=30, allow_redirects=True)
        cl = resp.headers.get('Content-Length', '0')
        if cl and cl != '0':
            print(f'[OK] {name}: HTTP {resp.status_code}, Size: {int(cl)/1024/1024:.1f} MB')
        else:
            print(f'[OK] {name}: HTTP {resp.status_code}')
    except Exception as e:
        print(f'[FAIL] {name}: {type(e).__name__}: {e}')