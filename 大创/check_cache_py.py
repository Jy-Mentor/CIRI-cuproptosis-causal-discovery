import pyreadr
cache = pyreadr.read_r('C:/Users/Jy-Mentor-7/Desktop/你信我 这不是重蹈覆辙/大创/analysis_cache.rds')
print('Keys:', list(cache.keys()))
for k in cache.keys():
    df = cache[k]
    if hasattr(df, 'shape'):
        print(k, df.shape)
    else:
        print(k, 'len:', len(df) if hasattr(df, '__len__') else type(df))
        if hasattr(df, 'head'):
            print('  First few values:', df.head())
