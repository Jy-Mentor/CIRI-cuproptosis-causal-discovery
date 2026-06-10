import pandas as pd
f = r'c:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\results\stage5_ppi_mcode\string_ppi.tsv'
df = pd.read_csv(f, sep='\t')
print(f'Score range: [{df["score"].min()}, {df["score"].max()}]')
print(f'>=400: {len(df[df["score"]>=400])}')
print(f'>=0.4: {len(df[df["score"]>=0.4])}')
