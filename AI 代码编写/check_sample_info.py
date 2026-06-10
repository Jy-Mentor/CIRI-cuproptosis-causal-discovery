# 查看系列矩阵文件的样本信息
with open(r'C:\Users\Jy-Mentor-7\Downloads\GSE16561_series_matrix (1).txt', 'r') as f:
    lines = []
    for i, line in enumerate(f):
        lines.append(line)
        if i >= 500:
            break
    
    # 打印包含样本标题的行
    print("样本标题:")
    for line in lines:
        if '!Sample_title' in line:
            print(line.strip())
    
    # 打印包含样本ID的行
    print("\n样本ID:")
    for line in lines:
        if '!Sample_geo_accession' in line:
            print(line.strip())