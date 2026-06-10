import pandas as pd
import time
import random
from biomart import BiomartServer

# 读取大鼠基因文件
data = pd.read_csv('C:\\Users\\Jy-Mentor-7\\Desktop\\新建 文本文档.txt', sep='→', engine='python')
rats_genes = data['GENE_SYMBOL'].tolist()
rats_genes = [gene for gene in rats_genes if gene != 'GENE_SYMBOL']  # 移除表头

# 连接到Ensembl Biomart服务器
server = BiomartServer('http://www.ensembl.org/biomart')

# 选择rat基因数据集
rat_dataset = server.datasets['rnorvegicus_gene_ensembl']

# 定义查询参数
attributes = ['external_gene_name', 'hsapiens_homolog_associated_gene_name']

# 分批处理基因列表，每批50个基因
batch_size = 50
homology_map = {}

for i in range(0, len(rats_genes), batch_size):
    batch_genes = rats_genes[i:i+batch_size]
    batch_num = i//batch_size + 1
    total_batches = (len(rats_genes) + batch_size - 1)//batch_size
    print(f"处理批次 {batch_num}/{total_batches}")
    
    # 重试逻辑
    max_retries = 3
    retry_count = 0
    success = False
    
    while retry_count < max_retries and not success:
        try:
            # 执行查询
            response = rat_dataset.search({
                'filters': {'external_gene_name': batch_genes},
                'attributes': attributes
            })
            
            # 解析结果
            for line in response.iter_lines():
                line = line.decode('utf-8')
                if line:
                    parts = line.split('\t')
                    if len(parts) == 2 and parts[1]:
                        homology_map[parts[0]] = parts[1]
            
            success = True
            print(f"  批次 {batch_num} 处理成功")
            
        except Exception as e:
            retry_count += 1
            print(f"  批次 {batch_num} 处理失败，尝试重试 {retry_count}/{max_retries}")
            print(f"  错误: {str(e)}")
            if retry_count < max_retries:
                # 随机延迟2-5秒后重试
                delay = random.uniform(2, 5)
                print(f"  等待 {delay:.1f} 秒后重试...")
                time.sleep(delay)
            else:
                print(f"  批次 {batch_num} 达到最大重试次数，跳过")
    
    # 每批处理后添加随机延迟，避免服务器拒绝请求
    if batch_num < total_batches:
        delay = random.uniform(1, 3)
        time.sleep(delay)

# 创建结果DataFrame
result_df = pd.DataFrame({
    'Rat_Gene': list(homology_map.keys()),
    'Human_Homolog': list(homology_map.values())
})

# 保存结果到CSV文件
output_path = 'C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\AI 代码编写\\rat_to_human_homology.csv'
result_df.to_csv(output_path, index=False, encoding='utf-8-sig')

print(f"同源映射完成！结果已保存到: {output_path}")
print(f"成功映射的基因数量: {len(result_df)}")
print(f"未能映射的基因数量: {len(rats_genes) - len(result_df)}")
