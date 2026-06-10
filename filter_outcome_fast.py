import pandas as pd 
import time 

# 1. 读取暴露SNP列表（从你的Excel提取后保存） 
snps = set() 
with open("D:/EQTL/mr_results_p5e-06/exposure_snplist.txt", "r") as f: 
    for line in f: 
        snps.add(line.strip()) 

print(f"暴露SNP数量: {len(snps)}") 

# 2. 分块过滤50GB结局文件 
# chunksize=200万行，减少I/O次数 
start = time.time() 
chunks = pd.read_csv( 
    "D:/EQTL/eqtlgen_ieu_outcome.csv", 
    chunksize=2000000, 
    low_memory=False, 
    dtype=str  # 全部读字符串，避免类型推断拖慢 
) 

filtered = [] 
for i, chunk in enumerate(chunks): 
    # 用isin过滤，pandas底层优化过 
    sub = chunk[chunk['SNP'].isin(snps)] 
    if len(sub) > 0: 
        filtered.append(sub) 
    
    if i % 5 == 0: 
        elapsed = time.time() - start 
        print(f"已处理 {i*2}00万行, 耗时 {elapsed:.0f}秒, 匹配到 {sum(len(x) for x in filtered)} 行") 

# 3. 合并保存 
result = pd.concat(filtered, ignore_index=True) 
result.to_csv("D:/EQTL/mr_results_p5e-06/outcome_filtered_fast.csv", index=False) 

print(f"\n完成! 过滤后 {len(result)} 行, 总耗时 {(time.time()-start)/60:.1f} 分钟")
