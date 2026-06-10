"""
创建人-小鼠同源映射本地缓存 (v2)
使用 mygene 批量查询 ENSMUSG → mouse symbol + 人类基因同名匹配
与 toxirna_feature_extractor.py 已验证的策略一致
"""
import os, sys, time
import pandas as pd
import mygene

FPKM_DIR = r"D:\反向网络药理学\GAT拓展维度\Toxi\rna_fpkm"
CACHE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "ml_output_v4", "mouse_to_human_orthologs.csv")

def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)

def main():
    os.makedirs(os.path.dirname(CACHE_PATH), exist_ok=True)

    # 读取一个 FPKM 文件获取所有 ENSMUSG ID
    tsv_files = [f for f in os.listdir(FPKM_DIR) if f.endswith('.tsv')]
    if not tsv_files:
        log(f"ERROR: 在 {FPKM_DIR} 中未找到 FPKM 文件")
        sys.exit(1)

    sample_file = os.path.join(FPKM_DIR, tsv_files[0])
    log(f"读取: {sample_file}")
    df = pd.read_csv(sample_file, sep='\t')
    id_col = df.columns[0]
    all_ids = [str(x).strip() for x in df[id_col].values]
    all_ids = [x for x in all_ids if x.startswith("ENSMUSG")]
    # 去掉版本号后缀 (.4, .15 等)
    all_ids = sorted(set(x.split('.')[0] for x in all_ids))
    log(f"共 {len(all_ids)} 个唯一 ENSMUSG ID（去版本后），使用 mygene 批量查询...")

    mg = mygene.MyGeneInfo()
    batch_size = 500
    ensmusg_to_symbol = {}
    total = len(all_ids)
    retries = 3

    for i in range(0, total, batch_size):
        batch = all_ids[i:i+batch_size]
        for attempt in range(retries):
            try:
                results = mg.querymany(batch, scopes='ensembl.gene',
                                       fields='symbol', species='mouse',
                                       returnall=True, as_dataframe=True)
                if 'out' in results and results['out'] is not None:
                    outdf = results['out']
                    for _, row in outdf.iterrows():
                        eid = str(row.get('query', ''))
                        if 'symbol' in row and isinstance(row['symbol'], str):
                            ensmusg_to_symbol[eid] = row['symbol'].upper()
                    log(f"  {i+len(batch)}/{total} → {len(ensmusg_to_symbol)} mapped")
                    break
            except Exception as e:
                log(f"  attempt {attempt+1}/{retries} 失败: {e}")
                time.sleep(5)
        # 每个 batch 间短暂暂停
        time.sleep(1)

    log(f"mygene 查询完成: {len(ensmusg_to_symbol)}/{total}")

    # 读取已知人类基因列表用于匹配
    script_dir = os.path.dirname(os.path.abspath(__file__))
    human_genes = set()
    for fname in [
        os.path.join(script_dir, "..", "..", "Desktop/GAT/drug_targets.txt"),
        os.path.join(script_dir, "..", "..", "Desktop/GAT/disease_genes.txt"),
        os.path.join(script_dir, "toxirna_enhanced_features.csv"),
    ]:
        try:
            if not os.path.exists(fname):
                continue
            if fname.endswith('.csv'):
                tdf = pd.read_csv(fname)
                if 'gene_symbol' in tdf.columns:
                    for g in tdf['gene_symbol'].values:
                        if isinstance(g, str) and g.strip():
                            human_genes.add(g.strip().upper())
            else:
                with open(fname) as f:
                    for line in f:
                        g = line.strip().upper()
                        if g:
                            human_genes.add(g)
        except Exception as e:
            log(f"  读取 {fname} 失败: {e}")

    log(f"人类基因参考集: {len(human_genes)} 个")

    # 同名匹配 + 保存缓存
    cache_rows = []
    matched = 0
    for eid, msym in ensmusg_to_symbol.items():
        hsym = msym if msym in human_genes else ''
        if hsym:
            matched += 1
        cache_rows.append({
            'ensmusg': eid,
            'mouse_symbol': msym,
            'human_symbol': hsym
        })

    cache_df = pd.DataFrame(cache_rows)
    cache_df.to_csv(CACHE_PATH, index=False, encoding='utf-8-sig')
    log(f"缓存已保存: {CACHE_PATH}")
    log(f"  ENSMUSG 映射: {len(cache_df)}")
    log(f"  匹配人类基因: {matched}")
    log(f"  未匹配: {len(cache_df) - matched}")

if __name__ == "__main__":
    main()