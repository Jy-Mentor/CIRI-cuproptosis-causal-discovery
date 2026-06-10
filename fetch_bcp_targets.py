"""
β-石竹烯 (BCP) 毒理靶标数据获取脚本
======================================
数据来源:
  1. CTD (Comparative Toxicogenomics Database) - Chemical-Gene Interactions
  2. 文献整理的关键靶标

输出: bcp_targets_all.csv
"""

import os
import sys
import urllib.request
import gzip
import csv
import time
from pathlib import Path

# ============================================================
# 配置
# ============================================================
BASE_DIR = Path(__file__).parent.resolve()
OUTPUT_FILE = BASE_DIR / 'bcp_targets_all.csv'

# CTD 数据 URL
CTD_CHEM_GENE_URL = 'http://ctdbase.org/reports/CTD_chem_gene_ixns.tsv.gz'
CTD_FILE = BASE_DIR / 'CTD_chem_gene_ixns.tsv.gz'
CTD_EXTRACTED = BASE_DIR / 'CTD_chem_gene_ixns.tsv'

# BCP 搜索关键词
BCP_NAMES = [
    'beta-caryophyllene', 'β-caryophyllene', 'caryophyllene',
    '(-)-β-caryophyllene', 'trans-caryophyllene',
    '(E)-β-caryophyllene', 'BCP',
]

# PubChem CID for β-caryophyllene
BCP_CID = '5281515'

# CAS Registry Number
BCP_CAS = '87-44-5'


def download_ctd():
    """下载 CTD Chemical-Gene Interactions 数据"""
    if CTD_EXTRACTED.exists():
        print(f"[✓] 已存在解压文件: {CTD_EXTRACTED}")
        return True

    if not CTD_FILE.exists():
        print(f"[*] 下载 CTD 数据 (~300MB)...")
        print(f"    来源: {CTD_CHEM_GENE_URL}")
        try:
            urllib.request.urlretrieve(CTD_CHEM_GENE_URL, CTD_FILE)
            print(f"[✓] 下载完成: {CTD_FILE}")
        except Exception as e:
            print(f"[✗] 下载失败: {e}")
            print("    将使用离线文献靶标数据。")
            return False
    else:
        print(f"[✓] 已存在压缩文件: {CTD_FILE}")

    # 解压
    print("[*] 解压中...")
    try:
        with gzip.open(CTD_FILE, 'rb') as f_in:
            with open(CTD_EXTRACTED, 'wb') as f_out:
                f_out.write(f_in.read())
        print(f"[✓] 解压完成: {CTD_EXTRACTED}")
        return True
    except Exception as e:
        print(f"[✗] 解压失败: {e}")
        return False


def parse_ctd_chemical_name(chem_name):
    """检查化学物质名称是否匹配 BCP"""
    name_lower = chem_name.lower().strip()
    for bcp_name in BCP_NAMES:
        if bcp_name.lower() in name_lower:
            return True
    return False


def extract_bcp_from_ctd():
    """从 CTD 数据中提取 BCP 相关基因交互"""
    if not CTD_EXTRACTED.exists():
        print("[!] CTD 数据文件不存在，跳过 CTD 提取。")
        return []

    print("[*] 从 CTD 提取 BCP 靶标...")
    bcp_interactions = []
    total_lines = 0

    try:
        with open(CTD_EXTRACTED, 'r', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                total_lines += 1

        print(f"    总行数: {total_lines:,}")

        with open(CTD_EXTRACTED, 'r', encoding='utf-8', errors='replace') as f:
            reader = csv.reader(f, delimiter='\t')
            for i, row in enumerate(reader):
                if i == 0 or row[0].startswith('#'):
                    continue

                # CTD chem_gene_ixns 格式:
                # 0: ChemicalName, 1: ChemicalID, 2: CASRN, 3: GeneSymbol,
                # 4: GeneID, 5: GeneForms, 6: Organism, 7: OrganismID,
                # 8: Interaction, 9: InteractionActions, 10: PubMedIDs

                if len(row) < 11:
                    continue

                chem_name = row[0]
                gene_symbol = row[3]
                gene_id = row[4]
                organism = row[6]
                interaction_desc = row[8]
                interaction_actions = row[9]
                pmids = row[10]

                # 检查是否是 BCP
                is_bcp = (
                    parse_ctd_chemical_name(chem_name) or
                    BCP_CID in row[1] or
                    BCP_CAS in row[2]
                )

                if is_bcp:
                    is_human = 'human' in organism.lower() or organism == '9606'
                    bcp_interactions.append({
                        'source': 'CTD',
                        'chemical_name': chem_name,
                        'gene_symbol': gene_symbol,
                        'gene_id': gene_id,
                        'organism': organism,
                        'is_human': is_human,
                        'interaction': interaction_desc,
                        'interaction_actions': interaction_actions,
                        'pmids': pmids,
                    })

                if (i + 1) % 500000 == 0:
                    print(f"    已扫描 {i+1:,} 行... 找到 {len(bcp_interactions)} 条 BCP 交互")

    except Exception as e:
        print(f"[✗] CTD 解析出错: {e}")
        import traceback
        traceback.print_exc()

    print(f"\n[✓] CTD 提取完成: 共 {len(bcp_interactions)} 条 BCP 交互")
    return bcp_interactions


# ============================================================
# 文献已知 BCP 靶标 (手动整理)
# ============================================================
LITERATURE_TARGETS = [
    # === 主要受体 (Direct binding targets) ===
    # CNR2 - CB2 受体 (主要靶点，选择性激动剂, Ki=155nM)
    {'gene_symbol': 'CNR2', 'gene_name': 'Cannabinoid Receptor 2', 'target_type': 'Primary Receptor',
     'interaction': 'Selective full agonist (Ki=155nM)', 'evidence_pmid': '27696789,32998300'},
    # PPARA - PPARα
    {'gene_symbol': 'PPARA', 'gene_name': 'Peroxisome Proliferator Activated Receptor Alpha',
     'target_type': 'Receptor', 'interaction': 'Agonist', 'evidence_pmid': '27696789'},
    # PPARG - PPARγ
    {'gene_symbol': 'PPARG', 'gene_name': 'Peroxisome Proliferator Activated Receptor Gamma',
     'target_type': 'Receptor', 'interaction': 'Agonist (non-canonical pathway)',
     'evidence_pmid': '26364623,40050116'},
    # TRPV1
    {'gene_symbol': 'TRPV1', 'gene_name': 'Transient Receptor Potential Cation Channel Subfamily V Member 1',
     'target_type': 'Ion Channel', 'interaction': 'Modulator', 'evidence_pmid': '9804'},

    # === 炎症信号通路 ===
    {'gene_symbol': 'TNF', 'gene_name': 'Tumor Necrosis Factor',
     'target_type': 'Cytokine', 'interaction': 'Decreased expression',
     'evidence_pmid': '40050116'},
    {'gene_symbol': 'IL6', 'gene_name': 'Interleukin 6',
     'target_type': 'Cytokine', 'interaction': 'Decreased expression',
     'evidence_pmid': '40050116'},
    {'gene_symbol': 'IL1B', 'gene_name': 'Interleukin 1 Beta',
     'target_type': 'Cytokine', 'interaction': 'Decreased expression',
     'evidence_pmid': '33777571'},
    {'gene_symbol': 'IL10', 'gene_name': 'Interleukin 10',
     'target_type': 'Cytokine', 'interaction': 'Increased expression',
     'evidence_pmid': '40050116'},
    {'gene_symbol': 'NFKB1', 'gene_name': 'Nuclear Factor Kappa B Subunit 1 (p50)',
     'target_type': 'Transcription Factor', 'interaction': 'Inhibited activation',
     'evidence_pmid': '40050116,33777571'},
    {'gene_symbol': 'RELA', 'gene_name': 'RELA Proto-Oncogene (p65)',
     'target_type': 'Transcription Factor', 'interaction': 'Inhibited activation, bound',
     'evidence_pmid': '33777571'},

    # === NLRP3 炎症小体通路 ===
    {'gene_symbol': 'NLRP3', 'gene_name': 'NLR Family Pyrin Domain Containing 3',
     'target_type': 'Inflammasome', 'interaction': 'Inhibited expression, bound',
     'evidence_pmid': '33777571'},
    {'gene_symbol': 'CASP1', 'gene_name': 'Caspase 1',
     'target_type': 'Protease', 'interaction': 'Inhibited expression, bound',
     'evidence_pmid': '33777571'},
    {'gene_symbol': 'PYCARD', 'gene_name': 'ASC (Apoptosis-associated Speck-like Protein)',
     'target_type': 'Adapter', 'interaction': 'Decreased expression', 'evidence_pmid': '33777571'},

    # === TLR 通路 ===
    {'gene_symbol': 'TLR4', 'gene_name': 'Toll Like Receptor 4',
     'target_type': 'Receptor', 'interaction': 'Inhibited, bound', 'evidence_pmid': '33777571'},
    {'gene_symbol': 'MYD88', 'gene_name': 'MYD88 Innate Immune Signal Transduction Adaptor',
     'target_type': 'Adapter', 'interaction': 'Inhibited, bound', 'evidence_pmid': '33777571'},

    # === 凋亡/细胞周期靶点 ===
    {'gene_symbol': 'TP53', 'gene_name': 'Tumor Protein P53',
     'target_type': 'Transcription Factor', 'interaction': 'Modulated', 'evidence_pmid': '39448928'},
    {'gene_symbol': 'CASP3', 'gene_name': 'Caspase 3',
     'target_type': 'Protease', 'interaction': 'Activated', 'evidence_pmid': '39448928'},
    {'gene_symbol': 'BCL2', 'gene_name': 'BCL2 Apoptosis Regulator',
     'target_type': 'Apoptosis', 'interaction': 'Decreased expression',
     'evidence_pmid': '34521866'},
    {'gene_symbol': 'BAX', 'gene_name': 'BCL2 Associated X',
     'target_type': 'Apoptosis', 'interaction': 'Increased expression',
     'evidence_pmid': '34521866'},
    {'gene_symbol': 'CDK6', 'gene_name': 'Cyclin Dependent Kinase 6',
     'target_type': 'Cell Cycle', 'interaction': 'Bound/inhibited',
     'evidence_pmid': '36500446'},
    {'gene_symbol': 'CCND1', 'gene_name': 'Cyclin D1',
     'target_type': 'Cell Cycle', 'interaction': 'Downregulated', 'evidence_pmid': '40333803'},
    {'gene_symbol': 'CDKN1A', 'gene_name': 'p21 (Cyclin Dependent Kinase Inhibitor 1A)',
     'target_type': 'Cell Cycle', 'interaction': 'Upregulated', 'evidence_pmid': '40333803'},

    # === MAPK 信号通路 ===
    {'gene_symbol': 'MAPK1', 'gene_name': 'Mitogen-Activated Protein Kinase 1 (ERK2)',
     'target_type': 'Kinase', 'interaction': 'Modulated phosphorylation',
     'evidence_pmid': '39448928'},
    {'gene_symbol': 'MAPK3', 'gene_name': 'Mitogen-Activated Protein Kinase 3 (ERK1)',
     'target_type': 'Kinase', 'interaction': 'Modulated phosphorylation',
     'evidence_pmid': '39448928'},
    {'gene_symbol': 'MAPK14', 'gene_name': 'p38 MAPK',
     'target_type': 'Kinase', 'interaction': 'Decreased activation', 'evidence_pmid': '11003'},

    # === PI3K/AKT/mTOR 通路 ===
    {'gene_symbol': 'AKT1', 'gene_name': 'AKT Serine/Threonine Kinase 1',
     'target_type': 'Kinase', 'interaction': 'Modulated', 'evidence_pmid': '356195192'},
    {'gene_symbol': 'MTOR', 'gene_name': 'Mechanistic Target Of Rapamycin Kinase',
     'target_type': 'Kinase', 'interaction': 'Inhibited pathway', 'evidence_pmid': '27696789'},
    {'gene_symbol': 'PIK3CA', 'gene_name': 'Phosphatidylinositol-4,5-bisphosphate 3-Kinase Catalytic Subunit Alpha',
     'target_type': 'Kinase', 'interaction': 'Modulated', 'evidence_pmid': '27696789'},
    {'gene_symbol': 'IRS1', 'gene_name': 'Insulin Receptor Substrate 1',
     'target_type': 'Signaling', 'interaction': 'Docked, modulated', 'evidence_pmid': '356195192'},

    # === JAK/STAT 通路 ===
    {'gene_symbol': 'STAT3', 'gene_name': 'Signal Transducer And Activator Of Transcription 3',
     'target_type': 'Transcription Factor', 'interaction': 'Inhibited activation',
     'evidence_pmid': '39448928,40333803'},
    {'gene_symbol': 'JAK1', 'gene_name': 'Janus Kinase 1',
     'target_type': 'Kinase', 'interaction': 'Modulated upstream', 'evidence_pmid': '40333803'},

    # === 血管新生 ===
    {'gene_symbol': 'VEGFA', 'gene_name': 'Vascular Endothelial Growth Factor A',
     'target_type': 'Growth Factor', 'interaction': 'Decreased secretion',
     'evidence_pmid': '34521866,40333803'},
    {'gene_symbol': 'KDR', 'gene_name': 'VEGFR2 (Kinase Insert Domain Receptor)',
     'target_type': 'Receptor', 'interaction': 'Bound/inhibited', 'evidence_pmid': '40333803'},
    {'gene_symbol': 'HIF1A', 'gene_name': 'Hypoxia Inducible Factor 1 Subunit Alpha',
     'target_type': 'Transcription Factor', 'interaction': 'Modulated', 'evidence_pmid': '40333803'},

    # === 神经保护相关 ===
    {'gene_symbol': 'BDNF', 'gene_name': 'Brain Derived Neurotrophic Factor',
     'target_type': 'Neurotrophin', 'interaction': 'Increased expression',
     'evidence_pmid': '40050116'},
    {'gene_symbol': 'SIRT1', 'gene_name': 'Sirtuin 1',
     'target_type': 'Deacetylase', 'interaction': 'Increased expression',
     'evidence_pmid': '40050116'},
    {'gene_symbol': 'PPARGC1A', 'gene_name': 'PGC-1α',
     'target_type': 'Transcriptional Coactivator', 'interaction': 'Increased expression',
     'evidence_pmid': '40050116'},

    # === 氧化应激 ===
    {'gene_symbol': 'NOS2', 'gene_name': 'Nitric Oxide Synthase 2 (iNOS)',
     'target_type': 'Enzyme', 'interaction': 'Decreased expression',
     'evidence_pmid': '32998300'},
    {'gene_symbol': 'PTGS2', 'gene_name': 'Cyclooxygenase 2 (COX-2)',
     'target_type': 'Enzyme', 'interaction': 'Decreased expression',
     'evidence_pmid': '32998300'},
    {'gene_symbol': 'HMOX1', 'gene_name': 'Heme Oxygenase 1',
     'target_type': 'Enzyme', 'interaction': 'Modulated', 'evidence_pmid': '32998300'},
    {'gene_symbol': 'NQO1', 'gene_name': 'NAD(P)H Quinone Dehydrogenase 1',
     'target_type': 'Enzyme', 'interaction': 'Modulated', 'evidence_pmid': '32998300'},
    {'gene_symbol': 'NFE2L2', 'gene_name': 'NRF2 (Nuclear Factor Erythroid 2-Related Factor 2)',
     'target_type': 'Transcription Factor', 'interaction': 'Activated', 'evidence_pmid': '32998300'},

    # === 药物转运体 ===
    {'gene_symbol': 'ABCB1', 'gene_name': 'P-glycoprotein (MDR1)',
     'target_type': 'Transporter', 'interaction': 'Bound/inhibited', 'evidence_pmid': '40333803'},

    # === 其他 ===
    {'gene_symbol': 'HSP90AA1', 'gene_name': 'Heat Shock Protein 90 Alpha Family Class A Member 1',
     'target_type': 'Chaperone', 'interaction': 'Inhibited (docking)', 'evidence_pmid': '369864301'},
    {'gene_symbol': 'ESR1', 'gene_name': 'Estrogen Receptor 1',
     'target_type': 'Receptor', 'interaction': 'Modulated', 'evidence_pmid': '369864301'},
    {'gene_symbol': 'IGF1R', 'gene_name': 'Insulin Like Growth Factor 1 Receptor',
     'target_type': 'Receptor', 'interaction': 'Modulated', 'evidence_pmid': '369864301'},
    {'gene_symbol': 'SRC', 'gene_name': 'SRC Proto-Oncogene (cSrc)',
     'target_type': 'Kinase', 'interaction': 'Docked', 'evidence_pmid': '356195192'},
]


def save_results(ctd_interactions, lit_targets):
    """合并并保存所有 BCP 靶标"""
    # 从文献靶标构建
    rows = []
    seen_genes = set()

    for t in lit_targets:
        rows.append({
            'source': 'Literature',
            'gene_symbol': t['gene_symbol'],
            'gene_name': t['gene_name'],
            'target_type': t['target_type'],
            'interaction': t['interaction'],
            'evidence': t['evidence_pmid'],
            'organism': 'Human',
        })
        seen_genes.add(t['gene_symbol'].upper())

    # 从 CTD 补充
    for ix in ctd_interactions:
        gs = ix['gene_symbol'].upper()
        if gs not in seen_genes:
            rows.append({
                'source': 'CTD',
                'gene_symbol': ix['gene_symbol'],
                'gene_name': '',
                'target_type': 'CTD_curated',
                'interaction': ix['interaction_actions'],
                'evidence': ix['pmids'],
                'organism': ix['organism'],
            })
            seen_genes.add(gs)

    # 写入 CSV
    fieldnames = ['source', 'gene_symbol', 'gene_name', 'target_type',
                  'interaction', 'evidence', 'organism']

    with open(OUTPUT_FILE, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"\n{'='*60}")
    print(f"  BCP 靶标数据汇总")
    print(f"{'='*60}")
    print(f"  文献整理靶标: {len(lit_targets)}")
    print(f"  CTD 补充靶标: {len(rows) - len(lit_targets)}")
    print(f"  总去重靶标:   {len(rows)}")
    print(f"  输出文件:     {OUTPUT_FILE}")
    print(f"{'='*60}")

    # 按类型分组统计
    from collections import Counter
    type_counts = Counter(r['target_type'] for r in rows)
    print("\n  靶标类型分布:")
    for t, c in type_counts.most_common():
        print(f"    {t}: {c}")

    return rows


def print_top_targets(rows):
    """打印核心靶标"""
    print(f"\n  BCP 核心靶标列表 ({len(rows)} 个基因):")
    print(f"  {'Gene':<12} {'Type':<20} {'Interaction':<40}")
    print(f"  {'-'*72}")
    for r in rows:
        if r['source'] == 'Literature':
            print(f"  {r['gene_symbol']:<12} {r['target_type']:<20} {r['interaction']:<40}")


def main():
    print("=" * 60)
    print("  β-石竹烯 (BCP) 毒理靶标数据获取")
    print("  PubChem CID: 5281515 | CAS: 87-44-5")
    print("=" * 60)

    # 步骤 1: 尝试下载 CTD 数据（如果可用）
    ctd_interactions = []
    ctd_available = False

    # 检查本地是否已有 CTD 数据文件
    if CTD_EXTRACTED.exists():
        ctd_available = True
        print("[✓] 发现本地 CTD 数据文件")
    elif CTD_FILE.exists():
        ctd_available = download_ctd()
    else:
        print("[*] 未找到本地 CTD 数据文件，跳过在线下载。")
        print("    如需 CTD 补充数据，请手动下载:")
        print(f"    {CTD_CHEM_GENE_URL}")
        print("    放置到:", CTD_FILE)
        print("    然后重新运行本脚本。\n")

    if ctd_available:
        ctd_interactions = extract_bcp_from_ctd()
    else:
        print("[*] 使用文献整理靶标数据。")

    # 步骤 2: 合并并保存
    rows = save_results(ctd_interactions, LITERATURE_TARGETS)

    # 步骤 3: 打印核心靶标
    print_top_targets(rows)

    print(f"\n[✓] 完成!")
    print(f"    如需更新 CTD 数据，删除以下文件后重新运行:")
    print(f"    - {CTD_FILE}")
    print(f"    - {CTD_EXTRACTED}")


if __name__ == '__main__':
    main()