# -*- coding: utf-8 -*-
"""
蛋白质特征提取脚本：基因Symbol → UniProt ID → 蛋白质序列 → ESM-2特征向量 → CSV

功能说明：
1. 从 gene_list.txt 读取基因Symbol列表（每行一个）
2. 通过 MyGene.info API 将基因Symbol映射为 UniProt ID（Swiss-Prot reviewed）
3. 从 UniProt REST API 获取蛋白质序列（优先Swiss-Prot reviewed）
4. 使用 HuggingFace transformers 加载 ESM-2 模型提取特征
5. 对最后一层隐藏状态做平均池化，得到特征向量
6. 可选 PCA 降维至256维
7. 结果保存为 gene_features.csv
8. 备选方案：AAC（氨基酸组成）20维特征

依赖安装：
    pip install torch transformers biopython pandas numpy scikit-learn requests

模型 auto-download 后会缓存至 %USERPROFILE%/.cache/huggingface/hub
"""

import os
import sys
import time
import logging
import warnings
from typing import Optional

# HuggingFace镜像设置（国内网络无法直连huggingface.co时使用）
os.environ.setdefault("HF_ENDPOINT", "https://hf-mirror.com")

import numpy as np
import pandas as pd
import requests

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)

# ==================== 配置区 ====================

GENE_LIST_FILE = "gene_list.txt"           # 输入基因列表
OUTPUT_CSV = "gene_features.csv"           # 输出特征文件
ESM2_MODEL_NAME = "facebook/esm2_t6_8M_UR50D"  # ESM-2模型（320维）
# 其他可选模型（按需切换）:
# "facebook/esm2_t12_35M_UR50D"   -> 480维
# "facebook/esm2_t30_150M_UR50D"  -> 640维
# "facebook/esm2_t33_650M_UR50D"  -> 1280维
# "facebook/esm2_t36_3B_UR50D"    -> 2560维

DO_PCA = True                             # 是否进行PCA降维
PCA_N_COMPONENTS = 256                    # 降维目标维度（需<=原始维度）
USE_AAC_FALLBACK = True                   # 序列获取失败时是否使用AAC特征
MAX_SEQUENCE_LENGTH = 1022                # ESM-2最大序列长度（超出则截断）
REQUEST_TIMEOUT = 30                      # HTTP请求超时（秒）
REQUEST_DELAY = 0.3                       # API请求间隔（避免限流）


# ==================== 第一步：读取基因列表 ====================

def read_gene_list(filepath: str) -> list[str]:
    """从文件读取基因Symbol列表，每行一个，跳过空行和注释行"""
    if not os.path.exists(filepath):
        logger.error(f"基因列表文件不存在: {filepath}")
        logger.info("请创建 %s 文件，每行一个基因Symbol（如 TLR4）", filepath)
        sys.exit(1)
    genes = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                genes.append(line.upper())
    if not genes:
        logger.error("基因列表为空，请检查 %s", filepath)
        sys.exit(1)
    logger.info("读取到 %d 个基因: %s", len(genes), genes[:10])
    if len(genes) > 10:
        logger.info("  ... 共 %d 个基因", len(genes))
    return genes


# ==================== 第二步：基因Symbol → UniProt ID ====================

def map_symbol_to_uniprot(gene_symbol: str) -> Optional[str]:
    """
    通过MyGene.info API将基因Symbol映射为UniProt ID (Swiss-Prot)
    返回第一个reviewed的UniProt ID，失败则返回None
    """
    url = "https://mygene.info/v3/query"
    params = {
        "q": gene_symbol,
        "scopes": "symbol",
        "fields": "uniprot",
        "species": "human",
        "size": 5,
    }
    try:
        resp = requests.get(url, params=params, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        hits = data.get("hits", [])
        for hit in hits:
            uniprot_info = hit.get("uniprot", {})
            if isinstance(uniprot_info, dict):
                swissprot = uniprot_info.get("Swiss-Prot")
                if swissprot:
                    if isinstance(swissprot, list):
                        return swissprot[0]
                    return swissprot
            elif isinstance(uniprot_info, list):
                for entry in uniprot_info:
                    if isinstance(entry, dict) and "Swiss-Prot" in entry:
                        val = entry["Swiss-Prot"]
                        return val[0] if isinstance(val, list) else val
        # 如果没找到Swiss-Prot，尝试取任何UniProt
        for hit in hits:
            uniprot_info = hit.get("uniprot", {})
            if isinstance(uniprot_info, dict):
                for key in ("TrEMBL",):
                    val = uniprot_info.get(key)
                    if val:
                        return val[0] if isinstance(val, list) else val
        logger.warning("  [%s] MyGene.info未找到UniProt映射", gene_symbol)
        return None
    except requests.exceptions.RequestException as e:
        logger.warning("  [%s] MyGene.info请求失败: %s", gene_symbol, e)
        return None


# ==================== 第三步：获取蛋白质序列 ====================

def fetch_protein_sequence(uniprot_id: str) -> Optional[str]:
    """
    通过UniProt REST API获取蛋白质序列
    优先获取reviewed (Swiss-Prot) 序列
    返回: 氨基酸序列字符串，失败返回None
    """
    # 优先使用Swiss-Prot格式
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    try:
        resp = requests.get(url, timeout=REQUEST_TIMEOUT)
        resp.raise_for_status()
        lines = resp.text.strip().splitlines()
        # 跳过FASTA头行（以>开头）
        seq_lines = [l.strip() for l in lines if not l.startswith(">") and l.strip()]
        sequence = "".join(seq_lines)
        if sequence:
            return sequence
        logger.warning("  [%s] 序列为空", uniprot_id)
        return None
    except requests.exceptions.RequestException as e:
        logger.warning("  [%s] UniProt请求失败: %s", uniprot_id, e)
        return None


# ==================== 第四步：AAC特征（备选方案） ====================

AMINO_ACIDS = [
    "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
    "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V",
]

def compute_aac(sequence: str) -> np.ndarray:
    """
    计算氨基酸组成（Amino Acid Composition）
    返回20维向量，每个元素为对应氨基酸的频率
    """
    seq = sequence.upper()
    total = len(seq)
    if total == 0:
        return np.zeros(len(AMINO_ACIDS), dtype=np.float32)
    aa_to_idx = {aa: i for i, aa in enumerate(AMINO_ACIDS)}
    counts = np.zeros(len(AMINO_ACIDS), dtype=np.float32)
    for aa in seq:
        if aa in aa_to_idx:
            counts[aa_to_idx[aa]] += 1
    return counts / total


# ==================== 第五步：ESM-2特征提取 ====================

def load_esm2_model(model_name: str):
    """
    加载ESM-2模型和tokenizer
    返回: (model, tokenizer, hidden_size)
    """
    logger.info("正在加载ESM-2模型: %s ...", model_name)
    try:
        from transformers import AutoModel, AutoTokenizer
    except ImportError:
        logger.error(
            "缺少transformers库，请安装: pip install transformers torch"
        )
        sys.exit(1)

    try:
        tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=False)
        model = AutoModel.from_pretrained(model_name, local_files_only=False)
    except Exception as e:
        logger.error("模型加载失败: %s", e)
        logger.info(
            "请先运行以下命令下载模型（或在有网络环境下首次运行自动下载）:\n"
            "    python -c \"from transformers import AutoModel, AutoTokenizer; "
            "AutoModel.from_pretrained('%s'); "
            "AutoTokenizer.from_pretrained('%s')\"",
            model_name, model_name,
        )
        raise RuntimeError(f"ESM-2模型加载失败: {e}")

    hidden_size = model.config.hidden_size
    logger.info(
        "模型加载成功: %s | 参数: %.1fM | 隐藏层维度: %d | 层数: %d",
        model_name,
        model.config.num_hidden_layers * model.config.hidden_size
        * model.config.intermediate_size // 1e6,
        hidden_size,
        model.config.num_hidden_layers,
    )
    model.eval()
    return model, tokenizer, hidden_size


def extract_esm2_embedding(
    sequence: str,
    model,
    tokenizer,
    max_length: int = 1022,
) -> Optional[np.ndarray]:
    """
    使用ESM-2模型提取蛋白质序列特征向量
    采用最后一层隐藏状态的平均池化
    返回: (hidden_size,) 维numpy数组
    """
    import torch

    # 截断超长序列
    seq = sequence[:max_length]

    # tokenize: ESM-2在序列首尾自动添加 <cls> 和 <eos>
    inputs = tokenizer(
        seq,
        return_tensors="pt",
        add_special_tokens=True,
        truncation=True,
        max_length=max_length + 2,
    )

    with torch.no_grad():
        outputs = model(**inputs)
        # last_hidden_state: (1, seq_len, hidden_size)
        last_hidden = outputs.last_hidden_state.squeeze(0)  # (seq_len, hidden_size)

    # 平均池化（排除特殊token或包含特殊token均可，此处对所有token平均）
    embedding = last_hidden.mean(dim=0).cpu().numpy().astype(np.float32)
    return embedding


# ==================== 第六步：主流程 ====================

def main():
    logger.info("=" * 60)
    logger.info("蛋白质特征提取管道 - ESM-2")
    logger.info("=" * 60)

    # 1. 读取基因列表
    genes = read_gene_list(GENE_LIST_FILE)

    # 2. 基因Symbol → UniProt ID
    logger.info("\n[步骤1/4] 基因Symbol → UniProt ID ...")
    gene_to_uniprot = {}
    for i, gene in enumerate(genes, 1):
        logger.info("  (%d/%d) %s", i, len(genes), gene)
        uniprot_id = map_symbol_to_uniprot(gene)
        if uniprot_id:
            gene_to_uniprot[gene] = uniprot_id
            logger.info("    -> %s", uniprot_id)
        else:
            logger.warning("    -> 未找到UniProt映射，跳过")
        time.sleep(REQUEST_DELAY)

    if not gene_to_uniprot:
        logger.error("所有基因均未找到UniProt映射，程序终止")
        sys.exit(1)
    logger.info("成功映射 %d / %d 个基因", len(gene_to_uniprot), len(genes))

    # 3. UniProt ID → 蛋白质序列
    logger.info("\n[步骤2/4] 获取蛋白质序列 ...")
    gene_to_seq = {}
    failed_genes = []
    for i, (gene, uniprot_id) in enumerate(gene_to_uniprot.items(), 1):
        logger.info("  (%d/%d) %s -> %s", i, len(gene_to_uniprot), gene, uniprot_id)
        seq = fetch_protein_sequence(uniprot_id)
        if seq:
            gene_to_seq[gene] = seq
            logger.info("    序列长度: %d aa", len(seq))
        else:
            logger.warning("    序列获取失败")
            failed_genes.append(gene)
        time.sleep(REQUEST_DELAY)

    if not gene_to_seq:
        logger.error("所有基因序列获取失败，程序终止")
        sys.exit(1)
    logger.info("成功获取 %d / %d 个序列", len(gene_to_seq), len(gene_to_uniprot))

    # 4. 特征提取
    # 4a. 尝试加载ESM-2（若失败且USE_AAC_FALLBACK=True则用AAC）
    use_esm2 = True
    model = tokenizer = hidden_size = None
    try:
        model, tokenizer, hidden_size = load_esm2_model(ESM2_MODEL_NAME)
    except Exception as e:
        logger.warning("ESM-2模型加载失败: %s", e)
        if USE_AAC_FALLBACK:
            logger.info("切换至AAC（氨基酸组成）备选方案（20维特征）")
            use_esm2 = False
        else:
            logger.error(
                "ESM-2加载失败且AAC备选未启用（USE_AAC_FALLBACK=False），程序终止"
            )
            sys.exit(1)

    logger.info("\n[步骤3/4] 提取特征向量 ...")
    feature_dict = {}
    for i, (gene, seq) in enumerate(gene_to_seq.items(), 1):
        logger.info("  (%d/%d) %s ...", i, len(gene_to_seq), gene)
        try:
            if use_esm2:
                emb = extract_esm2_embedding(seq, model, tokenizer, MAX_SEQUENCE_LENGTH)
                logger.info("    特征维度: %d", len(emb))
            else:
                emb = compute_aac(seq)
                logger.info("    AAC特征维度: %d", len(emb))
            feature_dict[gene] = emb
        except Exception as e:
            logger.error("    %s 特征提取失败: %s", gene, e)
            if USE_AAC_FALLBACK and use_esm2:
                logger.info("    使用AAC备选特征")
                feature_dict[gene] = compute_aac(seq)

    if not feature_dict:
        logger.error("特征提取全部失败，程序终止")
        sys.exit(1)

    # 5. 组装DataFrame
    logger.info("\n[步骤4/4] 保存结果 ...")
    gene_list_ordered = [g for g in genes if g in feature_dict]
    all_embs = [feature_dict[g] for g in gene_list_ordered]
    dim = all_embs[0].shape[0]
    feature_array = np.array(all_embs, dtype=np.float32)
    logger.info("特征矩阵形状: %s", feature_array.shape)

    # 6. 可选PCA降维
    apply_pca = DO_PCA and use_esm2 and PCA_N_COMPONENTS < dim
    if apply_pca:
        logger.info("正在进行PCA降维: %d -> %d ...", dim, PCA_N_COMPONENTS)
        from sklearn.decomposition import PCA as PCATransformer
        pca = PCATransformer(n_components=min(PCA_N_COMPONENTS, len(gene_list_ordered) - 1))
        feature_array = pca.fit_transform(feature_array)
        logger.info("PCA降维完成，形状: %s", feature_array.shape)
        logger.info("PCA解释方差比: %.4f (前%d个主成分)",
                     pca.explained_variance_ratio_.sum(), PCA_N_COMPONENTS)

    # 7. 保存CSV
    final_dim = feature_array.shape[1]
    columns = ["gene_symbol"] + [f"feat_{i+1:04d}" for i in range(final_dim)]
    df = pd.DataFrame(feature_array, columns=columns[1:])
    df.insert(0, "gene_symbol", gene_list_ordered)
    df.to_csv(OUTPUT_CSV, index=False, encoding="utf-8")
    logger.info("结果已保存: %s", OUTPUT_CSV)
    logger.info("  - 基因数量: %d", len(gene_list_ordered))
    logger.info("  - 特征维度: %d", final_dim)
    logger.info("  - 文件大小: %.1f KB", os.path.getsize(OUTPUT_CSV) / 1024)

    # 8. 打印失败列表
    all_failed = (
        [g for g in genes if g not in gene_to_uniprot]
        + [g for g in failed_genes if g in gene_to_uniprot]
        + [g for g in genes if g in gene_to_uniprot and g not in feature_dict]
    )
    if all_failed:
        logger.warning("以下基因处理失败: %s", all_failed)

    logger.info("\n全部完成！")


# ==================== 入口 ====================

if __name__ == "__main__":
    main()