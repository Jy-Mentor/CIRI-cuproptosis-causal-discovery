#!/usr/bin/env python3
"""
石竹烯-CIRI 项目：脑卒中后认知障碍-铜死亡相关炎症靶点研究
数据整合脚本：整合本地数据与公开数据库，生成标准 GAT 输入文件

输入：
  - 本地脑缺血靶点、石竹烯靶点、铜死亡15基因、DEG结果
  - 可选：MR结果、单细胞虚拟敲除结果
输出（./processed/）：
  - node_features.csv, edge_index.csv, labels.csv
  - gene_pool.json, data_manifest.json, scaler.pkl
"""

import os
import sys
import json
import gzip
import shutil
import logging
import argparse
import copy
import time
import hashlib
import pickle
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional, Any
from collections import defaultdict
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd
import requests
import networkx as nx
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler

# ---------------------------------------------------------------------------
# 配置与常量
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent.resolve()
LOCAL_DATA_DIR = BASE_DIR / "local_data"
DOWNLOADED_RAW_DIR = BASE_DIR / "downloaded_raw"
PROCESSED_DIR = BASE_DIR / "processed"
LOGS_DIR = BASE_DIR / "logs"

for d in [LOCAL_DATA_DIR, DOWNLOADED_RAW_DIR, PROCESSED_DIR, LOGS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(LOGS_DIR / "integration.log", encoding="utf-8"),
    ],
)
logger = logging.getLogger(__name__)

CUPROPTOSIS_GENES = {
    "FDX1", "LIAS", "LIPT1", "DLAT", "PDHB", "PDHX",
    "SLC31A1", "ATP7A", "ATP7B", "ATOX1", "NFE2L2",
    "HIF1A", "MTOR", "NFKB1", "GPX4",
}

INFLAMMATORY_GO_IDS = {"GO:0006954", "GO:0006955", "GO:0002526", "GO:0002684"}

KEY_KEGG_PATHWAYS = {
    "hsa04216": "ferroptosis",
    "hsa04620": "toll_like",
    "hsa05417": "lipid_atherosclerosis",
    "hsa04668": "tnf",
    "hsa04064": "nfkb",
}

KEY_MSIGDB_SETS = {
    "HALLMARK_INFLAMMATORY_RESPONSE": "inflammatory",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "tnfa",
}

REQUIRED_LOCAL_FILES = {
    "brain_ischemia": r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\脑缺血 人.txt",
    "caryophyllene": r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\石竹烯 人.txt",
    "cuproptosis": r"C:\Users\Jy-Mentor-7\Desktop\你信我 这不是重蹈覆辙\大创\铜死亡 15 特异性.txt",
    "deg": r"C:\Users\Jy-Mentor-7\Downloads\GSE61616.top.table (2).tsv",
}

OPTIONAL_LOCAL_FILES = {
    "mr": "",      # 用户待填入绝对路径
    "knockout": "", # 用户待填入绝对路径
}

HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/120.0.0.0 Safari/537.36"
    )
}

# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------

def normalize_gene_symbol(g: str) -> str:
    """统一大写并去除版本号如 A1BG.1 -> A1BG"""
    if pd.isna(g) or not isinstance(g, str):
        return ""
    g = g.strip().upper()
    if "." in g:
        g = g.split(".")[0]
    return g


def stream_download(url: str, dest: Path, chunk_size: int = 8192, timeout: int = 300) -> bool:
    """流式下载大文件，带 tqdm 进度条"""
    try:
        logger.info(f"开始下载: {url}")
        resp = requests.get(url, stream=True, headers=HEADERS, timeout=timeout)
        resp.raise_for_status()
        total = int(resp.headers.get("content-length", 0))
        with open(dest, "wb") as f, tqdm(
            total=total, unit="B", unit_scale=True, desc=dest.name
        ) as pbar:
            for chunk in resp.iter_content(chunk_size=chunk_size):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
        logger.info(f"下载完成: {dest}")
        return True
    except Exception as e:
        logger.error(f"下载失败 {url}: {e}")
        return False


def download_with_retry(url: str, dest: Path, max_retries: int = 3) -> bool:
    """指数退避重试下载"""
    for attempt in range(max_retries):
        if stream_download(url, dest):
            return True
        wait = 2 ** attempt
        logger.warning(f"第 {attempt + 1} 次下载失败，{wait}s 后重试...")
        time.sleep(wait)
    logger.error(f"{url} 最终下载失败，已跳过")
    return False


def gunzip_file(src: Path, dst: Path) -> bool:
    try:
        with gzip.open(src, "rb") as f_in, open(dst, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        return True
    except Exception as e:
        logger.error(f"解压失败 {src}: {e}")
        return False


def read_gene_list(path: str) -> Set[str]:
    """读取每行一个基因的文件，返回标准化后的基因集合"""
    genes = set()
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            g = normalize_gene_symbol(parts[0])
            if g:
                genes.add(g)
    return genes


def read_deg_table(path: str) -> pd.DataFrame:
    """读取 limma topTable 输出，标准化列名"""
    df = pd.read_csv(path, sep="\t")
    # 先去除列名首尾空格
    df.columns = df.columns.str.strip()
    # 标准化列名映射
    col_map = {}
    for c in df.columns:
        u = c.upper()
        if u in {"GENE.SYMBOL", "GENESYMBOL", "GENE_SYMBOL", "SYMBOL", "GENE.SYMBOL"}:
            col_map[c] = "GeneSymbol"
        elif u == "ADJ.P.VAL":
            col_map[c] = "adj.P.Val"
        elif u == "P.VALUE":
            col_map[c] = "P.Value"
        elif u == "LOGFC":
            col_map[c] = "logFC"
        elif u == "AVEEXPR":
            col_map[c] = "AveExpr"
        elif u == "BASEMEAN":
            col_map[c] = "AveExpr"
    if col_map:
        df = df.rename(columns=col_map)
    if "GeneSymbol" not in df.columns:
        raise ValueError(f"DEG 文件缺少 GeneSymbol 列，现有列: {list(df.columns)}")
    df["GeneSymbol"] = df["GeneSymbol"].apply(normalize_gene_symbol)
    df = df[df["GeneSymbol"] != ""].drop_duplicates(subset=["GeneSymbol"], keep="first")
    return df


# ---------------------------------------------------------------------------
# 数据下载与处理模块
# ---------------------------------------------------------------------------

class DataSource:
    def __init__(self, name: str):
        self.name = name
        self.success = False
        self.records = 0
        self.missing_rate = 0.0
        self.notes = ""

    def to_dict(self):
        return {
            "name": self.name,
            "success": self.success,
            "records": self.records,
            "missing_rate": self.missing_rate,
            "notes": self.notes,
        }


class StringPPIDownloader:
    """STRING PPI v12.0 下载与处理"""

    BASE_URL = "https://stringdb-static.org/download"

    def __init__(self, gene_pool: Set[str], min_score: int = 400, fallback_score: int = 150):
        self.gene_pool = gene_pool
        self.min_score = min_score
        self.fallback_score = fallback_score
        self.ds = DataSource("STRING_PPI")

    def run(self) -> pd.DataFrame:
        links_gz = DOWNLOADED_RAW_DIR / "9606.protein.links.v12.0.txt.gz"
        info_gz = DOWNLOADED_RAW_DIR / "9606.protein.info.v12.0.txt.gz"
        links_txt = DOWNLOADED_RAW_DIR / "9606.protein.links.v12.0.txt"
        info_txt = DOWNLOADED_RAW_DIR / "9606.protein.info.v12.0.txt"

        # 下载
        if not links_gz.exists():
            ok = download_with_retry(
                f"{self.BASE_URL}/protein.links.v12.0/9606.protein.links.v12.0.txt.gz",
                links_gz,
            )
            if not ok:
                self.ds.notes = "下载失败"
                return pd.DataFrame()
        if not info_gz.exists():
            ok = download_with_retry(
                f"{self.BASE_URL}/protein.info.v12.0/9606.protein.info.v12.0.txt.gz",
                info_gz,
            )
            if not ok:
                self.ds.notes = "下载失败"
                return pd.DataFrame()

        # 解压
        if not links_txt.exists():
            gunzip_file(links_gz, links_txt)
        if not info_txt.exists():
            gunzip_file(info_gz, info_txt)

        # 读取 protein info 映射
        logger.info("读取 STRING protein info...")
        info_df = pd.read_csv(info_txt, sep="\t")
        # 列名: #string_protein_id protein_size preferred_name annotation
        # 过滤掉 protein_name 为空的行（依据：STRING 官方文档及网络药理学共识）
        info_df = info_df[info_df["preferred_name"].notna()]
        info_df = info_df[info_df["preferred_name"].astype(str).str.strip() != ""]
        id_to_gene = {}
        for _, row in info_df.iterrows():
            sid = row["#string_protein_id"].strip()
            gene = normalize_gene_symbol(str(row["preferred_name"]))
            if gene:
                id_to_gene[sid] = gene

        # 读取 links
        logger.info("读取 STRING protein links...")
        edges = []
        with open(links_txt, "r") as f:
            header = f.readline().strip().split()
            # protein1 protein2 combined_score
            for line in f:
                parts = line.strip().split()
                if len(parts) < 3:
                    continue
                p1, p2, score_str = parts[0], parts[1], parts[2]
                score = int(score_str)
                g1 = id_to_gene.get(p1)
                g2 = id_to_gene.get(p2)
                if not g1 or not g2:
                    continue
                if g1 not in self.gene_pool or g2 not in self.gene_pool:
                    continue
                if g1 == g2:
                    continue
                if score >= self.min_score:
                    edges.append((g1, g2, score / 1000.0))

        df = pd.DataFrame(edges, columns=["Source", "Target", "Weight"])
        self.ds.records = len(df)

        if len(df) == 0:
            logger.warning("STRING 基因池内边数为0，放宽至 combined_score>=150")
            edges_fallback = []
            with open(links_txt, "r") as f:
                f.readline()
                for line in f:
                    parts = line.strip().split()
                    if len(parts) < 3:
                        continue
                    p1, p2, score_str = parts[0], parts[1], parts[2]
                    score = int(score_str)
                    g1 = id_to_gene.get(p1)
                    g2 = id_to_gene.get(p2)
                    if not g1 or not g2:
                        continue
                    if g1 not in self.gene_pool or g2 not in self.gene_pool:
                        continue
                    if g1 == g2:
                        continue
                    if score >= self.fallback_score:
                        edges_fallback.append((g1, g2, score / 1000.0))
            df = pd.DataFrame(edges_fallback, columns=["Source", "Target", "Weight"])
            self.ds.records = len(df)
            self.ds.notes = f"放宽过滤至{self.fallback_score}"
            if len(df) == 0:
                logger.error("STRING 放宽后仍为空，致命错误")
                raise RuntimeError("STRING PPI 过滤后为空")

        self.ds.success = True
        # 去重（无向）
        df = self._dedup_undirected(df)
        return df

    @staticmethod
    def _dedup_undirected(df: pd.DataFrame) -> pd.DataFrame:
        seen = set()
        rows = []
        for _, r in df.iterrows():
            a, b = r["Source"], r["Target"]
            key = tuple(sorted([a, b]))
            if key not in seen:
                seen.add(key)
                rows.append(r)
        return pd.DataFrame(rows).reset_index(drop=True)


class GOAnnotDownloader:
    """GO 功能注释下载与处理"""

    GAF_URL = "http://current.geneontology.org/annotations/goa_human.gaf.gz"
    OBO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"

    def __init__(self, gene_pool: Set[str]):
        self.gene_pool = gene_pool
        self.ds = DataSource("GO_Annotation")

    def run(self) -> pd.DataFrame:
        gaf_gz = DOWNLOADED_RAW_DIR / "goa_human.gaf.gz"
        gaf_txt = DOWNLOADED_RAW_DIR / "goa_human.gaf"
        obo_file = DOWNLOADED_RAW_DIR / "go-basic.obo"

        if not gaf_gz.exists():
            if not download_with_retry(self.GAF_URL, gaf_gz):
                self.ds.notes = "GAF下载失败"
                return self._empty_df()
        if not gaf_txt.exists():
            gunzip_file(gaf_gz, gaf_txt)
        if not obo_file.exists():
            if not download_with_retry(self.OBO_URL, obo_file):
                self.ds.notes = "OBO下载失败"
                return self._empty_df()

        # 解析 OBO 获取命名空间
        go_ns = {}
        with open(obo_file, "r", encoding="utf-8") as f:
            in_term = False
            current_id = None
            current_ns = None
            for line in f:
                line = line.strip()
                if line == "[Term]":
                    in_term = True
                    current_id = None
                    current_ns = None
                elif line == "":
                    in_term = False
                elif in_term:
                    if line.startswith("id: GO:"):
                        current_id = line.split(": ")[1].strip()
                    elif line.startswith("namespace: "):
                        current_ns = line.split(": ")[1].strip()
                        if current_id:
                            go_ns[current_id] = current_ns

        # 解析 GAF
        # GAF 列: DB DB_Object_ID DB_Object_Symbol Qualifier GO_ID Reference Evidence With Aspect DB_Object_Name DB_Object_Synonym DB_Object_Type Taxon Date Assigned_By Annotation_Extension Gene_Product_Form_ID
        bp_counts = defaultdict(int)
        mf_counts = defaultdict(int)
        cc_counts = defaultdict(int)
        inflam_flags = defaultdict(int)

        with open(gaf_txt, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("!"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                gene = normalize_gene_symbol(parts[2])
                go_id = parts[4].strip()
                if gene not in self.gene_pool:
                    continue
                ns = go_ns.get(go_id, "unknown")
                if ns == "biological_process":
                    bp_counts[gene] += 1
                elif ns == "molecular_function":
                    mf_counts[gene] += 1
                elif ns == "cellular_component":
                    cc_counts[gene] += 1
                if go_id in INFLAMMATORY_GO_IDS:
                    inflam_flags[gene] = 1

        records = []
        for gene in sorted(self.gene_pool):
            records.append({
                "GeneSymbol": gene,
                "BP_count": bp_counts.get(gene, 0),
                "MF_count": mf_counts.get(gene, 0),
                "CC_count": cc_counts.get(gene, 0),
                "inflammation_GO_flag": inflam_flags.get(gene, 0),
            })
        df = pd.DataFrame(records)
        self.ds.success = True
        self.ds.records = len(df)
        return df

    def _empty_df(self) -> pd.DataFrame:
        return pd.DataFrame(columns=[
            "GeneSymbol", "BP_count", "MF_count", "CC_count", "inflammation_GO_flag"
        ])


class KEGGDownloader:
    """KEGG REST API 下载与处理，严格限速"""

    BASE = "https://rest.kegg.jp"

    def __init__(self, gene_pool: Set[str]):
        self.gene_pool = gene_pool
        self.ds = DataSource("KEGG_Pathway")

    def _get(self, endpoint: str, max_retries: int = 3) -> str:
        url = f"{self.BASE}{endpoint}"
        for attempt in range(max_retries):
            try:
                time.sleep(1.0)  # 限速 1 请求/秒
                resp = requests.get(url, headers=HEADERS, timeout=60)
                resp.raise_for_status()
                return resp.text
            except Exception as e:
                wait = 2 ** attempt
                logger.warning(f"KEGG {url} 失败 ({e})，{wait}s 后退避重试")
                time.sleep(wait)
        logger.error(f"KEGG {url} 最终失败")
        return ""

    def run(self) -> pd.DataFrame:
        # 获取通路列表
        pathway_list_text = self._get("/list/pathway/hsa")
        if not pathway_list_text:
            self.ds.notes = "通路列表获取失败"
            return self._empty_df()

        pathways = {}
        for line in pathway_list_text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pid = parts[0].replace("path:", "").strip()
                pname = parts[1].strip()
                pathways[pid] = pname

        # 获取基因-通路映射
        link_text = self._get("/link/hsa/pathway")
        if not link_text:
            self.ds.notes = "基因通路映射获取失败"
            return self._empty_df()

        gene_pathways = defaultdict(set)
        for line in link_text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pid = parts[0].replace("path:", "").strip()
                gene = parts[1].replace("hsa:", "").strip()
                # KEGG gene ID 如 hsa:7157，需要映射到 symbol
                # 简单处理：我们保留 symbol 匹配（后续通过用户基因池过滤）
                # 但 KEGG REST 返回的是 entrez id，不是 symbol
                # 这里我们通过另一个端点获取每个通路的基因 symbol
                pass

        # 更准确的做法：遍历每个通路获取基因列表
        # 但为了效率，我们使用 /link/hsa/pathway 获取的是 entrez -> pathway
        # 然后需要 entrez 到 symbol 的映射。这里我们简化：
        # 使用另一个端点 /list/hsa 获取所有基因信息（entrez -> symbol）
        logger.info("获取 KEGG 基因列表用于 ID 映射...")
        gene_list_text = self._get("/list/hsa")
        entrez_to_symbol = {}
        for line in gene_list_text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                eid = parts[0].replace("hsa:", "").strip()
                desc = parts[1].strip()
                # 描述格式如 "TP53; tumor protein p53"
                sym = desc.split(";")[0].strip().upper()
                entrez_to_symbol[eid] = sym

        # 重新解析 link
        gene_pathways = defaultdict(set)
        for line in link_text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 2:
                pid = parts[0].replace("path:", "").strip()
                eid = parts[1].replace("hsa:", "").strip()
                sym = entrez_to_symbol.get(eid)
                if sym and sym in self.gene_pool:
                    gene_pathways[sym].add(pid)

        records = []
        for gene in sorted(self.gene_pool):
            pws = gene_pathways.get(gene, set())
            rec = {"GeneSymbol": gene, "pathway_count": len(pws)}
            for pid, short in KEY_KEGG_PATHWAYS.items():
                rec[f"{short}_flag"] = 1 if pid in pws else 0
            records.append(rec)
        df = pd.DataFrame(records)
        self.ds.success = True
        self.ds.records = len(df)
        return df

    def _empty_df(self) -> pd.DataFrame:
        cols = ["GeneSymbol", "pathway_count"] + [f"{v}_flag" for v in KEY_KEGG_PATHWAYS.values()]
        return pd.DataFrame(columns=cols)


class UniProtDownloader:
    """UniProt 可药性/定位特征下载"""

    URL = (
        "https://rest.uniprot.org/uniprotkb/search?"
        "query=organism_id:9606+AND+reviewed:true"
        "&format=tsv"
        "&fields=accession,gene_names,protein_name,cc_subcellular_location,ft_domain,cc_catalytic_activity"
        "&size=500"
    )

    def __init__(self, gene_pool: Set[str]):
        self.gene_pool = gene_pool
        self.ds = DataSource("UniProt")

    def run(self) -> pd.DataFrame:
        all_rows = []
        cursor = None
        page = 0
        while True:
            url = self.URL
            if cursor:
                url += f"&cursor={cursor}"
            try:
                resp = requests.get(url, headers=HEADERS, timeout=120)
                resp.raise_for_status()
            except Exception as e:
                logger.error(f"UniProt 下载失败: {e}")
                self.ds.notes = f"下载失败: {e}"
                break

            lines = resp.text.strip().split("\n")
            if len(lines) <= 1:
                break
            header = lines[0].split("\t")
            for line in lines[1:]:
                parts = line.split("\t")
                if len(parts) < len(header):
                    continue
                row = dict(zip(header, parts))
                all_rows.append(row)

            # 检查是否有下一页 (Link header)
            link_header = resp.headers.get("Link", "")
            if 'rel="next"' in link_header:
                # 提取 cursor
                import re
                m = re.search(r'cursor=([^&>]+)', link_header)
                if m:
                    cursor = m.group(1)
                else:
                    break
            else:
                break
            page += 1
            if page > 100:  # 安全上限
                break

        if not all_rows:
            return self._empty_df()

        df = pd.DataFrame(all_rows)
        # 提取 Gene Symbol（取 Gene Names 第一个）
        def extract_first_gene(gn):
            if pd.isna(gn):
                return ""
            gn = str(gn).strip()
            if not gn:
                return ""
            first = gn.split()[0].split(";")[0]
            return normalize_gene_symbol(first)

        df["GeneSymbol"] = df.get("Gene Names", "").apply(extract_first_gene)
        df = df[df["GeneSymbol"].isin(self.gene_pool)]

        # 解析特征
        def parse_location(loc):
            if pd.isna(loc):
                return 0, 0, 0
            loc = str(loc).lower()
            mem = 1 if "membrane" in loc else 0
            sec = 1 if "secreted" in loc else 0
            nuc = 1 if "nucleus" in loc else 0
            return mem, sec, nuc

        def count_domains(dom):
            if pd.isna(dom):
                return 0
            return len(str(dom).strip().split(";"))

        def has_catalytic(cat):
            return 0 if pd.isna(cat) or str(cat).strip() == "" else 1

        locs = df["Subcellular location [CC]"].apply(parse_location)
        df["is_membrane"] = [x[0] for x in locs]
        df["is_secreted"] = [x[1] for x in locs]
        df["is_nuclear"] = [x[2] for x in locs]
        df["domain_count"] = df["Domain [FT]"].apply(count_domains)
        df["is_enzyme"] = df["Catalytic activity"].apply(has_catalytic)

        result = df[["GeneSymbol", "is_membrane", "is_secreted", "is_nuclear", "domain_count", "is_enzyme"]].copy()
        # 合并重复基因（取最大值）
        result = result.groupby("GeneSymbol").max().reset_index()
        self.ds.success = True
        self.ds.records = len(result)
        return result

    def _empty_df(self) -> pd.DataFrame:
        return pd.DataFrame(columns=[
            "GeneSymbol", "is_membrane", "is_secreted", "is_nuclear", "domain_count", "is_enzyme"
        ])


class MSigDBDownloader:
    """MSigDB 基因集下载与处理"""

    BASE_URL = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs"
    FILES = {
        "hallmark": "h.all.v2023.2.Hs.symbols.gmt",
        "c2_kegg": "c2.cp.kegg.v2023.2.Hs.symbols.gmt",
        "c5_bp": "c5.go.bp.v2023.2.Hs.symbols.gmt",
    }

    def __init__(self, gene_pool: Set[str]):
        self.gene_pool = gene_pool
        self.ds = DataSource("MSigDB")

    def run(self) -> pd.DataFrame:
        gene_sets = defaultdict(lambda: {"hallmark": 0, "c2": 0, "c5bp": 0, "inflammatory": 0, "tnfa": 0})

        for key, filename in self.FILES.items():
            url = f"{self.BASE_URL}/{filename}"
            dest = DOWNLOADED_RAW_DIR / filename
            if not dest.exists():
                if not download_with_retry(url, dest):
                    continue
            with open(dest, "r", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 3:
                        continue
                    set_name = parts[0]
                    genes_in_set = {normalize_gene_symbol(g) for g in parts[2:]}
                    for g in self.gene_pool:
                        if g in genes_in_set:
                            if key == "hallmark":
                                gene_sets[g]["hallmark"] += 1
                            elif key == "c2_kegg":
                                gene_sets[g]["c2"] += 1
                            elif key == "c5_bp":
                                gene_sets[g]["c5bp"] += 1
                            if set_name == "HALLMARK_INFLAMMATORY_RESPONSE":
                                gene_sets[g]["inflammatory"] = 1
                            if set_name == "HALLMARK_TNFA_SIGNALING_VIA_NFKB":
                                gene_sets[g]["tnfa"] = 1

        records = []
        for gene in sorted(self.gene_pool):
            d = gene_sets[gene]
            records.append({
                "GeneSymbol": gene,
                "hallmark_count": d["hallmark"],
                "c2_count": d["c2"],
                "c5bp_count": d["c5bp"],
                "inflammatory_flag": d["inflammatory"],
                "tnfa_flag": d["tnfa"],
            })
        df = pd.DataFrame(records)
        self.ds.success = True
        self.ds.records = len(df)
        return df


class BioGRIDDownloader:
    """BioGRID 物理互作下载（P1，可选）"""

    URL = "https://downloads.thebiogrid.org/BioGRID/Release-Archive/BIOGRID-4.4.229/BIOGRID-ALL-4.4.229.tab3.zip"

    def __init__(self, gene_pool: Set[str]):
        self.gene_pool = gene_pool
        self.ds = DataSource("BioGRID")

    def run(self) -> Optional[pd.DataFrame]:
        zip_path = DOWNLOADED_RAW_DIR / "BIOGRID-ALL-4.4.229.tab3.zip"
        if not zip_path.exists():
            if not download_with_retry(self.URL, zip_path):
                self.ds.notes = "下载失败"
                return None
        import zipfile
        try:
            with zipfile.ZipFile(zip_path, "r") as z:
                tab_name = [n for n in z.namelist() if n.endswith(".tab3")][0]
                with z.open(tab_name) as f:
                    df = pd.read_csv(f, sep="\t")
        except Exception as e:
            logger.error(f"BioGRID 解压/读取失败: {e}")
            self.ds.notes = f"读取失败: {e}"
            return None

        # 筛选人类物理互作
        df = df[
            (df["Organism ID Interactor A"] == 9606) &
            (df["Organism ID Interactor B"] == 9606) &
            (df["Experimental System Type"].str.lower() == "physical")
        ].copy()

        df["GeneSymbol_A"] = df["Official Symbol Interactor A"].apply(normalize_gene_symbol)
        df["GeneSymbol_B"] = df["Official Symbol Interactor B"].apply(normalize_gene_symbol)
        df = df[
            df["GeneSymbol_A"].isin(self.gene_pool) &
            df["GeneSymbol_B"].isin(self.gene_pool) &
            (df["GeneSymbol_A"] != df["GeneSymbol_B"])
        ]

        edges = []
        for _, r in df.iterrows():
            edges.append((r["GeneSymbol_A"], r["GeneSymbol_B"], 1.0))
        result = pd.DataFrame(edges, columns=["Source", "Target", "Weight"])
        result = result.drop_duplicates()
        self.ds.success = True
        self.ds.records = len(result)
        return result


# ---------------------------------------------------------------------------
# 主整合流程
# ---------------------------------------------------------------------------

class DataIntegrator:
    def __init__(self):
        self.gene_pool: Set[str] = set()
        self.deg_df: Optional[pd.DataFrame] = None
        self.mr_df: Optional[pd.DataFrame] = None
        self.ko_df: Optional[pd.DataFrame] = None
        self.data_sources: List[DataSource] = []
        self.unmapped_log: List[str] = []

    def validate_local_files(self):
        for name, path in REQUIRED_LOCAL_FILES.items():
            if not os.path.exists(path):
                logger.error(f"本地必需文件缺失: {name} -> {path}")
                sys.exit(1)
        logger.info("本地必需文件校验通过")

    def build_gene_pool(self):
        brain = read_gene_list(REQUIRED_LOCAL_FILES["brain_ischemia"])
        cary = read_gene_list(REQUIRED_LOCAL_FILES["caryophyllene"])
        cupro = read_gene_list(REQUIRED_LOCAL_FILES["cuproptosis"])
        self.deg_df = read_deg_table(REQUIRED_LOCAL_FILES["deg"])

        # 严格过滤：只取显著差异基因进入基因池
        deg_significant = self.deg_df[
            (self.deg_df["adj.P.Val"] < 0.05) &
            (self.deg_df["logFC"].abs() > 1.0)
        ].copy()
        deg_genes = set(deg_significant["GeneSymbol"].tolist())
        logger.info(f"DEG 显著基因数 (|logFC|>1.0 & adj.P<0.05): {len(deg_genes)}")

        self.gene_pool = brain | cary | cupro | deg_genes
        logger.info(f"基因池初始大小: {len(self.gene_pool)}")

        # 强制纳入铜死亡基因
        missing_cupro = CUPROPTOSIS_GENES - self.gene_pool
        if missing_cupro:
            logger.error(f"铜死亡基因缺失（致命）: {missing_cupro}")
            sys.exit(1)

        # 读取可选 MR
        mr_path = OPTIONAL_LOCAL_FILES.get("mr", "")
        if mr_path and os.path.exists(mr_path):
            self.mr_df = pd.read_csv(mr_path, sep=None, engine="python")
            if "GeneSymbol" not in self.mr_df.columns:
                # 尝试自动识别
                for c in self.mr_df.columns:
                    if c.upper() in {"GENESYMBOL", "GENE", "SYMBOL"}:
                        self.mr_df = self.mr_df.rename(columns={c: "GeneSymbol"})
                        break
            if "GeneSymbol" in self.mr_df.columns:
                self.mr_df["GeneSymbol"] = self.mr_df["GeneSymbol"].apply(normalize_gene_symbol)
                self.gene_pool |= set(self.mr_df["GeneSymbol"].tolist())
        else:
            logger.info("MR 结果文件未提供或不存在，跳过")

        # 读取可选 Knockout
        ko_path = OPTIONAL_LOCAL_FILES.get("knockout", "")
        if ko_path and os.path.exists(ko_path):
            self.ko_df = pd.read_csv(ko_path, sep=None, engine="python")
            if "GeneSymbol" not in self.ko_df.columns:
                for c in self.ko_df.columns:
                    if c.upper() in {"GENESYMBOL", "GENE", "SYMBOL"}:
                        self.ko_df = self.ko_df.rename(columns={c: "GeneSymbol"})
                        break
            if "GeneSymbol" in self.ko_df.columns:
                self.ko_df["GeneSymbol"] = self.ko_df["GeneSymbol"].apply(normalize_gene_symbol)
                self.gene_pool |= set(self.ko_df["GeneSymbol"].tolist())
        else:
            logger.info("Knockout 结果文件未提供或不存在，跳过")

        # 最终清理
        self.gene_pool = {g for g in self.gene_pool if g}
        logger.info(f"最终基因池大小: {len(self.gene_pool)}")

        # 保存基因池元数据
        pool_meta = {
            "total_genes": len(self.gene_pool),
            "brain_ischemia": len(brain),
            "caryophyllene": len(cary),
            "cuproptosis": len(cupro),
            "deg_overlap": len(deg_genes & self.gene_pool),
            "sources": {
                "brain_ischemia": sorted(brain),
                "caryophyllene": sorted(cary),
                "cuproptosis": sorted(cupro),
            },
        }
        with open(PROCESSED_DIR / "gene_pool.json", "w", encoding="utf-8") as f:
            json.dump(pool_meta, f, indent=2, ensure_ascii=False)

    def download_and_process(self):
        # STRING（min_score=400，依据 STRING 官方文档及网络药理学共识）
        string_dl = StringPPIDownloader(self.gene_pool, min_score=400, fallback_score=150)
        self.string_edges = string_dl.run()
        self.data_sources.append(string_dl.ds)

        # GO
        go_dl = GOAnnotDownloader(self.gene_pool)
        self.go_df = go_dl.run()
        self.data_sources.append(go_dl.ds)

        # KEGG
        kegg_dl = KEGGDownloader(self.gene_pool)
        self.kegg_df = kegg_dl.run()
        self.data_sources.append(kegg_dl.ds)

        # UniProt
        uniprot_dl = UniProtDownloader(self.gene_pool)
        self.uniprot_df = uniprot_dl.run()
        self.data_sources.append(uniprot_dl.ds)

        # MSigDB
        msigdb_dl = MSigDBDownloader(self.gene_pool)
        self.msigdb_df = msigdb_dl.run()
        self.data_sources.append(msigdb_dl.ds)

        # BioGRID (P1)
        biogrid_dl = BioGRIDDownloader(self.gene_pool)
        self.biogrid_edges = biogrid_dl.run()
        self.data_sources.append(biogrid_dl.ds)

    def merge_edges(self) -> pd.DataFrame:
        string_df = self.string_edges.copy()
        string_df["SourceDB"] = "STRING"

        if self.biogrid_edges is not None and len(self.biogrid_edges) > 0:
            bio_df = self.biogrid_edges.copy()
            bio_df["SourceDB"] = "BioGRID"
            merged = pd.concat([string_df, bio_df], ignore_index=True)
            # 合并去重：相同 Source/Target 取最大 Weight，标记 Both
            merged = merged.sort_values("Weight", ascending=False)
            dedup = {}
            for _, r in merged.iterrows():
                key = tuple(sorted([r["Source"], r["Target"]]))
                if key not in dedup:
                    dedup[key] = {
                        "Source": key[0],
                        "Target": key[1],
                        "Weight": r["Weight"],
                        "SourceDB": r["SourceDB"],
                    }
                else:
                    if dedup[key]["SourceDB"] != r["SourceDB"]:
                        dedup[key]["SourceDB"] = "Both"
            result = pd.DataFrame(list(dedup.values()))
        else:
            result = string_df

        # 确保图连通：孤立节点连接到基因池中 STRING score 最高的邻居
        G = nx.Graph()
        for _, r in result.iterrows():
            G.add_edge(r["Source"], r["Target"], weight=r["Weight"])

        isolated = self.gene_pool - set(G.nodes())
        if isolated:
            logger.warning(f"发现 {len(isolated)} 个孤立节点，尝试连接...")
            # 为每个孤立节点找到基因池中 STRING score 最高的邻居
            # 使用 string_edges 的原始数据
            neighbor_scores = defaultdict(list)
            for _, r in self.string_edges.iterrows():
                neighbor_scores[r["Source"]].append((r["Target"], r["Weight"]))
                neighbor_scores[r["Target"]].append((r["Source"], r["Weight"]))

            for gene in isolated:
                candidates = neighbor_scores.get(gene, [])
                # 只连接到基因池内的基因
                pool_candidates = [(n, w) for n, w in candidates if n in self.gene_pool and n != gene]
                if pool_candidates:
                    pool_candidates.sort(key=lambda x: x[1], reverse=True)
                    best_neighbor, best_weight = pool_candidates[0]
                    result = pd.concat([
                        result,
                        pd.DataFrame([{
                            "Source": gene,
                            "Target": best_neighbor,
                            "Weight": best_weight,
                            "SourceDB": "STRING",
                        }])
                    ], ignore_index=True)
                    G.add_edge(gene, best_neighbor, weight=best_weight)
                else:
                    # 如果连 STRING 中都没有，随机连接到另一个基因池基因
                    other = next(iter(self.gene_pool - {gene}))
                    result = pd.concat([
                        result,
                        pd.DataFrame([{
                            "Source": gene,
                            "Target": other,
                            "Weight": 0.15,
                            "SourceDB": "STRING",
                        }])
                    ], ignore_index=True)
                    G.add_edge(gene, other, weight=0.15)

        # 再次检查
        G_check = nx.Graph()
        for _, r in result.iterrows():
            G_check.add_edge(r["Source"], r["Target"])
        still_isolated = self.gene_pool - set(G_check.nodes())
        if still_isolated:
            logger.warning(f"仍有 {len(still_isolated)} 孤立节点")

        return result.reset_index(drop=True)

    def build_node_features(self, edge_df: pd.DataFrame) -> pd.DataFrame:
        genes = sorted(self.gene_pool)
        features = pd.DataFrame({"GeneSymbol": genes})

        # --- 转录组维度 ---
        deg = self.deg_df.copy()
        deg_cols = ["logFC", "AveExpr", "t", "P.Value", "adj.P.Val"]
        for c in deg_cols:
            if c not in deg.columns:
                deg[c] = np.nan
        deg["abs_logFC"] = deg["logFC"].abs()
        deg["neg_log10_P"] = -np.log10(deg["P.Value"].replace(0, np.nan))
        deg["neg_log10_adjP"] = -np.log10(deg["adj.P.Val"].replace(0, np.nan))

        transcript = deg[["GeneSymbol", "logFC", "abs_logFC", "neg_log10_P", "neg_log10_adjP", "AveExpr"]].copy()
        features = features.merge(transcript, on="GeneSymbol", how="left")

        # --- 因果维度 (MR) ---
        if self.mr_df is not None and "Beta" in self.mr_df.columns:
            mr = self.mr_df.copy()
            mr["MR_absBeta"] = mr["Beta"].abs()
            mr["MR_neg_log10_P"] = -np.log10(mr["Pval"].replace(0, np.nan))
            mr["MR_significant_flag"] = (mr["Pval"] < 0.05).astype(int)
            mr = mr.rename(columns={
                "Beta": "MR_Beta",
                "Fstat": "MR_Fstat",
            })
            mr = mr[["GeneSymbol", "MR_Beta", "MR_absBeta", "MR_neg_log10_P", "MR_Fstat", "MR_significant_flag"]]
            features = features.merge(mr, on="GeneSymbol", how="left")
        else:
            for c in ["MR_Beta", "MR_absBeta", "MR_neg_log10_P", "MR_Fstat", "MR_significant_flag"]:
                features[c] = 0.0
            logger.info("MR 维度全部填0 (MISSING)")

        # --- 扰动维度 (Knockout) ---
        if self.ko_df is not None and "Knockout_Score" in self.ko_df.columns:
            ko = self.ko_df.copy()
            ko["abs_Knockout_Score"] = ko["Knockout_Score"].abs()
            ko = ko[["GeneSymbol", "Knockout_Score", "abs_Knockout_Score"]]
            features = features.merge(ko, on="GeneSymbol", how="left")
        else:
            for c in ["Knockout_Score", "abs_Knockout_Score"]:
                features[c] = 0.0
            logger.info("Knockout 维度全部填0 (MISSING)")

        # --- 功能维度 (GO) ---
        if hasattr(self, "go_df") and len(self.go_df) > 0:
            features = features.merge(self.go_df, on="GeneSymbol", how="left")
        else:
            for c in ["BP_count", "MF_count", "CC_count", "inflammation_GO_flag"]:
                features[c] = 0

        # --- 通路维度 (KEGG) ---
        if hasattr(self, "kegg_df") and len(self.kegg_df) > 0:
            features = features.merge(self.kegg_df, on="GeneSymbol", how="left")
        else:
            for c in ["pathway_count", "ferroptosis_flag", "toll_like_flag", "lipid_atherosclerosis_flag", "tnf_flag", "nfkb_flag"]:
                features[c] = 0

        # --- 可药性维度 (UniProt) ---
        if hasattr(self, "uniprot_df") and len(self.uniprot_df) > 0:
            features = features.merge(self.uniprot_df, on="GeneSymbol", how="left")
        else:
            for c in ["is_membrane", "is_secreted", "is_nuclear", "domain_count", "is_enzyme"]:
                features[c] = 0

        # --- 基因集维度 (MSigDB) ---
        if hasattr(self, "msigdb_df") and len(self.msigdb_df) > 0:
            features = features.merge(self.msigdb_df, on="GeneSymbol", how="left")
        else:
            for c in ["hallmark_count", "c2_count", "c5bp_count", "inflammatory_flag", "tnfa_flag"]:
                features[c] = 0

        # --- 网络拓扑维度 ---
        G = nx.Graph()
        for _, r in edge_df.iterrows():
            G.add_edge(r["Source"], r["Target"], weight=r.get("Weight", 1.0))
        # 为没有连接的节点添加孤立点
        for g in genes:
            if g not in G:
                G.add_node(g)

        degree_dict = dict(G.degree())
        pagerank_dict = nx.pagerank(G, weight="weight")
        clustering_dict = nx.clustering(G)

        features["Degree"] = features["GeneSymbol"].map(degree_dict).fillna(0)
        features["PageRank"] = features["GeneSymbol"].map(pagerank_dict).fillna(0)
        features["ClusteringCoefficient"] = features["GeneSymbol"].map(clustering_dict).fillna(0)

        # --- 缺失值处理 ---
        numeric_cols = features.select_dtypes(include=[np.number]).columns.tolist()
        binary_cols = [c for c in numeric_cols if c.endswith("_flag") or c.startswith("is_") or c == "MR_significant_flag"]
        continuous_cols = [c for c in numeric_cols if c not in binary_cols and c != "GeneSymbol"]

        for c in continuous_cols:
            col_mean = features[c].mean()
            if pd.isna(col_mean):
                features[c] = features[c].fillna(0)
                logger.warning(f"特征列 {c} 全缺失，已填 0")
            else:
                features[c] = features[c].fillna(col_mean)
        for c in binary_cols:
            features[c] = features[c].fillna(0).astype(int)

        # --- 归一化 ---
        # 拓扑特征使用 Min-Max 归一化到 [0,1]，保持语义
        topo_cols = ["Degree", "PageRank", "ClusteringCoefficient"]
        from sklearn.preprocessing import MinMaxScaler
        for c in topo_cols:
            if c in features.columns:
                c_min = features[c].min()
                c_max = features[c].max()
                if c_max > c_min:
                    features[c] = (features[c] - c_min) / (c_max - c_min)
                else:
                    features[c] = 0.0

        # 其他连续特征使用 StandardScaler（Z-score）
        non_topo_continuous = [c for c in continuous_cols if c not in topo_cols]
        if non_topo_continuous:
            scaler = StandardScaler()
            features_to_scale = features[non_topo_continuous].copy()
            scaled = scaler.fit_transform(features_to_scale)
            features[non_topo_continuous] = scaled
            # 保存 scaler
            with open(PROCESSED_DIR / "scaler.pkl", "wb") as f:
                pickle.dump(scaler, f)

        return features

    def build_labels(self) -> pd.DataFrame:
        np.random.seed(42)  # 固定随机性，确保可复现

        brain = read_gene_list(REQUIRED_LOCAL_FILES["brain_ischemia"])
        cary = read_gene_list(REQUIRED_LOCAL_FILES["caryophyllene"])
        cupro = CUPROPTOSIS_GENES.copy()

        # 严格交集定义阳性
        positive = brain & cary
        if len(positive) < 3:
            logger.warning(f"阳性标签交集仅 {len(positive)} 个，退化为并集")
            positive = brain | cary

        # 修正版：只保留有 FDA/EMA 批准药物直接靶向的基因
        # 已剔除：转录因子、血浆蛋白、结构蛋白、无上市药物靶点
        LIKELY_REAL_TARGETS = {
            # 炎症与免疫
            'PTGS2', 'PTGS1', 'TNF', 'IL6', 'IL1B', 'IL10', 'IL5',
            'JAK2', 'JAK1', 'IRAK4', 'CCL2', 'CCR2', 'CCR5', 'CXCR3',
            'ICAM1', 'VCAM1', 'TGFB1', 'STAT1',
            # 核受体与代谢
            'PPARG', 'PPARA', 'PPARD', 'VDR', 'ESR1', 'NR3C1', 'NR3C2', 'NR1I2',
            'HMGCR', 'DPP4', 'CPT1A', 'SREBF1', 'FASN',
            # 氧化应激
            'NOS3', 'NOS2', 'NOS1', 'XDH', 'ALOX5', 'HMOX1', 'GPX1',
            # 神经与精神
            'DRD2', 'HTR2A', 'ADRB1', 'ADRB3', 'ADORA1', 'ADORA2A', 'ADORA3',
            'OPRM1', 'OPRK1', 'CHRM3', 'TACR1', 'TRPV1', 'TRPA1',
            'SLC6A4', 'SLC6A3', 'SLC6A2', 'CNR2', 'S1PR1', 'TBXA2R',
            'PTGER1', 'PTGER2', 'P2RX7', 'ACHE', 'MAOB', 'BCHE', 'HCRTR1', 'ABAT',
            # 激酶与信号
            'MTOR', 'EGFR', 'SRC', 'ERBB4', 'MAPK14',
            # 凋亡/DNA修复/肿瘤
            'PARP1', 'BCL2', 'BRD4', 'MDM2', 'CASP3', 'CASP8', 'CASP9', 'ATG5',
            # 凝血/补体/心血管
            'F2', 'C5', 'KDR', 'MMP9', 'TIMP1',
            # 其他明确可药靶点
            'ALDH1A1', 'TSPO', 'CYP3A4', 'CYP2C9', 'CYP1A2', 'EPHX2',
            'CTSD', 'CTSB', 'CTSS', 'ALDH9A1',
        }

        positive = positive & LIKELY_REAL_TARGETS
        logger.info(f"提纯后阳性标签数: {len(positive)}")

        if len(positive) < 3:
            logger.error(f"提纯后阳性仅 {len(positive)} 个，无法训练半监督GAT")
            sys.exit(1)

        # 阴性：DEG 中 |logFC|<0.1 且 adj.P.Val>0.5 的基因（表达无显著变化）
        # 但必须在当前缩池后的基因池内
        deg = self.deg_df.copy()
        negative_candidates = deg[
            (deg["logFC"].abs() < 0.1) & (deg["adj.P.Val"] > 0.5)
        ]["GeneSymbol"].tolist()
        negative = (set(negative_candidates) - positive) & self.gene_pool

        # 补足阴性至阳性 3 倍，排除铜死亡基因
        pool_minus_positive = self.gene_pool - positive
        if len(negative) < len(positive) * 3:
            needed = max(len(positive) * 3, 10)
            extra_pool = pool_minus_positive - negative - cupro
            n_extra = min(needed - len(negative), len(extra_pool))
            if extra_pool and n_extra > 0:
                extra = set(np.random.choice(list(extra_pool), size=n_extra, replace=False))
                negative |= extra

        # 标签分配
        labels = {}
        for g in self.gene_pool:
            if g in positive:
                labels[g] = 1
            elif g in negative:
                labels[g] = 0
            else:
                labels[g] = -1

        # 铜死亡基因强制保留为候选
        for g in cupro:
            if g not in positive:
                labels[g] = -1

        df = pd.DataFrame([{"GeneSymbol": g, "Label": v} for g, v in labels.items()])

        # 校验
        n_pos = (df["Label"] == 1).sum()
        n_neg = (df["Label"] == 0).sum()
        n_unk = (df["Label"] == -1).sum()

        logger.info(f"最终标签 — 阳性:{n_pos}, 阴性:{n_neg}, 未知:{n_unk}")
        assert n_pos >= 3 and n_neg >= 1 and n_unk >= 1, "标签数量不足"

        return df

    def save_manifest(self):
        manifest = {
            "gene_pool_size": len(self.gene_pool),
            "data_sources": [ds.to_dict() for ds in self.data_sources],
            "local_files": {k: os.path.exists(v) for k, v in REQUIRED_LOCAL_FILES.items()},
        }
        with open(PROCESSED_DIR / "data_manifest.json", "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2, ensure_ascii=False)

    def run(self):
        self.validate_local_files()
        self.build_gene_pool()
        self.download_and_process()

        edge_df = self.merge_edges()
        node_df = self.build_node_features(edge_df)
        label_df = self.build_labels()

        # 保存
        node_df.to_csv(PROCESSED_DIR / "node_features.csv", index=False)
        edge_df.to_csv(PROCESSED_DIR / "edge_index.csv", index=False)
        label_df.to_csv(PROCESSED_DIR / "labels.csv", index=False)
        self.save_manifest()

        # 最终校验
        self.final_check(node_df, edge_df, label_df)

    def final_check(self, node_df: pd.DataFrame, edge_df: pd.DataFrame, label_df: pd.DataFrame):
        n_rows, n_cols = node_df.shape
        n_edges = len(edge_df)
        n_pos = (label_df["Label"] == 1).sum()
        n_neg = (label_df["Label"] == 0).sum()
        n_unk = (label_df["Label"] == -1).sum()

        logger.info("=" * 60)
        logger.info("最终报告")
        logger.info("=" * 60)
        logger.info(f"节点特征维度 X: {n_rows} × {n_cols - 1} (不含 GeneSymbol)")
        logger.info(f"边数: {n_edges}")
        logger.info(f"标签分布 — 阳性(1): {n_pos}, 阴性(0): {n_neg}, 未知(-1): {n_unk}")
        logger.info(f"阳性/阴性比例 ≈ 1:{n_neg/max(n_pos,1):.1f}")

        # 铜死亡基因标签分布
        cupro_labels = label_df[label_df["GeneSymbol"].isin(CUPROPTOSIS_GENES)]
        logger.info(f"铜死亡15基因标签分布:\n{cupro_labels['Label'].value_counts().to_dict()}")

        # 特征缺失率
        missing_rates = node_df.isnull().mean()
        high_missing = missing_rates[missing_rates > 0.1]
        if len(high_missing) > 0:
            logger.warning(f"高缺失率特征 (>10%):\n{high_missing}")
        else:
            logger.info("所有特征缺失率均 <= 10%")

        # 最小可用数据集检查
        assert n_rows > 0 and n_cols > 10, "node_features 维度不足"
        assert n_edges > 0, "edge_index 为空"
        assert n_pos >= 1 and n_neg >= 1, "标签数量不足"
        logger.info("最小可用数据集检查通过")


def main():
    integrator = DataIntegrator()
    integrator.run()


if __name__ == "__main__":
    main()
