# -*- coding: utf-8 -*-
"""
GAT Embedding Extractor v3 — Full Model (All Edges)
  Single GAT model trained on all edges (DT + TD + PPI) simultaneously.
  Gene embeddings (64-dim) → PCA (12-dim) used by BOTH DT and DG classifiers.
  Leak protection is handled at RF level (removing distance features).
"""
import os, sys, warnings
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import HeteroData
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore")

GAT_DIR = r"C:\Users\Jy-Mentor-7\Desktop\GAT"
RF_DIR = r"C:\Users\Jy-Mentor-7\Desktop\随机森林"
FEATURE_TABLE = os.path.join(RF_DIR, "gene_features_table_with_gat.csv")

HIDDEN_DIM = 128
OUT_DIM = 64
GAT_HEADS = 4
DROPOUT = 0.3
LR = 1e-3
WEIGHT_DECAY = 1e-5
EPOCHS = 300
PATIENCE = 40
RANDOM_SEED = 42
N_PCA = 12

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"[Init] Device: {device}")

import random
def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

set_seed(RANDOM_SEED)

# ── Data Loading ──

def load_drug_targets(path):
    if not os.path.exists(path): raise FileNotFoundError(f"Not found: {path}")
    for enc in ("utf-8","gbk","latin1"):
        try:
            with open(path,"r",encoding=enc) as fh: lines=[l.strip() for l in fh if l.strip()]
            break
        except UnicodeDecodeError: continue
    return list(dict.fromkeys([g.upper() for g in lines]))

def load_disease_genes(path):
    if not os.path.exists(path): raise FileNotFoundError(f"Not found: {path}")
    for enc in ("utf-8","gbk","latin1"):
        try:
            with open(path,"r",encoding=enc) as fh: lines=[l.strip() for l in fh if l.strip()]
            break
        except UnicodeDecodeError: continue
    if "," in lines[0] or "\t" in lines[0]:
        genes=[]
        for ln in lines:
            parts=ln.replace("\t",",").split(",")
            g=parts[0].strip().upper()
            if g and all(c.isalnum() or c=="-" for c in g): genes.append(g)
        return list(dict.fromkeys(genes))
    return list(dict.fromkeys([g.upper() for g in lines if g]))

def load_ppi(path):
    ppi_edges=[]
    for enc in ("utf-8","gbk","latin1"):
        try:
            with open(path,"r",encoding=enc) as fh:
                for line in fh:
                    line=line.strip()
                    if not line or "protein" in line.lower(): continue
                    parts=line.replace("\t",",").split(",")
                    if len(parts)>=3:
                        a,b=parts[0].strip().upper(),parts[1].strip().upper()
                        if a and b: ppi_edges.append((a,b))
            break
        except UnicodeDecodeError: continue
    return ppi_edges

def load_gene_features(path):
    ft=pd.read_csv(path)
    if "gene_symbol" in ft.columns:
        ft=ft.set_index("gene_symbol")
    arr=ft.values.astype(np.float32)
    names=ft.index.tolist()
    arr=np.nan_to_num(arr,0.0)
    return arr,names

def load_drug_fingerprint(path):
    df=pd.read_csv(path)
    fp=df.drop(columns=df.columns[:1]).values.astype(np.float32).flatten()
    return fp

def load_all_genes(path):
    with open(path) as f: return [l.strip().upper() for l in f if l.strip()]

# ── Graph Build ──

def build_graph(dt_list,dg_list,ppi_edges,gene_feat,gene_names,drug_fp,all_genes):
    data=HeteroData()

    all_genes_set=set(g.upper() for g in all_genes)
    gene_to_idx={g:i for i,g in enumerate(all_genes)}
    n_genes=len(all_genes)

    gene_feat_dict={g.upper():gene_feat[i] for i,g in enumerate(gene_names) if i<len(gene_feat)}
    gx=np.zeros((n_genes,gene_feat.shape[1]),dtype=np.float32)
    for g in all_genes:
        if g in gene_feat_dict: gx[gene_to_idx[g]]=gene_feat_dict[g]
    data["gene"].x=torch.from_numpy(gx).float()
    data["gene"].num_nodes=n_genes

    data["drug"].x=torch.from_numpy(drug_fp.reshape(1,-1)).float()
    data["drug"].num_nodes=1

    data["disease"].x=torch.ones((1,32),dtype=torch.float32)
    data["disease"].num_nodes=1

    dt_src,dt_dst=[],[]
    for g in dt_list:
        if g in gene_to_idx:
            dt_src.append(0); dt_dst.append(gene_to_idx[g])
    data["drug","targets","gene"].edge_index=torch.tensor([dt_src,dt_dst],dtype=torch.long)

    td_src,td_dst=[],[]
    for g in dg_list:
        if g in gene_to_idx:
            td_src.append(gene_to_idx[g]); td_dst.append(0)
    data["gene","associated_with","disease"].edge_index=torch.tensor([td_src,td_dst],dtype=torch.long)

    ppi_src,ppi_dst=[],[]
    for a,b in ppi_edges:
        if a in gene_to_idx and b in gene_to_idx:
            ppi_src.append(gene_to_idx[a]); ppi_dst.append(gene_to_idx[b])
    data["gene","interacts","gene"].edge_index=torch.tensor([ppi_src,ppi_dst],dtype=torch.long)

    print(f"[Build] DT:{len(dt_src)} TD:{len(td_src)} PPI:{len(ppi_src)} Genes:{n_genes}")
    return data,gene_to_idx,all_genes

# ── HeteroGAT ──

from torch_geometric.nn import GATConv, HeteroConv

class HeteroGAT(torch.nn.Module):
    def __init__(self, metadata, drug_dim, gene_dim, disease_dim,
                 hidden_dim, out_dim, heads, dropout):
        super().__init__()
        self.dropout_rate=dropout
        self.drug_proj=torch.nn.Linear(drug_dim,hidden_dim)
        self.gene_proj=torch.nn.Linear(gene_dim,hidden_dim)
        self.disease_proj=torch.nn.Linear(disease_dim,hidden_dim)

        convs1={}
        for etype in metadata[1]:
            convs1[etype]=GATConv((-1,-1),hidden_dim//heads,heads=heads,concat=True,
                                   dropout=dropout,add_self_loops=(etype[0]=="gene" and etype[2]=="gene"))
        self.conv1=HeteroConv(convs1,aggr="mean")

        convs2={}
        for etype in metadata[1]:
            convs2[etype]=GATConv((-1,-1),out_dim,heads=1,concat=False,
                                   dropout=dropout,add_self_loops=(etype[0]=="gene" and etype[2]=="gene"))
        self.conv2=HeteroConv(convs2,aggr="mean")

        self.dt_decoder=torch.nn.Bilinear(out_dim,out_dim,1)
        self.td_decoder=torch.nn.Bilinear(out_dim,out_dim,1)
        self.drug_out=torch.nn.Linear(hidden_dim,out_dim)

    def forward(self,x_dict,edge_index_dict):
        x_dict={k:v for k,v in x_dict.items()}
        proj={
            "drug":self.drug_proj(x_dict["drug"]),
            "gene":self.gene_proj(x_dict["gene"]),
            "disease":self.disease_proj(x_dict["disease"]),
        }
        out1=self.conv1(proj,edge_index_dict)
        for k in proj:
            if k not in out1: out1[k]=proj[k]
        out1={k:v.relu() for k,v in out1.items()}
        out1={k:F.dropout(v,p=self.dropout_rate,training=self.training) for k,v in out1.items()}
        out2=self.conv2(out1,edge_index_dict)
        for k in out1:
            if k not in out2: out2[k]=out1[k]
        out2["drug"]=self.drug_out(out2["drug"])
        return out2

# ── Negative Sampler ──

class NegSampler:
    def __init__(self,candidates,n_genes,seed=42):
        self.candidates=torch.from_numpy(candidates)
        self.rng=torch.Generator()
        self.rng.manual_seed(seed)
    def sample(self,pos_edge,is_dt):
        n=pos_edge.size(1)
        idx=torch.randint(0,len(self.candidates),(n,),generator=self.rng)
        neg=self.candidates[idx].to(pos_edge.device)
        if is_dt: return torch.stack([pos_edge[0],neg])
        else: return torch.stack([neg,pos_edge[1]])

# ── Train FULL Model (all edges) ──

def train_model(model,data,dt_edges,td_edges,dt_sampler,td_sampler,n_genes):
    opt=torch.optim.Adam(model.parameters(),lr=LR,weight_decay=WEIGHT_DECAY)
    best_loss,wait=float("inf"),0
    for ep in range(EPOCHS):
        model.train()
        loss_dt=torch.tensor(0.0,device=device)
        if dt_edges is not None:
            pos=dt_edges.to(device)
            neg=dt_sampler.sample(pos,is_dt=True)
            e=torch.cat([pos,neg],dim=1)
            lbl=torch.cat([torch.ones(pos.size(1)),torch.zeros(neg.size(1))]).to(device)
            z=model(data.x_dict,data.edge_index_dict)
            logits=model.dt_decoder(z["drug"][e[0]],z["gene"][e[1]]).squeeze(-1)
            loss_dt=F.binary_cross_entropy_with_logits(logits,lbl)

        loss_td=torch.tensor(0.0,device=device)
        if td_edges is not None:
            pos=td_edges.to(device)
            neg=td_sampler.sample(pos,is_dt=False)
            e=torch.cat([pos,neg],dim=1)
            lbl=torch.cat([torch.ones(pos.size(1)),torch.zeros(neg.size(1))]).to(device)
            z=model(data.x_dict,data.edge_index_dict)
            logits=model.td_decoder(z["gene"][e[0]],z["disease"][e[1]]).squeeze(-1)
            loss_td=F.binary_cross_entropy_with_logits(logits,lbl)

        loss=loss_dt+loss_td
        opt.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(),1.0)
        opt.step()

        if ep%50==0: print(f"  E{ep:3d}: {loss.item():.4f}")
        if loss.item()<best_loss-1e-4: best_loss=loss.item();wait=0
        else: wait+=1
        if wait>=PATIENCE: print(f"  Stop@{ep}");break

# ── Build negative candidates ──

def build_neg_candidates(ppi_ei,pos_indices,n_genes):
    ppi_src=ppi_ei[0].cpu().numpy();ppi_dst=ppi_ei[1].cpu().numpy()
    nb={i:set() for i in range(n_genes)}
    for s,d in zip(ppi_src,ppi_dst): nb[s].add(d);nb[d].add(s)
    cand=set()
    for idx in pos_indices: cand.update(nb.get(idx,set()))
    cand-=set(pos_indices)
    return np.array(sorted(cand),dtype=np.int64)

# ── Main ──

def main():
    print("="*60)
    print("  GAT Embedding v3 — FULL Model (All Edges)")
    print("  Single model trained on DT + TD + PPI simultaneously.")
    print("  Gene embeddings used for BOTH DT and DG classifiers.")
    print("="*60)

    dt_genes=load_drug_targets(os.path.join(GAT_DIR,"drug_targets.txt"))
    dg_genes=load_disease_genes(os.path.join(GAT_DIR,"disease_genes.txt"))
    ppi_edges=load_ppi(os.path.join(GAT_DIR,"ppi_subgraph.csv"))
    gene_feat,gene_names=load_gene_features(os.path.join(os.path.dirname(os.path.abspath(__file__)),"subgraph_embeddings.csv"))
    drug_fp=load_drug_fingerprint(os.path.join(GAT_DIR,"drug_fingerprint.csv"))
    all_genes=load_all_genes(os.path.join(GAT_DIR,"subgraph_genes.txt"))

    data,gene_to_idx,gene_list=build_graph(dt_genes,dg_genes,ppi_edges,gene_feat,gene_names,drug_fp,all_genes)
    n_genes=len(gene_list)
    data=data.to(device)

    full_dt=data["drug","targets","gene"].edge_index
    full_td=data["gene","associated_with","disease"].edge_index
    full_ppi=data["gene","interacts","gene"].edge_index

    dt_idx=[gene_to_idx[g] for g in dt_genes if g in gene_to_idx]
    dg_idx=[gene_to_idx[g] for g in dg_genes if g in gene_to_idx]
    dt_cand=build_neg_candidates(full_ppi,dt_idx,n_genes)
    dg_cand=build_neg_candidates(full_ppi,dg_idx,n_genes)

    dt_sampler=NegSampler(dt_cand,n_genes,RANDOM_SEED)
    td_sampler=NegSampler(dg_cand,n_genes,RANDOM_SEED)

    print("\n"+"="*50)
    print("  FULL GAT Model: DT + TD + PPI edges")
    print("="*50)

    set_seed(RANDOM_SEED)
    model=HeteroGAT(
        metadata=data.metadata(),
        drug_dim=data["drug"].x.size(1),gene_dim=data["gene"].x.size(1),
        disease_dim=data["disease"].x.size(1),hidden_dim=HIDDEN_DIM,
        out_dim=OUT_DIM,heads=GAT_HEADS,dropout=DROPOUT,
    ).to(device)

    train_model(model,data,full_dt,full_td,dt_sampler,td_sampler,n_genes)

    model.eval()
    with torch.no_grad():
        z=model(data.x_dict,data.edge_index_dict)
        emb=z["gene"].cpu().numpy()

    pca=PCA(n_components=N_PCA,random_state=RANDOM_SEED)
    emb_pca=pca.fit_transform(emb)
    var=pca.explained_variance_ratio_.sum()
    print(f"\n[PCA] Full model: 64d→{N_PCA}d, var={var:.3f}")

    cols=[f"gat_full_{i}" for i in range(N_PCA)]

    df_emb=pd.DataFrame(emb_pca,columns=cols,index=gene_list).reset_index()
    df_emb.columns=["gene_symbol"]+cols
    df_emb["gene_symbol"]=df_emb["gene_symbol"].str.upper()

    ft=pd.read_csv(FEATURE_TABLE)
    ft["gene_symbol"]=ft["gene_symbol"].str.upper()
    for c in cols:
        if c in ft.columns: ft=ft.drop(columns=[c])

    merged=ft.merge(df_emb,on="gene_symbol",how="left")
    for c in cols:
        merged[c]=merged[c].fillna(0.0)

    out_path=os.path.join(RF_DIR,"gene_features_table_with_gat_emb.csv")
    merged.to_csv(out_path,index=False,encoding="utf-8-sig")
    n_matched=(merged[cols].var(axis=1)>1e-8).sum()
    print(f"\n[Save] → {out_path}  ({merged.shape[0]}×{merged.shape[1]})")
    print(f"[Match] Full GAT: {n_matched}/{n_genes}")
    print("\n[Done] Full GAT embeddings ready for both classifiers.")

if __name__=="__main__":
    main()