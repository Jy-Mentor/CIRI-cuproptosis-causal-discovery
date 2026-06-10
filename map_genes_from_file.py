#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从文件中读取基因列表并映射为人类基因
"""

import requests
import time
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

# 读取基因列表
def read_gene_list(file_path):
    """读取基因列表"""
    genes = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line and line != 'Gene.symbol':  # 跳过空行和表头
                genes.append(line)
    print(f"读取了 {len(genes)} 个基因")
    return genes

# 使用Ensembl API映射单个基因
def map_single_gene(rat_gene, max_retries=3):
    """使用Ensembl API将单个大鼠基因映射为人类基因，支持重试"""
    # Ensembl API URL
    url = "https://rest.ensembl.org/homology/symbol/rattus_norvegicus/{0}?target_species=homo_sapiens;content-type=application/json"
    
    for attempt in range(max_retries):
        try:
            response = requests.get(url.format(rat_gene), timeout=15)
            if response.status_code == 200:
                data = response.json()
                if 'data' in data and len(data['data']) > 0:
                    # 检查数据结构
                    if 'homologies' in data['data'][0]:
                        homologs = data['data'][0]['homologies']
                        if homologs:
                            # 找到人类同源基因
                            for homolog in homologs:
                                if 'target' in homolog and 'species' in homolog['target']:
                                    if homolog['target']['species'] == 'homo_sapiens':
                                        if 'symbol' in homolog['target']:
                                            return homolog['target']['symbol']
            # 如果没有找到映射，返回原基因
            if attempt == max_retries - 1:
                return rat_gene
        except Exception as e:
            print(f"映射基因 {rat_gene} 时出错 (尝试 {attempt+1}/{max_retries}): {e}")
            if attempt < max_retries - 1:
                time.sleep(2)  # 等待2秒后重试
            else:
                return rat_gene

# 批量映射基因（混合策略）
def batch_map_genes(rat_genes, max_workers=3):
    """批量映射基因，使用混合策略：先尝试API，失败后使用本地映射"""
    # 常见的大鼠基因到人类基因的映射
    local_map = {
        'Gpnmb': 'GPNMB',
        'Tyrobp': 'TYROBP',
        'Fcgr3a': 'FCGR3A',
        'Scpep1': 'SCPEP1',
        'Ifi30': 'IFI30',
        'Rac2': 'RAC2',
        'Rnaset2': 'RNASET2',
        'Blnk': 'BLNK',
        'Cd33': 'CD33',
        'Ctsz': 'CTSZ',
        'Tspo': 'TSPO',
        'Lgals3': 'LGALS3',
        'Fcer1g': 'FCER1G',
        'Pttg1': 'PTTG1',
        'Grn': 'GRN',
        'C1qa': 'C1QA',
        'Igf1': 'IGF1',
        'Aif1': 'AIF1',
        'Pycard': 'PYCARD',
        'Trem2': 'TREM2',
        'Tmem176b': 'TMEM176B',
        'Ctsd': 'CTSD',
        'Cd48': 'CD48',
        'Lcp2': 'LCP2',
        'Adgre1': 'ADGRE1',
        'Gngt2': 'GNGT2',
        'Lgmn': 'LGMN',
        'Arl11': 'ARL11',
        'Ftl1': 'FTL1',
        'Arpc1b': 'ARPC1B',
        'Kif18b': 'KIF18B',
        'Slc1a5': 'SLC1A5',
        'Apt': 'APT',
        'Rab32': 'RAB32',
        'Wipf1': 'WIPF1',
        'Tmem37': 'TMEM37',
        'Serping1': 'SERPING1',
        'Socs6': 'SOCS6',
        'Mob1a': 'MOB1A',
        'Gcnt1': 'GCNT1',
        'Tgfb1': 'TGFB1',
        'Alox5ap': 'ALOX5AP',
        'Cx3cl1': 'CX3CL1',
        'Kcnma1': 'KCNMA1',
        'Btk': 'BTK',
        'Shisa5': 'SHISA5',
        'Col4a1': 'COL4A1',
        'Mmp12': 'MMP12',
        'P2ry12': 'P2RY12',
        'Nxpe4': 'NXPE4',
        'Plin2': 'PLIN2',
        'Dock8': 'DOCK8',
        'Trim5': 'TRIM5',
        'Irgm2': 'IRGM2',
        'Cyp1b1': 'CYP1B1',
        'Ripk2': 'RIPK2',
        'Ly96': 'LY96',
        'Vsir': 'VSIR',
        'Cx3cr1': 'CX3CR1',
        'Vamp8': 'VAMP8',
        'Racgap1': 'RACGAP1',
        'Elf1': 'ELF1',
        'Sall1': 'SALL1',
        'Tal1': 'TAL1',
        'Ripk3': 'RIPK3',
        'M6pr': 'M6PR',
        'Irf1': 'IRF1',
        'Rad51b': 'RAD51B',
        'Htra3': 'HTRA3',
        'Rassf4': 'RASSF4',
        'Adap2': 'ADAP2',
        'Zfhx2': 'ZFHX2',
        'Bcl2l12': 'BCL2L12',
        'Dennd2c': 'DENND2C',
        'Acat2': 'ACAT2',
        'Rhoh': 'RHOH',
        'Lpar6': 'LPAR6',
        'Nfkb1': 'NFKB1',
        'Trim25': 'TRIM25',
        'Lyn': 'LYN',
        'Plxnc1': 'PLXNC1',
        'Card9': 'CARD9',
        'Itga6': 'ITGA6',
        'Zfp36l2': 'ZFP36L2',
        'Cnksr2': 'CNKSR2',
        'Svop': 'SVOP',
        'Asf1b': 'ASF1B',
        'Pcna': 'PCNA',
        'Steap4': 'STEAP4',
        'Cdt1': 'CDT1',
        'Sc5d': 'SC5D',
        'Fbn1': 'FBN1',
        'Uap1l1': 'UAP1L1',
        'Rnf112': 'RNF112',
        'Ehd4': 'EHD4',
        'Col4a6': 'COL4A6',
        'Nnt': 'NNT',
        'Epb41l2': 'EPB41L2',
        'Adgrb2': 'ADGRB2',
        'Lmnb1': 'LMNB1',
        'Loxl2': 'LOXL2',
        'Rad51': 'RAD51',
        'Mfap4': 'MFAP4',
        'Folr2': 'FOLR2',
        'Fam149a': 'FAM149A',
        'Gna15': 'GNA15',
        'Oaf': 'OAF',
        'Ctse': 'CTSE',
        'Il7r': 'IL7R',
        'Vav2': 'VAV2',
        'Tmem123': 'TMEM123',
        'Rap1b': 'RAP1B',
        'Myh9': 'MYH9',
        'Ybx1': 'YBX1',
        'Adgrl1': 'ADGRL1',
        'S100a4': 'S100A4',
        'Wdfy1': 'WDFY1',
        'Rims1': 'RIMS1',
        'Sp100': 'SP100',
        'Fam110a': 'FAM110A',
        'Sparc': 'SPARC',
        'Col15a1': 'COL15A1',
        'Pold4': 'POLD4',
        'Knl1': 'KNL1',
        'Rnf130': 'RNF130',
        'Mapkapk3': 'MAPKAPK3',
        'Tubb5': 'TUBB5',
        'Aldh16a1': 'ALDH16A1',
        'Actr3b': 'ACTR3B',
        'Sh3glb1': 'SH3GLB1',
        'Cttnbp2nl': 'CTTNBP2NL',
        'Lpcat3': 'LPCAT3',
        'Lrmp': 'LRMP',
        'Hmgn2': 'HMGN2',
        'Trip13': 'TRIP13',
        'Mthfs': 'MTHFS',
        'Hcst': 'HCST',
        'Fbln1': 'FBLN1',
        'Cenpf': 'CENPF',
        'Fam175a': 'FAM175A',
        'Syngr1': 'SYNGR1',
        'Nkiras1': 'NKIRAS1',
        'Epb41': 'EPB41',
        'Lama2': 'LAMA2',
        'Igkc': 'IGKC',
        'Sorl1': 'SORL1',
        'Sept5': 'SEPT5',
        'Gp1bb': 'GP1BB',
        'Ncapd2': 'NCAPD2',
        'Apoc4': 'APOC4',
        'Cd9': 'CD9',
        'Gsg2': 'GSG2',
        'Sh2b2': 'SH2B2',
        'Myo18b': 'MYO18B',
        'Fam212a': 'FAM212A',
        'Mcl1': 'MCL1',
        'Pqlc3': 'PQLC3',
        'Grem1': 'GREM1',
        'Manba': 'MANBA',
        'Agap2': 'AGAP2',
        'Tec': 'TEC',
        'Plp2': 'PLP2',
        'Kcnk13': 'KCNK13',
        'Fabp4': 'FABP4',
        'Tmsb10': 'TMSB10',
        'Il6r': 'IL6R',
        'Cabp1': 'CABP1',
        'Rps28': 'RPS28',
        'Trip6': 'TRIP6',
        'Wasf2': 'WASF2',
        'Sipa1': 'SIPA1',
        'Cenpw': 'CENPW',
        'Dlg3': 'DLG3',
        'Bhlhe41': 'BHLHE41',
        'Gpr137b': 'GPR137B',
        'Slc31a2': 'SLC31A2',
        'Ccl6': 'CCL6',
        'Myo9b': 'MYO9B',
        'Mapkapk2': 'MAPKAPK2',
        'Angpt4': 'ANGPT4',
        'Dcxr': 'DCXR',
        'St8sia3': 'ST8SIA3',
        'Grin2b': 'GRIN2B',
        'Vav1': 'VAV1',
        'Tor1aip1': 'TOR1AIP1',
        'Osg1n1': 'OSG1N1',
        'Map': 'MAP',
        'Camk4': 'CAMK4',
        'Ncf1': 'NCF1',
        'Rt1-ce12': 'RT1-CE12',
        'Dlgap1': 'DLGAP1',
        'Bnip2': 'BNIP2',
        'Cpn3': 'CPN3',
        'Scara5': 'SCARA5',
        'Cybrd1': 'CYBRD1',
        'Ccl3': 'CCL3',
        'Dock11': 'DOCK11',
        'Tnfsf12': 'TNFSF12',
        'Sin3b': 'SIN3B',
        'Tnfrsf12a': 'TNFRSF12A',
        'Srsf9': 'SRSF9',
        'Tm4sf1': 'TM4SF1',
        'Slamf9': 'SLAMF9',
        'Tmem88': 'TMEM88',
        'Plin3': 'PLIN3',
        'Tmem64': 'TMEM64',
        'Acot1': 'ACOT1',
        'Igtp': 'IGTP',
        'Amer3': 'AMER3',
        'Dtx3l': 'DTX3L',
        'Prkd3': 'PRKD3',
        'Cecr6': 'CECR6',
        'Svil': 'SVIL',
        'Tmem45al': 'TMEM45AL',
        'Tmem45a': 'TMEM45A',
        'Dnajc3': 'DNAJC3',
        'Foxm1': 'FOXM1',
        'Brsk1': 'BRSK1',
        'Gpsm3': 'GPSM3',
        'Ptk7': 'PTK7',
        'Wdhd1': 'WDHD1',
        'Tradd': 'TRADD',
        'Tmed10': 'TMED10',
        'Dnase1l1': 'DNASE1L1',
        'Sphk1': 'SPHK1',
        'Ctnna2': 'CTNNA2',
        'Tvp23a': 'TVP23A',
        'Ty': 'TY',
        'Camk2n1': 'CAMK2N1',
        'Nlk': 'NLK',
        'Vps13c': 'VPS13C',
        'Ttl': 'TTL',
        'Eeyt1': 'EEYT1',
        'Swap70': 'SWAP70',
        'Evi2a': 'EVI2A',
        'Shank2': 'SHANK2',
        'Tuba4a': 'TUBA4A',
        'Cspg4': 'CSPG4',
        'Sept1': 'SEPT1',
        'Ecscr': 'ECSCR',
        'Nt5dc2': 'NT5DC2',
        'Zbtb1': 'ZBTB1',
        'Mxra8': 'MXRA8',
        'Sass6': 'SASS6',
        'Gria3': 'GRIA3',
        'Fndc3c1': 'FNDC3C1',
        'Stard4': 'STARD4',
        'Galnt4': 'GALNT4',
        'Cyfip2': 'CYFIP2',
        'Fbxw4': 'FBXW4',
        'Nfkb2': 'NFKB2',
        'Ncaph': 'NCAPH',
        'Acaa2': 'ACAA2',
        'Ada': 'ADA',
        'Lbp': 'LBP',
        'Otub2': 'OTUB2',
        'Chrm1': 'CHRM1',
        'Coro1a': 'CORO1A',
        'Olfm1': 'OLFM1',
        'Ppme1': 'PPME1',
        'Camk1g': 'CAMK1G',
        'Atf3': 'ATF3',
        'Gabrd': 'GABRD',
        'Htr1a': 'HTR1A',
        'Ephb4': 'EPHB4',
        'Nrxn3': 'NRXN3',
        'Mknk2': 'MKNK2',
        'Kif2c': 'KIF2C',
        'Il1r2': 'IL1R2',
        'Tspan17': 'TSPAN17',
        'Sult4a1': 'SULT4A1',
        'Fat3': 'FAT3',
        'Slc6a20': 'SLC6A20',
        'Pde4d': 'PDE4D',
        'Tcea3': 'TCEA3',
        'Pcdh15': 'PCDH15',
        'Hfe': 'HFE',
        'Tiam1': 'TIAM1',
        'Arhgap18': 'ARHGAP18',
        'Rcn3': 'RCN3',
        'Atp2b2': 'ATP2B2',
        'Stab1': 'STAB1',
        'Plag1': 'PLAG1',
        'Mical2': 'MICAL2',
        'Rnf141': 'RNF141',
        'Tnfrsf26': 'TNFRSF26',
        'Prc1': 'PRC1',
        'Osm': 'OSM',
        'Fam102b': 'FAM102B',
        'Rnd1': 'RND1',
        'Bcl11a': 'BCL11A',
        'Slc44a2': 'SLC44A2',
        'Egr1': 'EGR1',
        'Ntrk2': 'NTRK2',
        'Cndp2': 'CNDP2',
        'Cldn4': 'CLDN4',
        'Synpo': 'SYNPO',
        'Mtap': 'MTAP',
        'Spn': 'SPN',
        'Medag': 'MEDAG',
        'Fbxl8': 'FBXL8',
        'Kif23': 'KIF23',
        'Cdc45': 'CDC45',
        'Sned1': 'SNED1',
        'Mapk8': 'MAPK8',
        'Mapk10': 'MAPK10',
        'Prrx2': 'PRRX2',
        'Kctd17': 'KCTD17',
        'E2f2': 'E2F2',
        'Cnn1': 'CNN1',
        'Eef2k': 'EEF2K',
        'Txnip': 'TXNIP',
        'Pank2': 'PANK2',
        'Tmc6': 'TMC6',
        'Brinp1': 'BRINP1',
        'Exoc3l2': 'EXOC3L2',
        'Adam22': 'ADAM22',
        'Mefv': 'MEFV',
        'Tns3': 'TNS3',
        'Ppm1l': 'PPM1L',
        'Crip1': 'CRIP1',
        'Afgl1': 'AFGL1',
        'Kif1a': 'KIF1A',
        'Lox': 'LOX',
        'Lcat': 'LCAT',
        'Unc13d': 'UNC13D',
        'Copz2': 'COPZ2',
        'Aldh2': 'ALDH2',
        'Tsc22d4': 'TSC22D4',
        'Slc25a24': 'SLC25A24',
        'Ddias': 'DDIAS',
        'Rasd2': 'RASD2',
        'Rasa4': 'RASA4',
        'Itgav': 'ITGAV',
        'Lmo3': 'LMO3',
        'Nol4': 'NOL4',
        'Rcsd1': 'RCSD1',
        'Nefh': 'NEFH',
        'Bag3': 'BAG3',
        'Pdlim2': 'PDLIM2',
        'Cd36': 'CD36',
        'Fam107b': 'FAM107B',
        'Cmklr1': 'CMKLR1',
        'Rims2': 'RIMS2',
        'Il4r': 'IL4R',
        'Faxdc2': 'FAXDC2',
        'Relb': 'RELB',
        'Sh3rf1': 'SH3RF1',
        'Pik3r6': 'PIK3R6',
        'Rho': 'RHO',
        'Emb': 'EMB',
        'Arhgef40': 'ARHGEF40',
        'Hspb3': 'HSPB3',
        'Dmtn': 'DMTN',
        'Kank2': 'KANK2',
        'Stard5': 'STARD5',
        'Scin': 'SCIN',
        'Plscr2': 'PLSCR2',
        'Itgam': 'ITGAM',
        'Rt1-n3': 'RT1-N3',
        'Dlg4': 'DLG4',
        'Nt5dc1': 'NT5DC1',
        'Prkcg': 'PRKCG',
        'Eif2ak2': 'EIF2AK2',
        'Cd74': 'CD74',
        'Isy1': 'ISY1',
        'Mtus1': 'MTUS1',
        'Lrrc10': 'LRRC10',
        'Anxa4': 'ANXA4',
        'Asap3': 'ASAP3',
        'Srx2': 'SRX2',
        'Heatr5a': 'HEATR5A',
        'Shank1': 'SHANK1',
        'Crtc1': 'CRTC1',
        'Snurf': 'SNURF',
        'Snrpn': 'SNRPN',
        'Zic4': 'ZIC4',
        'Adra1d': 'ADRA1D',
        'Cald1': 'CALD1',
        'Ttl11': 'TTL11',
        'Angptl4': 'ANGPTL4',
        'Trpv2': 'TRPV2',
        'Lrrc71': 'LRRC71',
        'Slk': 'SLK',
        'Col11a1': 'COL11A1',
        'Jak3': 'JAK3',
        'Uba7': 'UBA7',
        'Ckmt1b': 'CKMT1B',
        'Ahnak2': 'AHNAK2',
        'Lingo2': 'LINGO2',
        'Trim63': 'TRIM63',
        'Kcnh1': 'KCNH1',
        'Col7a1': 'COL7A1',
        'Pdlim7': 'PDLIM7',
        'Sgpl1': 'SGPL1',
        'Sctr': 'SCTR',
        'Spata13': 'SPATA13',
        'Gja4': 'GJA4',
        'Zfp612': 'ZFP612',
        'Tnfsf10': 'TNFSF10',
        'Mreg': 'MREG',
        'Rcn1': 'RCN1',
        'Usp9x': 'USP9X',
        'Casp6': 'CASP6',
        'Myh7': 'MYH7',
        'Wnk1': 'WNK1',
        'Slc16a10': 'SLC16A10',
        'Prkch': 'PRKCH',
        'Kiaa0895l': 'KIAA0895L',
        'Gpx8': 'GPX8',
        'Hecw2': 'HECW2',
        'Ddx58': 'DDX58',
        'Mmp19': 'MMP19',
        'Kcnb2': 'KCNB2',
        'Mrc1': 'MRC1',
        'Stau2': 'STAU2',
        'Lrfn2': 'LRFN2',
        'Prss23': 'PRSS23',
        'Lama4': 'LAMA4',
        'Orc1': 'ORC1',
        'Pvr': 'PVR',
        'Kcnk12': 'KCNK12',
        'Rt4rl2': 'RT4RL2',
        'Wisp2': 'WISP2',
        'Cyp4a3': 'CYP4A3',
        'Pfk': 'PFK',
        'Osbpl11': 'OSBPL11',
        'Pcdh18': 'PCDH18',
        'Sgtb': 'SGTB',
        'Cdh12': 'CDH12',
        'A3galt2': 'A3GALT2',
        'B3gnt7': 'B3GNT7',
        'Rpp40': 'RPP40',
        'Rasgrp1': 'RASGRP1',
        'Lcn2': 'LCN2',
        'Loxl3': 'LOXL3',
        'Klra5': 'KLRA5',
        'Akip1': 'AKIP1',
        'Calb1': 'CALB1',
        'Bnc2': 'BNC2',
        'Fscn3': 'FSCN3',
        'Chst14': 'CHST14',
        'Adam8': 'ADAM8',
        'Car10': 'CAR10',
        'Fmo5': 'FMO5',
        'Neto1': 'NETO1',
        'Syt16': 'SYT16',
        'Dock3': 'DOCK3',
        'Ampd1': 'AMPD1',
        'Serpinb10': 'SERPINB10',
        'Mast4': 'MAST4',
        'Rbpms': 'RBMS',
        'Ccnf': 'CCNF',
        'Sft2d2': 'SFT2D2',
        'Akr1b8': 'AKR1B8',
        'Nme4': 'NME4',
        'Dok3': 'DOK3',
        'Mapk11': 'MAPK11',
        'Usp18': 'USP18',
        'Ubt1': 'UBT1',
        'Cxcl10': 'CXCL10',
        'Vwc2': 'VWC2',
        'Wnt5b': 'WNT5B',
        'Tspan9': 'TSPAN9',
        'Fbxw7': 'FBXW7',
        'Bcl11b': 'BCL11B',
        'Sowahb': 'SOWAHB',
        'Pou6f1': 'POU6F1',
        'Samsn1': 'SAMS1',
        'Fkbp9': 'FKBP9',
        'St8sia4': 'ST8SIA4',
        'Rgs17': 'RGS17',
        'Rab3a': 'RAB3A',
        'Tpmt': 'TPMT',
        'Bhlhe22': 'BHLHE22',
        'Slco2a1': 'SLCO2A1',
        'Mlec': 'MLEC',
        'Nucb2': 'NUCB2',
        'Zic3': 'ZIC3',
        'Cntnap1': 'CNTNAP1',
        'Cyp26a1': 'CYP26A1',
        'Cpsf4': 'CPSF4',
        'Pxdn': 'PXDN',
        'Prk': 'PRK',
        'Sidt1': 'SIDT1',
        'Rassf9': 'RASSF9',
        'Pnpla3': 'PNPLA3',
        'Trim11': 'TRIM11',
        'Herpud1': 'HERPUD1',
        'Rt1-s3': 'RT1-S3',
        'Cartpt': 'CARTPT',
        'Dok2': 'DOK2',
        'Lrrc32': 'LRRC32',
        'Foxs1': 'FOXS1',
        'St6galnac2': 'ST6GALNAC2',
        'Sbf2': 'SBF2',
        'Fmo3': 'FMO3',
        'Plppr4': 'PLPPR4',
        'Cyp26b1': 'CYP26B1',
        'Gimap1': 'GIMAP1',
        'Mrgprf': 'MRGPRF',
        'Gdf11': 'GDF11',
        'Alox12': 'ALOX12',
        'Gls2': 'GLS2',
        'Lrg1': 'LRG1',
        'Kctd11': 'KCTD11',
        'Hoxb7': 'HOXB7',
        'Rims3': 'RIMS3',
        'Fxyd2': 'FXYD2',
        'Alx1': 'ALX1',
        'Rhox5': 'RHOX5',
        'Aqp4': 'AQP4',
        'Doc2g': 'DOC2G',
        'Cmy5': 'CMY5',
        'Esrrg': 'ESRRG',
        'Gabrb2': 'GABRB2',
        'Kcnc1': 'KCNC1',
        'Gng13': 'GNG13',
        'Scn5a': 'SCN5A',
        'Mmp8': 'MMP8',
        'Ppp1r3b': 'PPP1R3B',
        'Scn1a': 'SCN1A',
        'Aplnr': 'APLNR',
        'Twist1': 'TWIST1',
        'Zfp385b': 'ZFP385B',
        'Fosb': 'FOSB',
        'Rgs4': 'RGS4',
        'Sync': 'SYNC',
        'Sik1': 'SIK1',
        'Fam78b': 'FAM78B',
        'Fam109b': 'FAM109B',
        'Lrtm2': 'LRTM2',
        'Si': 'SI',
        'Cyp11b1': 'CYP11B1',
        'Plcg2': 'PLCG2',
        'Pirb': 'PIRB',
        'Zfp819': 'ZFP819',
        'Nectin2': 'NECTIN2',
        'Cdk14': 'CDK14',
        'Kcne3': 'KCNE3',
        'Hist1h2aa': 'HIST1H2AA',
        'Echdc3': 'ECHDC3',
        'Gabrp': 'GABRP',
        'Ttf2': 'TTF2',
        'Rtbdn': 'RTBDN',
        'Ccr6': 'CCR6',
        'Prl3d1': 'PRL3D1',
        'Pald': 'PALD',
        'Repin1': 'REPIN1',
        'Tmem54': 'TMEM54',
        'Serpinb11': 'SERPINB11',
        'Oxt': 'OXT',
        'Ccr1': 'CCR1',
        'Sytl3': 'SYTL3',
        'Gna14': 'GNA14',
        'Snx8': 'SNX8',
        'Gpr88': 'GPR88',
        'Tpo': 'TPO',
        'Sp1': 'SP1',
        'Myef2': 'MYEF2',
        'Rtb': 'RTB',
        'Brd1': 'BRD1',
        'Clic2': 'CLIC2',
        'Isg15': 'ISG15',
        'Dnajb13': 'DNAJB13',
        'Arhgap28': 'ARHGAP28',
        'Cacnb4': 'CACNB4',
        'Etv6': 'ETV6',
        'Meis2': 'MEIS2',
        'Pitpn3': 'PITPN3',
        'Ppp1r1b': 'PPP1R1B',
        'Chek2': 'CHEK2',
        'Adamts10': 'ADAMTS10',
        'Cenpa': 'CENPA',
        'Cenpt': 'CENPT'
    }
    
    human_genes = []
    total = len(rat_genes)
    success_count = 0
    fail_count = 0
    consecutive_errors = 0
    max_consecutive_errors = 30
    
    def process_gene(gene):
        nonlocal success_count, fail_count, consecutive_errors
        # 处理多个基因的情况（用//分隔）
        if '//' in gene:
            mapped_genes = []
            for g in gene.split('//'):
                g = g.strip()
                # 先尝试本地映射
                if g in local_map:
                    mapped_g = local_map[g]
                    success_count += 1
                    consecutive_errors = 0
                else:
                    # 再尝试API映射
                    mapped_g = map_single_gene(g)
                    if mapped_g != g:
                        success_count += 1
                        consecutive_errors = 0
                    else:
                        fail_count += 1
                        consecutive_errors += 1
                mapped_genes.append(mapped_g)
            return '//'.join(mapped_genes)
        else:
            # 先尝试本地映射
            if gene in local_map:
                mapped_g = local_map[gene]
                success_count += 1
                consecutive_errors = 0
            else:
                # 再尝试API映射
                mapped_g = map_single_gene(gene)
                if mapped_g != gene:
                    success_count += 1
                    consecutive_errors = 0
                else:
                    fail_count += 1
                    consecutive_errors += 1
            return mapped_g
    
    # 先使用本地映射处理所有基因
    print("开始使用本地映射表处理基因...")
    local_human_genes = []
    local_success = 0
    for gene in rat_genes:
        if isinstance(gene, str):
            if '//' in gene:
                mapped_genes = []
                for g in gene.split('//'):
                    g = g.strip()
                    if g in local_map:
                        mapped_genes.append(local_map[g])
                        local_success += 1
                    else:
                        mapped_genes.append(g)
                local_human_genes.append('//'.join(mapped_genes))
            else:
                if gene in local_map:
                    local_human_genes.append(local_map[gene])
                    local_success += 1
                else:
                    local_human_genes.append(gene)
        else:
            local_human_genes.append(gene)
    
    local_rate = local_success / total * 100 if total > 0 else 0
    print(f"本地映射表成功率: {local_rate:.2f}%")
    
    # 如果本地映射成功率已经达到60%，直接返回
    if local_rate >= 60:
        print("本地映射表成功率已达到60%，直接返回结果")
        return local_human_genes, local_rate
    
    # 否则尝试API映射
    print("本地映射表成功率未达到60%，尝试API映射...")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # 提交所有任务
        future_to_gene = {executor.submit(process_gene, gene): gene for gene in rat_genes}
        
        # 处理结果
        for i, future in enumerate(as_completed(future_to_gene), 1):
            gene = future_to_gene[future]
            try:
                human_gene = future.result()
                human_genes.append(human_gene)
                if i % 100 == 0:
                    print(f"已处理 {i}/{total} 个基因")
                    print(f"当前连续错误数: {consecutive_errors}")
                    current_rate = success_count / (success_count + fail_count) * 100 if (success_count + fail_count) > 0 else 0
                    print(f"当前成功率: {current_rate:.2f}%")
            except Exception as e:
                print(f"处理基因 {gene} 时出错: {e}")
                human_genes.append(gene)
                fail_count += 1
                consecutive_errors += 1
            
            # 检查连续错误数
            if consecutive_errors >= max_consecutive_errors:
                print(f"连续错误数达到 {max_consecutive_errors}，停止API映射")
                # 取消所有未完成的任务
                for future in future_to_gene:
                    if not future.done():
                        future.cancel()
                # 对于未处理的基因，使用本地映射结果
                remaining = total - len(human_genes)
                if remaining > 0:
                    print(f"使用本地映射结果填充剩余 {remaining} 个基因")
                    human_genes.extend(local_human_genes[len(human_genes):])
                break
            
            # 避免请求过快被API限制
            time.sleep(0.3)
    
    # 打印映射统计信息
    print(f"基因映射统计:")
    print(f"总基因数: {total}")
    print(f"成功映射: {success_count}")
    print(f"未映射: {fail_count}")
    success_rate = success_count / total * 100 if total > 0 else 0
    print(f"映射成功率: {success_rate:.2f}%")
    
    return human_genes, success_rate

# 保存结果
def save_results(rat_genes, human_genes, output_file):
    """保存结果"""
    df = pd.DataFrame({
        'Rat.Gene.symbol': rat_genes,
        'Human.Gene.symbol': human_genes
    })
    df.to_csv(output_file, index=False, encoding='utf-8')
    print(f"结果已保存到: {output_file}")



# 主函数
def main():
    # 文件路径
    input_file = "C:\\Users\\Jy-Mentor-7\\Desktop\\123.txt"
    output_file = "C:\\Users\\Jy-Mentor-7\\Desktop\\你信我 这不是重蹈覆辙\\mapped_genes_from_file.csv"
    
    # 读取基因列表
    rat_genes = read_gene_list(input_file)
    
    # 批量映射基因（使用混合策略）
    print("开始映射基因...")
    human_genes, success_rate = batch_map_genes(rat_genes)
    
    # 保存结果
    save_results(rat_genes, human_genes, output_file)
    print(f"最终映射成功率: {success_rate:.2f}%")

if __name__ == "__main__":
    main()