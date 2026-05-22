# CIRI-cuproptosis-causal-discovery

Causal Discovery and Emergent Coarsening Reveal Cuproptosis Regulatory Modules in Cerebral Ischemia-Reperfusion Injury: A Multi-Omics Framework for Ethnomedicine Repositioning

## Overview

This project implements a full-chain causal inference framework to decode cuproptosis regulatory mechanisms in cerebral ischemia-reperfusion injury (CIRI) and identify potential ethnomedicine candidates.

### Pipeline Architecture

```
L1: Phenotype Anchoring → scRNA-seq + Bulk differential analysis
L2a: Causal Discovery → PC-NOTEARS DAG construction
L2b: Causal Emergence Coarsening → Greedy-Spectral Coarsening
L2c: Causal Intervention Simulation → IDA + Linear SEM
L3: Genetic Causal Anchoring → TwoSampleMR dual certification
L4: Ideal Intervention Subgraph → VGAE generation
L5: Reverse Screening → MCI + Hypergraph screening
L6: Validation Loop → Triple-anchor validation
```

## Directory Structure

```
├── L1_phenotype_anchoring/       # scRNA-seq + Bulk differential analysis
├── L2a_causal_discovery/         # PC-NOTEARS DAG construction
├── L2b_causal_coarsening/        # Greedy-Spectral Coarsening
├── L2c_causal_intervention/      # IDA + Linear SEM simulation
├── L3_genetic_causal/            # TwoSampleMR analysis
├── L4_ideal_subgraph/            # VGAE generation
├── L5_reverse_screening/         # MCI + Hypergraph screening
├── L6_validation/                # Triple-anchor validation
├── data/                         # Raw and processed data
├── results/                      # Analysis outputs
├── figures/                      # Publication figures
└── config/                       # Configuration files
```

## Quick Start

### Environment Setup

```bash
# Python environment
pip install -r requirements.txt

# R environment
Rscript -e "renv::restore()"
```

### Run L1 Analysis

```bash
# scRNA-seq analysis
python L1_phenotype_anchoring/scrna_analysis.py

# Bulk validation
Rscript L1_phenotype_anchoring/bulk_validation.R
```

## Datasets

| Dataset | Type | Platform | Source |
|---------|------|----------|--------|
| GSE174574 | scRNA-seq | 10x Genomics | GEO |
| GSE23160 | Bulk RNA-seq | Illumina microarray | GEO |
| eQTLGen | eQTL | Blood | IEU Open GWAS |
| GTEx v8 | eQTL | Brain tissues | GTEx Portal |
| TCMSP | Compound-Target | Database | TCMSP |
| ETCM | Compound-Target | Database | ETCM |

## Cuproptosis Gene Set

### Core Genes (10)
FDX1, LIAS, LIPT1, DLD, DLAT, PDHA1, PDHB, MTF1, GLS, CDKN2A

### Extended Genes (6)
SIRT7, ATP7B, CTR1 (SLC31A1), COX17, ATOX1, CCS

## Citation

If you use this code, please cite:

```
[Your paper citation will be added here]
```

## License

MIT License - see LICENSE file for details.

## Contact

For questions or collaboration, please open an issue or contact the authors.