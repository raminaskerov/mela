# Multi-Layer Transcriptomic Analysis of Encorafenib Resistance in BRAF-Mutant Melanoma

**Reproduction and extension of:** Çolako\u011flu Bergel et al. (2025), *Scientific Reports*  
**Dataset:** [GEO GSE283251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283251)  
**Analysis outputs:** [Zenodo — *DOI pending*]

---

## Overview

This repository provides an independent R/Bioconductor reimplementation of the transcriptomic analysis from Çolakoğlu Bergel et al. (2025), which investigated intrinsic resistance mechanisms to the BRAF inhibitor Encorafenib in A375 melanoma cell lines. The original study used a Galaxy-based workflow (Arga & Gülfidan, not publicly available); this repository offers a fully scripted, reproducible alternative and extends the analysis with several additional layers: co-expression network analysis (WGCNA), transcription factor activity scoring, multi-layer convergence, and protein–protein interaction hub analysis.

The central biological focus is the role of **iron metabolism, ferritinophagy, and lysosomal biology** as resistance mechanisms, with the miRNA–mRNA regulatory axis (particularly hsa-miR-140-3p → IREB2) as a key reference interaction from the source paper.

---

## Experimental Design

Four conditions from GSE283251 (A375 human melanoma cells):

| Label | Condition |
|-------|-----------|
| A375_SC | Sensitive cells, vehicle control |
| A375_RC | Resistant cells, vehicle control |
| A375_S10 | Sensitive cells + 10 nM Encorafenib |
| A375_R10 | Resistant cells + 10 nM Encorafenib |

> **Note:** n = 1 per condition (no biological replicates). NOISeq was selected as the differential expression method because it is specifically designed for no-replicate RNA-seq data. This single-replicate design sets a ceiling on statistical depth; analytical emphasis is placed on network topology, pathway-level patterns, and multi-layer convergence rather than individual gene-level inference.

Four comparisons analyzed: SC vs RC · SC vs S10 · RC vs R10 · S10 vs R10

---

## Analysis Pipeline

| # | Layer | Method | Package(s) |
|---|-------|--------|-----------|
| 1 | Differential expression (mRNA + miRNA) | NOISeq | `NOISeq` |
| 2 | Pathway enrichment | GSEA on full ranked gene list | `fgsea`, `clusterProfiler` |
| 3 | miRNA–mRNA regulatory networks | Pathway-guided reverse multiMiR query | `multiMiR` 1.32.0 |
| 4 | Co-expression modules | WGCNA | `WGCNA` |
| 5 | Transcription factor activity | CollecTRI regulon fgsea + decoupleR ULM | `decoupleR`, `fgsea` |
| 6 | Multi-layer convergence | 4-layer intersection scoring | custom R |
| 7 | PPI hub analysis | STRING-based degree + WGCNA/TF overlap | `in progress` |
| 8 | ML generalization | Cross-dataset feature classification | `planned` |

**Pathway databases used:** KEGG, Reactome, MSigDB Hallmark, CollecTRI (TF regulons)

---

## Key Findings (SC vs RC Comparison)

### Convergence analysis

Six genes reach all four analytical layers (GSEA Hallmark leading edge + WGCNA hub + miRNA network target + TF pathway convergence):

| Gene | Pathway context | Key miRNA regulators |
|------|----------------|----------------------|
| CHEK1 | DNA Damage Response / G2M | let-7 family (8 members, validated) |
| CDC6 | Cell Cycle / G2M | let-7 family, miR-93-5p, miR-25-3p (validated) |
| FAS | Apoptosis | hsa-miR-326 |
| HLA-A | Interferon Signaling | hsa-miR-744-5p, hsa-miR-1229-3p |
| ATP6V1C1 | Lysosomal / V-ATPase | hsa-miR-744-5p (validated), hsa-miR-326 |
| TCIRG1 | Lysosomal / V-ATPase | hsa-miR-744-5p |

Four additional genes reach three layers: TFRC, FTH1, TWIST1, ABL1.

### Two-arm mechanistic model for resistance

**Arm 1 — Cell cycle suppression:**  
let-7 family upregulation in RC → depletion of CHEK1 and CDC6 → E2F transcriptional suppression → loss of G2M checkpoint and homologous recombination capacity.

**Arm 2 — Lysosomal and immune activation:**  
Downregulation of miR-326 and miR-744-5p in RC → release of EMT regulators (TWIST1), interferon-stimulated genes (HLA-A), and lysosomal V-ATPase subunits (ATP6V1C1, TCIRG1) → coordinated by HIF1A, JUN, RELA, and STAT1 activity.

**Most novel extension relative to source paper:**  
Validated suppression of V-ATPase subunits (ATP6V1C1, TCIRG1) by miR-744-5p, directly connecting the miRNA regulatory axis to lysosomal ferritinophagy capacity in resistant cells — a mechanistic link not described in the original study.

---

## Repository Structure

```
.
├── scripts/
│   ├── 01_preprocessing.R           # GEO download, count matrix assembly
│   ├── 02_DEG_NOISeq.R              # Differential expression, 4 comparisons
│   ├── 03_GSEA_pathways.R           # KEGG, Reactome, Hallmark GSEA
│   ├── 04_miRNA_networks.R          # Pathway-guided miRNA–mRNA networks
│   ├── 05_WGCNA.R                   # Co-expression module construction
│   ├── 06_TF_enrichment.R           # CollecTRI fgsea + decoupleR ULM scoring
│   ├── 07_convergence.R             # 4-layer convergence + scoring
│   ├── 08_PPI_network.R             # PPI hub + TF/mRNA network  [in progress]
│   └── 09_ML_extension.R            # Multi-dataset ML classification [planned]
├── data/
│   └── README.md                    # GEO download instructions (raw data not included)
├── figures/                         # Key output figures
│   ├── gsea_hallmark_SC_vs_RC.pdf
│   ├── wgcna_module_trait_heatmap.pdf
│   ├── miRNA_network_arm1_arm2.pdf
│   └── four_layer_convergence.pdf
├── outputs/                         # Full results tables — see Data Availability
│   └── .gitkeep
└── README.md
```

---

## Planned Extensions

### PPI / TF–mRNA Network *(in progress)*
Integration of STRING protein interaction data with WGCNA hub genes and TF candidate list to identify high-connectivity resistance hubs. Outputs include annotated hub table and Cytoscape-compatible network files (CX format).
Results got but no verification and analysis at biological level.

### ML Generalization *(planned)*
Multi-dataset machine learning to assess whether the resistance gene signature generalizes across cell line models:
- GSE283251 (A375, primary dataset)
- GSE202118 (M229 parental/resistant)
- GSE114443 (8 monolayer samples)
- PRJNA602782, GSE45558 (in queue)
- GSE148638 (processed counts, pending verification)

Goal: feature importance analysis to identify the most robust cross-dataset resistance markers.

---

## Setup & Dependencies

R ≥ 4.3, Bioconductor 3.20

```r
# Install Bioconductor packages
BiocManager::install(c(
  "NOISeq", "clusterProfiler", "fgsea", "WGCNA",
  "multiMiR", "decoupleR", "viper", "OmnipathR",
  "org.Hs.eg.db", "ReactomePA"
))
```

Visualization: `ggplot2`, `pheatmap`, `Cairo` (Windows Unicode fix for PDF output)  
Network export: Cytoscape (CX-format files)

---

## Data Availability

**Raw counts:** GEO [GSE283251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283251) — download instructions in `data/README.md`  
**Analysis outputs** (DEG tables, GSEA results, network edges, WGCNA module files): [Zenodo — *DOI pending after deposit*]

---

## Reference

Çolako\u011flu Bergel N, et al. (2025). *[Full title]*. *Scientific Reports*. DOI: *[add when available]*

---

## Notes on Reproduction

- NOISeq computes M = log2(A/B), so for the SC vs RC comparison, **positive M = higher in SC**. Sign corrected then.
- Histone gene filter applied (replication-dependent HIST1/HIST2 clusters removed before enrichment).
- BH correction is not applied to NOISeq probability scores; `prob ≥ 0.8` threshold used directly, consistent with NOISeq documentation.
- `fgsea` ties warning (~85% tied ranks) is expected from NOISeq's integer count ratios and does not indicate a pipeline error; S10 vs R10 GSEA results are the least reliable comparison for this reason.
