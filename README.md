# Multi-Layer Transcriptomic Analysis of Encorafenib Resistance in BRAF-Mutant Melanoma

**Reproduction and extension of:** Çolako\u011flu Bergel et al. (2025), *Scientific Reports*  
**Dataset:** [GEO GSE283251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283251)  
**Analysis outputs:** [Zenodo — *DOI pending*]

---

## Overview

This repository provides an independent R/Bioconductor reimplementation of the transcriptomic analysis from Çolakoğlu Bergel et al. (2025), which investigated intrinsic resistance mechanisms to the BRAF inhibitor Encorafenib in A375 melanoma cell lines. The original study used a Galaxy-based workflow (Arga & Gülfidan, not publicly available); this repository offers a fully scripted, reproducible alternative and extends the analysis with several additional layers: pathway enrichment, miRNA–mRNA network analysis, co-expression network analysis (WGCNA), transcription factor enrichment, and an optional legacy TF/WGCNA linkage stage.

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

The current R scripts are organized as a numbered checkpoint chain:

| Stage | Script | Input | Output |
|---|---|---|---|
| 01 | `01_preprocessing.R` | GEO count files | `checkpoint_01.RData` |
| 02 | `02_DEG_NOISeq.R` | `checkpoint_01.RData` | DEG tables + `checkpoint_02.RData` |
| 03 | `03_enrichment.R` | `checkpoint_02.RData` | pathway GSEA tables + `checkpoint_03.RData` |
| 04 | `04_visualizations.R` | `checkpoint_03.RData` | volcano plots, heatmaps, summary figures |
| 05 | `05_mirna_mrna_integration.R` | `checkpoint_03.RData` + `pathway_gene_df` | reverse multiMiR network outputs + `checkpoint_05.RData` |
| 06 | `06_TF_enrichment.R` | `checkpoint_03.RData` | TF fgsea, TF activity scores + `checkpoint_06.RData` |
| 07 | `07_WGCNA.R` | `checkpoint_03.RData` | WGCNA modules + `checkpoint_07.RData` |
| 08 | `08_legacy_wgcna_tf_linkage.R` | `checkpoint_06.RData` + `checkpoint_07.RData` | legacy TF/WGCNA linkage tables + `checkpoint_08.RData` |

`main pipeline 2.R` is now only a thin runner that sources these stage scripts in order.

**Pathway databases used:** KEGG, Reactome, MSigDB Hallmark, CollecTRI / DoRothEA (TF regulons)

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
├── r_analiz_melanoma/
│   ├── 01_preprocessing.R
│   ├── 02_DEG_NOISeq.R
│   ├── 03_enrichment.R
│   ├── 04_visualizations.R
│   ├── 05_mirna_mrna_integration.R
│   ├── 06_TF_enrichment.R
│   ├── 07_WGCNA.R
│   ├── 08_legacy_wgcna_tf_linkage.R
│   └── main pipeline 2.R
├── data/
│   └── README.md                    # GEO download instructions (raw data not included)
├── output/ or data/ilk/newoutput/    # analysis outputs, depending on the local setup
└── README.md
```

---

## Planned Extensions

### Legacy TF/WGCNA linkage
The old WGCNA `2.13` comparison logic was split out of the core module stage and moved into `08_legacy_wgcna_tf_linkage.R`. It is optional and exists mainly for compatibility with the earlier interpretation chain that tied TF leading-edge hits to a small legacy target-gene list.

### PPI / TF–mRNA Network *(in progress)*
Integration of STRING protein interaction data with WGCNA hub genes and TF candidate list to identify high-connectivity resistance hubs. Outputs include annotated hub table and Cytoscape-compatible network files.

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
  "multiMiR", "decoupleR", "org.Hs.eg.db",
  "ReactomePA", "msigdbr", "enrichplot"
))
```

Visualization: `ggplot2`, `pheatmap`, `Cairo` (Windows Unicode fix for PDF output)  
Network export: Cytoscape (CX-format files)

---

## Data Availability

**Raw counts:** GEO [GSE283251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283251) — download instructions in `data/README.md`  
**Analysis outputs** (DEG tables, pathway results, network edges, WGCNA module files): [Zenodo — *DOI pending after deposit*]

---

## Reference

Çolako\u011flu Bergel N, et al. (2025). *[Full title]*. *Scientific Reports*. DOI: *[add when available]*

---

## Notes on Reproduction

- NOISeq computes M = log2(A/B), so for the SC vs RC comparison, **positive M = higher in SC**. The stage scripts flip signs where needed so downstream tables consistently interpret positive log2FC as higher in the second comparison group.
- Histone gene filter applied (replication-dependent HIST1/HIST2 clusters removed before enrichment).
- BH correction is not applied to NOISeq probability scores; `prob ≥ 0.8` threshold used directly, consistent with NOISeq documentation.
- `fgsea` ties warning (~85% tied ranks) is expected from NOISeq's integer count ratios and does not indicate a pipeline error; S10 vs R10 GSEA results are the least reliable comparison for this reason.
- The checkpoint chain is explicit: `01 -> 02 -> 03 -> 04`, `03 -> 05`, `03 -> 06`, `03 -> 07`, `06 + 07 -> 08`.
