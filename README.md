# Neutrophil Orchestrates Pre-Metastatic Niche Formation and Immuno-suppressive Microenvironment Development for Prostate Cancer Bone Metastasis

This repository contains the source code and analysis pipeline for the single-cell RNA sequencing (scRNA-seq) data analysis presented in our manuscript.

## 📄 Project Overview

Bone metastasis represents a lethal stage of prostate cancer, yet the dynamic reprogramming of the bone marrow (BM) microenvironment remains poorly defined. In this study, we utilized single-cell profiling of BM samples across multiple disease stages to uncover a critical signaling axis orchestrated by neutrophils.

**Abstract:**

The dynamic reprogramming of the bone marrow (BM) microenvironment during metastatic progression is poorly defined. Through single-cell profiling of BM samples across multiple disease stages, we uncovered a primary tumor-to-BM signaling axis that orchestrates niche formation through neutrophil-mediated remodeling. Tumor-associated neutrophils (TANs) emerge during the pre-metastatic stage, persist throughout metastatic progression, correlate with poor patient prognosis, and are sufficient to establish an immunosuppressive niche that facilitates metastatic colonization.

Mechanistically, we demonstrate that a specific **tumor-derived growth factor (Factor X)** drives the differentiation of TANs. These TANs disrupt immune surveillance within the pre-metastatic niche (PMN) by suppressing CD8+ T cell effector functions and blocking the differentiation of anti-metastatic monocytes. Within established metastases, TANs induce immunosuppression by inducing CD8+ T cell hypoxia and promoting regulatory T cell (Treg) expansion.

Pharmacological inhibition of the **Factor X Receptor (Target R)** abrogated PMN formation and potentiated anti-PD-1 therapy, achieving enhanced suppression of bone metastasis. Translationally, our investigator-initiated clinical trial revealed that combined **Target R** and PD-1 blockade elicits clinically meaningful responses in metastatic castration-resistant prostate cancer (mCRPC) patients. Collectively, these findings delineate a pathogenic signaling cascade from primary **tumor-derived Factor X** to systemic immunosuppression via TANs, positioning this axis as an actionable immunotherapy target.

> **Note:** Specific gene names and therapeutic targets have been anonymized in this repository description pending peer review and publication.

## 📂 Repository Structure

The analysis is organized by cell type modules to ensure reproducibility and clarity. All scripts share a common setup file for environment consistency.

```text
PCa_BM/
├── analysis/                       # Core analytical pipelines and processing scripts
│   ├── atac_seq_analysis.sh        # Pipeline for ATAC-seq data analysis
│   ├── clustering.r                # Dimensionality reduction and initial clustering
│   ├── compass.r                   # Metabolic flux analysis using COMPASS
│   ├── infercnv.r                  # Copy Number Variation (CNV) inference pipeline
│   ├── preprocess.r                # Quality control, filtering, and normalization
│   ├── pySCENIC.script.py          # Python implementation for SCENIC regulon analysis
│   └── scenic.sh                   # Shell script wrapper for running SCENIC
│
├── Plotting/                       # Downstream analysis and figure generation by cell type
│   ├── 00_Setup_and_Functions.R    # Global environment setup and custom helper functions
│   ├── 01_Myeloid_Analysis.R       # Neutrophils, Monocytes, and Macrophages analysis
│   ├── 02_T_NK_Analysis.R          # CD8+, CD4+ T cells, and NK cells analysis
│   ├── 03_Stromal_Analysis.R       # Osteoblasts, Endothelial cells, and Fibroblasts analysis
│   ├── 04_Epithelial_Analysis.R    # Tumor cell identification and characterization
│   ├── 05_B_Cell_Analysis.R        # B cell lineage and tissue preference analysis
│   ├── Other_tcells.r              # Supplementary analysis for specific T cell subsets
│   └── health_BM.r                 # Comparative analysis with healthy bone marrow reference
│
└── data/                           # Input data directory (Raw counts, metadata, etc.)
```

## 🛠 System Requirements & Dependencies

The code was developed and tested using **R (v4.x)**. Key dependencies include:

  * **Seurat** (v4/v5) - Single-cell data analysis
  * **Monocle** (v2) / **scVelo** - Trajectory inference
  * **NicheNet** / **CellPhoneDB** - Cell-cell communication
  * **SCENIC** (pySCENIC) - Regulon analysis
  * **inferCNV** - Copy number variation inference
  * **COMPASS** - Metabolic flux analysis

Please refer to `00_Setup_and_Functions.R` for the complete list of libraries.

## 🚀 Usage Instructions

1.  **Setup Environment:** Run `00_Setup_and_Functions.R` first to load necessary libraries and define custom functions (e.g., `quick_fastmnn`, `dotp`).
2.  **Data Placement:** Ensure raw `.h5seurat` or count matrices are placed in the `data/raw/` directory as indicated in the scripts.
3.  **Step-by-Step Analysis:** Run the numbered scripts (01-05) sequentially or independently based on the cell type of interest.

## 📞 Contact

For questions regarding the code or data availability, please contact the corresponding author or open an issue in this repository.

-----

*Last Updated: Nov 2025*
