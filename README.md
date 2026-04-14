# Ketamine-Spatial-Analysis
Integrated R & Python pipeline for spatial reconstruction of prefrontal circuits.
## 1. Overview
This repository contains the spatial transcriptomics (ST) workflow used to validate the spatial organization and crosstalk potential of 8 key cell subtypes in the PFC.
## 2. Script Details
### A. Preprocessing & Integration (process.r)
Function: Standardizes merFISH data and manages cross-species ortholog mapping.  
Key Logic: Utilizes convertHumanToMouse (Tax 9606 to 10090) with a "sum" expression strategy.  
Interoperability: Converts Seurat objects to .h5ad format via srt_to_adata for Tangram compatibility.
### B. Spatial Mapping (run_tangram2.py)
Algorithm: Tangram (v1.0.4).
Stage 1: Maps snRNA-seq nuclei to spatial coordinates using 'cells' mode (Supporting Figure 6a-h).
Stage 2: Calculates microenvironment compositions and generates spatial pie charts (Supporting Figure 6i).
### 3. Requirements
R: Seurat (v3.1.4), homologene.  
Python: tangram-sc (v1.0.4).
