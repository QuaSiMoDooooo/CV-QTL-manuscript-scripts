# Code repository for the CV-QTL manuscript

This repository contains all the analysis scripts for the CV-QTL manuscript. The scripts are organized by sub-folder following the order presented in the manuscript.

## Project Structure

| Folder Number | Directory Name | Description | Key Components |
|----------------|----------------|-------------|----------------|
| 01 | `01_WGS/` | Whole Genome Sequencing data processing pipeline | Quality control, alignment, variant calling (SNP/SV), methylation calling |
| 02 | `02_RNAseq/` | RNA-seq data analysis pipeline | Alternative polyadenylation, splicing, gene expression quantification, alignment |
| 03 | `03_SeqDownstreams_and_Association/` | Downstream analysis of sequencing data | RNAseq quality check, data statistics, variant annotation, QTL/GWAS association |
| 04 | `04_CV_LD_pattern/` | Linkage Disequilibrium pattern analysis | PLINK to VCF conversion, LD analysis, SV QTL tag rate calculation |
| 05 | `05_HiFi_based_5mC_profiling/` | HiFi-based 5mC profiling analysis | Methylation annotation, probe mapping, comparison with external databases |
| 06 | `06_CV_QTL_characteristic/` | Complex Variant QTL characteristics analysis | SNP QTL comparison, cis-eQTL case analysis, trans-QTL hotspot identification |
| 07 | `07_clinical_measurements_GWAS/` | Clinical measurements GWAS analysis | Conditional analysis, fine mapping, Manhattan plot generation |
| 08 | `08_GWASCatalog_SV_signal/` | GWAS Catalog SV signal analysis | LD-based tagging, visualization, Enformer model utilities |
| 09 | `09_coloc_analysis/` | Colocalization analysis | GWAS-QTL colocalization, statistical analysis, network visualization |
| 10 | `10_mediation_analysis/` | Mediation analysis | QTL-GWAS overlap, batch mediation analysis, transcription factor enrichment |
| 11 | `11_SMR_analysis/` | Summary-data-based Mendelian Randomization analysis | BESD file generation, SMR analysis, result filtering and visualization |
