# Colocalization Analysis

This directory contains scripts for performing colocalization analysis between GWAS signals and molecular QTLs (expression, splicing, APA, methylation) to identify shared genetic associations.

## Directory Contents

### Analysis Pipeline

#### 01_get_chr1-22_GWAS_and_eQTL_sentinelVARs_range/
- **code.R**: Extract GWAS and eQTL sentinel variants within 1Mb ranges for colocalization analysis

#### 02_get_molphe_within_sentinelVARs_range/
- **01_merge_bed_files.sh**: Merge multiple BED files and extract genomic regions
- **02_get_molphe.R**: Filter molecular phenotypes within sentinel variant ranges

#### 03_fliter_molphe_matched_data/
- **code.R**: Filter molecular phenotype matched data for colocalization analysis

#### 04_rerun_QTL/
- **process_QTL.smk**: Snakemake workflow for QTL analysis with parallel processing
- **split_QTL_molphe.R**: Split QTL results by molecular phenotype for colocalization analysis
- **scripts/tools/**: Utility scripts for QTL analysis

#### 05_coloc/
- **01_01_get_overlap_phe.R**: Identify molecular phenotypes overlapping with GWAS sentinel variant ranges
- **01_02_gwas_coloc_all_qtl_each_file.R**: Perform colocalization analysis between GWAS and QTL data
- **02_01_get_overlap_phe.R**: Identify molecular phenotypes overlapping with QTL sentinel variant ranges
- **02_02_eqtl_qtl_each_file.R**: Perform colocalization analysis between eQTL and other QTL data
- **03_01_clean_res.R**: Clean and filter colocalization results
- **03_02_stats.R**: Statistical analysis of colocalization results
- **04_trait_or_gene_based_graph_vis.R**: Visualize colocalization results as network graphs