# Mediation Analysis

## Directory Contents

### Analysis Pipeline

#### 01_cis_QTL_overlap/
- **SNV_QTL.R**: Process cis-QTL data for SNP variants and filter variants with multiple molecular phenotypes
- **SV_QTL.R**: Process cis-QTL data for SV variants and filter variants with multiple molecular phenotypes

#### 01_QTL_GWAS_overlap/
- **code.R**: Identify overlapping variants between QTL and GWAS data

#### 02_format/
- **01_cis_QTL_overlap.R**: Format cis-QTL overlap data for mediation analysis
- **02_QTL_GWAS_overlap.R**: Format QTL and GWAS overlap data for mediation analysis

#### 03_batch_mediation_analysis/
- **02_cis_QTL_overlap_batch.R**: Batch mediation analysis for cis-QTL overlap data using medflex package
- **02_QTL_GWAS_overlap_batch.R**: Batch mediation analysis for QTL and GWAS overlap data using medflex package

#### 04_batch_clean/
- **01_summarize_med_res.R**: Summarize mediation analysis results from batch processing
- **02_med_res_all_stats.R**: Statistical analysis of mediation results including effect significance and mediation rates
- **03_med_rates_with_same_molGene_under_sig_TE_ME.R**: Calculate mediation rates for same molecular gene under significant total and mediated effects
- **04_med_rates_coloc_under_sig_TE_ME.R**: Calculate mediation rates under significant total and mediated effects of colocalization results
- **05_med_TF.R**: Transcription factor enrichment analysis for mediation results
