# Directory Contents

This directory contains scripts for Linkage Disequilibrium (LD) pattern analysis between complex variants (CVs) and single nucleotide variants (SNVs).

## Workflow Overview

The analysis pipeline consists of five main steps:

### 1. PLINK to VCF Conversion (`01_plink2vcf.sh`)
- **Purpose**: Convert PLINK binary format files (.bed/.bim/.fam) to VCF format

### 2. VCF Formatting (`02_vcf_formatter.py`)
- **Purpose**: Standardize variant IDs and genotypes across different variant types

### 3. VCF File Merging (`03_all_vartype_vcf_merge.sh`)
- **Purpose**: Merge multiple VCF files (SNP, InDel, MNV, SV) into a single VCF

### 4. LD Analysis (`04_ld_snakemake.smk`)
- **Purpose**: Perform linkage disequilibrium analysis using PLINK

### 5. SV QTL Tag Rate Calculation (`05_calc_SV_QTL_tag_rate.R`)
- **Purpose**: Calculate SV QTL tag rate and analyze SNP QTL overlap

## Data Flow

```
PLINK files → VCF conversion → Formatting → Merging → LD analysis → Tag rate calculation
```

## Dependencies

- **Tools**: PLINK, bcftools, bgzip, tabix, Snakemake
- **Languages**: Bash, Python, R
- **R Packages**: tidyverse, data.table, stringr

## Output Structure

- `01_flt_vcf/`: Converted VCF files
- `02_vcf_format_simplified/`: Formatted VCF files  
- `03_all_vartype_merged.vcf`: Merged VCF file
- `04_ld_result/`: Chromosome-specific LD results
- Statistical outputs from tag rate analysis

        