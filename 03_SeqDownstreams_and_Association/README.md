# Directory Contents

Upstream analysis includes all sequencing data processing and data generation directly identified from sequencing data, such as sequence alignment, variant identification, expression quantification, etc.

Downstream analysis includes using upstream analysis data for data mining to obtain biological conclusions, such as QTL association, GWAS, etc.

### RNAseq Quality Check

- **Background**: The sequencing company reported RNA degradation in RNAseq samples and classified samples into three grades (A, B, C) based on degradation status
- **Objective**: Check whether gene expression in A/B/C samples is affected by RNA degradation, and further examine factors affecting gene expression
- **Methods**:
  - Literature review indicates RNA degradation manifests as 3' end shortening. First check alignment of sequencing fragments to transcripts to examine 3' end shortening
  - Use different variables for correction, compare PCA results before and after correction to determine which factor(s) have the greatest impact on gene expression. Three variables used for correction with 1, 2, and 3 variables respectively:
    1. Sample quality grouping by sequencing company (A, B, C)
    2. Clustering based on transcriptome alignment results (1, 2, 3)
    3. Different sequencing result delivery batches (0114, 0429, 0531)
- **Results**:
  - From transcript alignment results, sample grading and RNA 3' shortening show no clear correlation. Each group has three distinct subgroups. Housekeeping gene expression tests based on clustering results show significant differences, while tests based on sample grouping show no significant differences
  - With single variable correction, both batch and quality grouping can eliminate inter-group differences
  - With two variables, using clustering results and delivery batch works well but is weaker than single variable
  - With three variables, the effect is similar to two variables

Overall, using clustering results for correction works best

### data_stat

This section includes overview statistics for all data, including:
- Sequencing data statistics and visualization
- Variant count statistics and visualization
- Phenotype data statistics and visualization
- Statistics and visualization comparisons with other datasets

### Variant Annotation
This section includes annotation of all variants, including:
- Regional enrichment of variants
- Network structure data for database construction
- Variant comparison
- pLOF statistics

### HiFi-based 5mC Annotation

HiFi sequencing 5mC annotation, including:
- 5mC chip ID conversion
- 5mC genomic region annotation
- Regional methylation level statistics

### phenotype_qtl

Association analysis between phenotypes and variant genotypes. Sequencing molecular phenotypes use QTL association methods, while physical examination biochemical indicators use GWAS methods.

Genotype quality control:
- SNP: Missing rate <5%, MAF >5%, HWE >1e-6, non-multiallelic variants
- MNV: Missing rate <5%, MAF >5%
- SV: Missing rate <5%, MAF >5%

Phenotype quality control:
- Gene expression: Average expression >0.1 (log2(TPM+1))
- APA: Missing rate <5%
- AS: Missing rate <5%
- Physical examination indicators: Only numerical values are associated, numerical values distinguish gender-independent and sex hormones

Covariates:
- QTL: Age, gender, 15 PEER, 5 PCA
- GWAS: Age, gender (non-sex hormones), 5 PCA

Association:
- QTL: Exclude Y and M chromosomes, cis window 1Mb
- GWAS: Exclude Y and M chromosomes