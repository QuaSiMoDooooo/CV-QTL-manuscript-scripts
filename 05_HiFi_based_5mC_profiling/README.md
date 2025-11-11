# Directory Contents

This directory contains scripts for HiFi-based 5mC profiling analysis, including methylation site annotation, probe mapping, and comparison with external databases.

## Workflow Overview

The analysis pipeline consists of three main components:

### 1. HiFi Methylation Annotation (`hifi_meth_annot/`)
- **Purpose**: Annotate HiFi-based 5mC sites with genomic context and map probe IDs
- **Scripts**:
  - `annotatr.R`: Annotate methylation sites with genomic features (CpG islands, promoters, exons, etc.) using annotatr package
  - `map_probe_id.R`: Map HM450K and EPIC array probe IDs to hg38 coordinates using easylift package

### 2. Methylation Levels Comparison (`levels_comparison/`)
- **Purpose**: Compare methylation levels between HZAU cohort and external methylation databases
- **Scripts**:
  - `methbank_CRMs.R`: Compare with MethBank CRMs database (whole blood 850K data), calculate average methylation levels and classify into high/medium/low categories
  - `GSE186458_blood_WGBS.R`: Compare with GSE186458 blood WGBS data, load and merge samples to calculate average methylation levels

### 3. Methylation Sites Comparison (`sites_comparison/`)
- **Purpose**: Analyze overlap between different methylation detection platforms
- **Scripts**:
  - `450k_850k_WGBS.R`: Compare site overlaps between HiFi-based 5mC, WGBS, EPIC array, and HM450 array