# Data Statistics

## Sample Distribution Statistics

- **Statistical Methods**: Age (under 40, 40-60, over 60), Gender, Region (map visualization)
- **Input Data**: `sudo ln -s /home/data_admin/Clinic/Sample_info/sample_info.txt data/sample_info.txt`
- **Plotting Script**: `scripts/plot/plot_sample_infos.R`
    - Sample region and gender distribution: `results/plot/cohort_sample_region.pdf`
    - Sample age distribution: `results/plot/cohort_sample_age.pdf`

## Sequencing Alignment Statistics

- **Statistical Methods**: Extract information from sequencing and alignment results using Python scripts, including sequencing depth, reference genome alignment coverage, and alignment rate
- **Input Data**:
    ```bash
    ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/alignment/GRCh38_241031.depth.summary.txt data/depth.summary.txt
    ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/alignment/GRCh38_241031.mappingrate.summary.txt data/mappingrate.summary.txt
    ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/alignment/GRCh38_241031.coverage.summary.txt data/coverage.summary.txt
    ```
- **Statistical Results**:
    - Sequencing depth: `data/depth.summary.txt`
    - Alignment rate: `data/mappingrate.summary.txt`
    - Alignment coverage: `data/coverage.summary.txt`
- **Plotting Script**: `scripts/plot/plot_sequencing_stat.R`
- **Plot Results**:
    - Sequencing depth: `results/plot/sampleSeq_sequencing_depth.pdf`
    - Alignment rate: `results/plot/sampleSeq_alignment_rate.pdf`
    - Alignment coverage: `results/plot/sampleSeq_alignment_coverage.pdf`

## Alignment Error Rate Statistics

- **Statistical Results**: `cp /home/wtian/project/HZAU_cohort_meth/wdy_assist/01_sequencing_error_rate/results/all_samples_stats.tsv data/error_stat.txt`
- **Plotting Script**: `scripts/plot/plot_error_stat.R`
- **Plot Results**: `results/plot/sampleSeq_mapping_error.pdf`

## SNP Statistics

- **Statistical Methods**: Generate SNP and InDel BED files from cohort SNP results using Python scripts, containing chr, ID, start, end, AF (MAF)
- **Input Data**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SNP/clair3_cohort/GRCh38_241031.cohort.vcf.gz data/cohort_SNP_InDel.vcf.gz`
- **Statistical Script**: `scripts/stat/stat_SNP_InDel_2Bed.py`
- **Statistical Results**:
    - SNP_all_bed: `results/bed/SNP.all.bed`
    - SNP_common_bed: `results/bed/SNP.common.bed`
    - SNP_rare_bed: `results/bed/SNP.rare.bed`
    - InDel_all_bed: `results/bed/InDel.all.bed`
    - InDel_common_bed: `results/bed/InDel.common.bed`
    - InDel_rare_bed: `results/bed/InDel.rare.bed`

## MNV Statistics

- **Statistical Methods**: Convert MNV raw data to BED format using Python scripts, containing chr, ID, start, end, AF (MAF), type
- **Input Data**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/MNV/clair3/GRCh38_241031.mnv.withID.txt data/cohort_MNV.txt`
- **Statistical Script**: `scripts/stat/stat_MNV2Bed.py`
- **Statistical Results**:
    - MNV_all_bed: `results/bed/MNV.all.bed`
    - MNV_common_bed: `results/bed/MNV.common.bed`
    - MNV_rare_bed: `results/bed/MNV.rare.bed`
- **Plotting Script**: `scripts/plot/plot_cohortMNV_stat.R`
- **Plot Results**: Different MNV frequency counts: `results/plot/cohort_MNV_freq_count.pdf`

## SV Statistics

### Software Overlap in SV Detection

- **Statistical Methods**: Obtain overlapping SV counts from each software in merged SV datasets
- **Input Data**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SV/sample_merged data/sample_SV_merged`
- **Statistical Script**: `scripts/stat/stat_sampleSV_software_overlap.py`
- **Statistical Results**: `results/stat/sampleSV_software_overlap.txt`
- **Plotting Script**: `scripts/plot/plot_sampleSV_software_overlap.R`
- **Plot Results**:
    - Euler diagram (area reflects count): `results/plot/sampleSV_software_overlap_euler.pdf`
    - Venn diagram (pure count): `results/plot/sampleSV_software_overlap_venn.pdf`

### Non-redundant SV Count Statistics

- **Statistical Methods**: Count SV types in non-redundant SV data for each sample
- **Input Data**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SV/sample_filtered data/sample_SV_filtered`
- **Statistical Script**: `scripts/stat/stat_sampleSV_type_count.py`
- **Statistical Results**: `results/stat/Variant_stat/sampleSV_type_count.txt`
- **Plotting Script**: `scripts/plot/plot_sampleSV_type_count.R`
- **Plot Results**:
    - SV counts by frequency: `results/plot/sampleSV_type_count.pdf`
    - Non-redundant SV statistics table: `results/plot/sampleSV_type_stat.html` (requires manual conversion to PDF)

### Tandem Repeats to BED Conversion

- **Input Data**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/TR/cohort_trgt/GRCh38_241031.normal.vcf.gz data/cohort_TR.vcf.gz`
- **Statistical Script**: `scripts/stat/stat_TR2Bed.py`
- **Statistical Results**: `results/bed/TR.all.bed`

### SV Filtering Step Remaining Counts

- **Statistical Methods**: Count remaining SVs after different filtering steps for each sample
- **Input Data**: `data/sample_SV_filtered`
- **Statistical Script**: `scripts/stat/stat_sampleSV_filter_remain.py`
- **Statistical Results**: `results/stat/Variant_stat/sampleSV_filter_remain.txt`
- **Plotting Script**: `scripts/plot/plot_sampleSV_filter_remain.R`
- **Plot Results**: `results/plot/sampleSV_filter_remain.pdf`

### Cohort Non-redundant SV Feature Statistics

- **Statistical Methods**: Convert cohort SV data to BED format using Python scripts, containing chr, ID, start, end, AF (MAF), type, len, support; visualize SV counts, length distribution, support numbers
- **Input Data**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SV/cohort_merged/GRCh38_241031.SV.vcf.gz data/cohort_SV.vcf.gz`
- **Statistical Script**: `scripts/stat/stat_SV2Bed.py`
- **Statistical Results**:
    - SV_all_bed: `results/bed/SV.all.bed`
    - SV_common_bed: `results/bed/SV.common.bed`
    - SV_rare_bed: `results/bed/SV.rare.bed`
- **Plotting Script**: `scripts/plot/plot_cohortSV_stat.R`
- **Plot Results**:
    - SV length distribution by type: `results/plot/cohortSV_length.pdf`
    - Mean support samples by SV type: `results/plot/cohortSV_meanSupport.pdf`

### Genome-wide Distribution of Various SVs

- **Input Data**:
    - Gene: `results/stat/circos_density/gene_density.txt`
    - Repeat: `results/stat/circos_density/repeat_density.txt`
    - SV-All: `results/stat/circos_density/SV_all_density.txt`
    - SV-INS: `results/stat/circos_density/SV_INS_density.txt`
    - SV-DEL: `results/stat/circos_density/SV_DEL_density.txt`
    - SV-INV: `results/stat/circos_density/SV_INV_density.txt`
    - SV-DUP: `results/stat/circos_density/SV_DUP_density.txt`
- **Plotting Script**: `scripts/plot/plot_cohortSV_circos.R`
- **Plot Results**: Circos plot of different SV types: `results/plot/cohortSV_circos.pdf`

## All Variant Type Statistics

### Variant Counts by Frequency

- **Input Data**:
    - SNP-Common: `results/bed/SNP.common.bed`
    - SNP-Rare: `results/bed/SNP.rare.bed`
    - InDel-Common: `results/bed/InDel.common.bed`
    - InDel-Rare: `results/bed/InDel.rare.bed`
    - MNV-Common: `results/bed/MNV.common.bed`
    - MNV-Rare: `results/bed/MNV.rare.bed`
    - SV-Common: `results/bed/SV.common.bed`
    - SV-Rare: `results/bed/SV.rare.bed`
- **Statistical Script**: `scripts/stat/stat_cohortVariant_amount_byFreq.py`
- **Statistical Results**: `results/stat/cohortVariant_amount_byFreq.txt`
- **Plotting Script**: `scripts/plot/plot_cohortVariant_amount_byFreq.R`
- **Plot Results**: `results/plot/cohortVariant_amount_byFreq.pdf`

### Variant Accumulation Curve with Sample Increase

- **Statistical Methods**: For SNPs, count unique sites with increasing samples; for SVs, merge cohort SVs progressively and count identified SVs at different sample sizes
- **Statistical Results**:
    - SNP: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SNP/SNP_sample_increase.stat.txt data/cohort_SNP_sample_increase.txt`
    - SV: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SV/cohort_step/GRCh38_241031_sample_increase_stat.txt data/cohort_SV_sample_increase.txt`
- **Plotting Script**: `scripts/plot/plot_cohortVariant_sample_increase.R`
- **Plot Results**: `results/plot/cohortVariant_sample_increase.pdf`

### Genome-wide Coverage of All Variant Types

- **Statistical Methods**: Use `bedtools merge -i <variant>.bed > <variant>.merged.bed` to merge variant BED files, import to R for coverage statistics
- **Input Data**: `scripts/stat/merge_variant_bed.sh`
    - Cohort SNP: all/common/rare merged BED files
    - Cohort InDel: all/common/rare merged BED files
    - Cohort MNV: all/common/rare merged BED files
    - Cohort SV: all/common/rare merged BED files
- **Statistical Script**: `scripts/stat/stat_cohortVariant_coverage.R`
- **Statistical Results**: `results/stat/cohortVariant_coverage.txt`

### Genome-wide Distribution of All Variant Types (Circos)

- **Input Data**:
    - Gene positions: Converted from GTF to BED
    - Gene, Repeat, SNP, InDel, MNV, SV, TR density files
- **Statistical Script**: `scripts/stat/calculate_density.R`
- **Statistical Results**: Density files for Gene, Repeat, SNP, InDel, MNV, SV, TR
- **Plotting Script**: `scripts/plot/plot_cohortVariant_circos.R`
- **Plot Results**: `results/plot/cohortVariant_circos.pdf`

## GWAS Association Results Statistics

### Merge and Filter GWAS Results by Variant Type

- **Statistical Methods**: Add MARKER column to each variant's original list, merge results, filter variants using 5e-8 threshold
- **Input Data**: GWAS results for SNP, InDel, MNV, SV
- **Statistical Script**: `scripts/stat/stat_merge_GWAS.py`
- **Statistical Results**: Filtered GWAS results for each variant type

## QTL Association Results Statistics

### Filter and Count QTL Results

- **Statistical Methods**: Filter using 5e-8 threshold
- **Input Data**: QTL results for various variant-phenotype combinations
- **Statistical Script**: `scripts/stat/stat_QTL_merge_and_filter_QTL.py`
- **Statistical Results**: Filtered QTL results for each variant-phenotype combination

### Visualize QTL Results

- **Input Data**: QTL results and variant genotype-phenotype data
- **Plotting Script**: `scripts/plot/plot_QTL_result.R` with parameters for QTL type, variant ID, phenotype ID
- **Plot Results**: QTL plots showing phenotype distribution by genotype

## Clinical Indicator Range Results Statistics

### Data Preparation

1. Indicator ranges: `ln -s /home/data_admin/Clinic/Biochemistry.all.level.txt data/Clinic/bio_level.txt`
2. Indicator risks: `ln -s /home/data_admin/Clinic/Biochemistry.all.risk.txt data/Clinic/bio_risk.txt`
3. Indicator information: `ln -s /home/data_admin/Clinic/Biochemistry/Biochemistry.all.txt data/Clinic/bio_info.txt`
