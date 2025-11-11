# HZAU Cohort WGS Data Processing Pipeline

## Introduction

This is the WGS data processing pipeline for the HZAU cohort study, which includes the following steps:

1. Sequencing file quality control (completed by sequencing company)
2. Reference genome alignment
   - Align to GRCh38 reference genome using pbmm2
   - Merge all Movie sequencing files for each sample
   - Sort and index using samtools
   - Perform depth and coverage statistics
3. Variant calling
   - SNP
      1. Call short variants (SNP, INDEL) using pepper-deepvariant
      2. Merge and filter population variants using glnexus
      3. Perform SNP annotation
   - SV
      1. Call structural variants using three software tools (to be determined)
      2. Merge results from three tools at individual sample level - sample variant set
      3. Merge identification results at population level - cohort variant set
4. Methylation calling (5mc)
   - Call methylation using jasmine (completed by sequencing company)
   - Calculate methylation probabilities using pb-cpg-tools - sample methylation results
   - Merge methylation probabilities from all samples - cohort methylation results

## Directory Structure

```plaintext
WGS
├── raw_data:        Not synchronized   Stores the raw sequencing data, this data should not be modified
├── benchmarks:      Not synchronized   Stores the resource statistics for each step, created automatically during runtime
├── configs:         Synchronized       Workflows configurations
├── profiles:        Synchronized       SnakeMake runtime parameter configurations
├── logs:            Not synchronized   Stores the output during runtime, created automatically during runtime
├── resources:       Not synchronized   Input data and reference data, recommended to use symbolic links
│   ├── fastq:       Not synchronized   Stores fastq format sequencing files, in general, your input is of this format
│   ├── reference:   Not synchronized   Reference files, such as reference genome
│   └── uBAM:        Not synchronized   Stores bam format sequencing files, in general, your input is of this format
├── results:         Not synchronized   Output results of the workflow
└── workflow:        Synchronized       All the workflows
   ├── envs:         Synchronized       Conda environment configurations required for each workflow
   ├── rules:        Synchronized       Workflows
   ├── scripts:      Synchronized       Some auxiliary scripts
   └── Snakefile:    Synchronized       Entry point of the workflow
```

## Usage

### Install Appropriate SnakeMake Version

This pipeline is written using SnakeMake v8.2.3. Since SnakeMake v7 to v8 introduced breaking changes, the SnakeMake version used to run this pipeline should not be lower than 8.

You can follow the official tutorial to install using conda:

```bash
conda create -n snakemake bioconda::snakemake
```

Then, you need to clone this repository to your local machine and enter the directory:

```bash
git clone https://gitee.com/eastsunw/pacbio-wgs-snakemake-pipeline.git
cd pacbio-wgs-snakemake-pipeline
```

### Prepare Data

Refer to the directory description below to prepare your data. Note that directories marked as `Not synchronized` will be automatically generated during runtime or need to be provided by you. You need to place sequencing files of the same format in one folder and organize the filenames to match the regular expression in the configuration. Your expression should contain at least the following elements:

```plaintext
<sample_id>__<movie_id>.bam
```

For example: `wholeblood_16A__m84128_231117_083350_s3_hifi_reads_bc2082.bam`. The movie_id can contain various information separated by other delimiters, but do not use `__`.

Next, you need to place your reference genome and other information in `resources/reference`. The directory structure here is not critical as you can specify the paths of required files in the configuration file.

### Modify Configuration Files

After preparing the data, you should go to the configs directory to modify your configuration files. The configuration files are in YAML format. There are two files in this directory: `config.yaml` and `config.yaml`. The former is used to configure your cohort information, and the latter contains software parameters for each workflow.

### Run the Pipeline

After all preparations are completed, run the following command in the directory to start the pipeline:

```bash
conda activate snakemake
snakemake --profile workflow/profiles
```

If you only want to test whether the pipeline works normally, you can use the `-n` or `--dry-run` parameter for dry-run:

```bash
conda activate snakemake
snakemake --profile workflow/profiles --dry-run
```

## Notes

1. This pipeline is only suitable for PacBio WGS sequencing data and not for other types of data.
2. Some environments are local and need to be manually configured by the user. These environments are typically required because the Python version and the software cannot coexist, while other environments are specified through yaml files and do not require manual configuration.

        