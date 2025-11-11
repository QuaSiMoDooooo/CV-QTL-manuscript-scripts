#!/usr/bin/env Rscript

# -------------
# FileName     : 07_smr_all_cis_flt_stats.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Statistical analysis of filtered cis SMR results
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR")
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

data = fread("05_smr_all_cis_flt.tsv")

length(unique(data$probeID))
length(unique(data$trait))

data %>%
  group_by(trait, topSNP) %>%
  summarise(unique_mol_phe_count = n_distinct(mol_phe)) %>%
  arrange(desc(unique_mol_phe_count))

length(unique(data$topSNP))

data %>%
  group_by(topSNP) %>%
  summarise(unique_trait_count = n_distinct(trait)) %>%
  filter(unique_trait_count > 1) %>%
  nrow()

data[data$topSNP=="rs3828807",]