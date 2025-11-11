#!/usr/bin/env Rscript

# -------------
# FileName     : 450k_850k_WGBS.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Compare methylation site overlaps between HiFi-based 5mC, WGBS, EPIC array, and HM450 array
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(pryr))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

hzau <- fread("01_all_3kw_meth_sites_uniq.tsv")
hzau_sites <- paste0(hzau$chr, "_", hzau$end)

WGBS <- fread("Blood_samples_average_methylation_sites.bed")
WGBS_sites <- paste0(WGBS$chr, "_", WGBS$start+1)

array850 <- fread("EPIC.manifest.hg38.tsv")
array850_sites <- paste0(array850$chr, "_", array850$hg38_start)

array450 <- fread("hm450.manifest.hg38.tsv")
array450_sites <- paste0(array450$chr, "_", array450$hg38_start)

print("Overlap between HZAU and WGBS:")
table(WGBS_sites %in% hzau_sites)

print("Overlap between HZAU and EPIC array:")
table(array850_sites %in% hzau_sites)

print("Overlap between HZAU and HM450 array:")
table(array450_sites %in% hzau_sites)