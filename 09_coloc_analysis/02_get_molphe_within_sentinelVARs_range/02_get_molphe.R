#!/usr/bin/env Rscript

# -------------
# FileName     : 02_get_molphe.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Filter molecular phenotypes within sentinel variant ranges for colocalization analysis
# -------------

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

options(stringsAsFactors = FALSE)
options(scipen = 999)
mc <- getOption('mc.cores', 36)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/02_get_molphe_within_sentinelVARs_range")

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/data/QTL_flt"
range_dir = "."
qtl_paths = list.files(qtl_dir, pattern = "*.tsv$", full.names = TRUE)
range_paths = list.files(range_dir, pattern = "*.bed$", full.names = TRUE)
qtl_var_type_list = str_split(str_split(basename(qtl_paths), "\\.", simplify = TRUE)[, 1], "_", simplify = TRUE)[, 3]
range_var_type_list = str_split(str_split(basename(range_paths), "\\.", simplify = TRUE)[, 1], "_", simplify = TRUE)[, 3]
var_type_list = intersect(qtl_var_type_list, range_var_type_list)

filter_molphe_by_range <- function(qtl_df_dist, range) {
    range_gr <- GRanges(
        seqnames = range$chr,
        ranges = IRanges(start = range$start, end = range$end)
    )
    qtl_gr <- GRanges(
        seqnames = qtl_df_dist$chr,
        ranges = IRanges(start = qtl_df_dist$start, end = qtl_df_dist$end)
    )
    hits <- findOverlaps(qtl_gr, range_gr)
    qtl_df_dist_overlap <- qtl_df_dist[queryHits(hits), ]
    return(qtl_df_dist_overlap)
}

for (j in 1:length(var_type_list)) {
    var_type = var_type_list[j]
    qtl_path = qtl_paths[qtl_var_type_list == var_type]
    range_path = range_paths[range_var_type_list == var_type]
    qtl = fread(qtl_path)
    
    tmp = str_split(str_split(qtl$Phenotype_position, ":", simplify = TRUE)[,2],"-", simplify = TRUE)
    start = tmp[,1]
    end = tmp[,2]
    qtl_df = data.frame(
        molphe = qtl$molphe,
        Phenotype = qtl$Phenotype,
        chr = str_split(qtl$Phenotype_position, ":", simplify = TRUE)[,1],
        start = as.numeric(start), 
        end = as.numeric(end)
    )
    
    qtl_df$Phenotype[qtl_df$molphe == "sQTL"] <- str_replace(qtl_df$Phenotype[qtl_df$molphe == "sQTL"], "^(chr\\d+):\\d+:\\d+:", "\\1:")
    
    qtl_df_dist = qtl_df %>% distinct()
    qtl_df_dist = qtl_df_dist %>% arrange(molphe, Phenotype, chr, start, end)
    
    range = fread(range_path)
    colnames(range) = c("chr", "start", "end")
    range$chr = paste0("chr", range$chr)
    
    qtl_df_dist_overlap = filter_molphe_by_range(qtl_df_dist, range)
    fwrite(qtl_df_dist_overlap, paste0("02_molphe_flt_", var_type, ".tsv"), sep = "\t")
}