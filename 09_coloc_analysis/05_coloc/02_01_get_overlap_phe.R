#!/usr/bin/env Rscript

# -------------
# FileName     : 02_01_get_overlap_phe.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Identify molecular phenotypes overlapping with QTL sentinel variant ranges for colocalization analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(xQTLbiolinks))
suppressPackageStartupMessages(library(parallel))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc")

mc <- getOption("mc.cores", 5)

range_dir <- "../01_get_chr1-22_GWAS_and_eQTL_sentinelVARs_range/"

qtl_flt_dir <- "../data/QTL_flt"
qtl_flt_files <- list.files(qtl_flt_dir, pattern = "\\.tsv$", full.names = TRUE)
qtl_var_type_list <- str_split(basename(qtl_flt_files), "_", simplify = TRUE)[, 3]
var_type_list <- c("SNP","SV")
qtl_flt_files <- qtl_flt_files[qtl_var_type_list %in% var_type_list]
qtl_var_type_list <- qtl_var_type_list[qtl_var_type_list %in% var_type_list]

qtl_dir <- "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/04_rerun_QTL/results"
molphe_names <- c("splicing","APA","methylation")

filter_molphe_by_range <- function(df, range) {
    range_gr <- GRanges(
        seqnames = range$chr,
        ranges = IRanges(start = range$start, end = range$end)
    )
    qtl_gr <- GRanges(
        seqnames = df$chr,
        ranges = IRanges(start = df$start, end = df$end)
    )
    hits <- findOverlaps(qtl_gr, range_gr)
    df_overlap <- df[queryHits(hits), ]
    return(df_overlap)
}

for (j in 1:length(qtl_var_type_list)) {
    var_type = qtl_var_type_list[j]
    print(var_type)

    for (m in 1:length(molphe_names)) {
        molphe_name = molphe_names[m]
        print(molphe_name)
        mc <- getOption("mc.cores", 60)
        if (molphe_name == "methylation"){ mc <- getOption("mc.cores", 6) }

        overlap_phe_dir = paste0("02_01_overlap_phe/",var_type,"-",molphe_name)
        if (!dir.exists(overlap_phe_dir)){ dir.create(overlap_phe_dir, recursive = TRUE)}

        qtl_path = paste0(qtl_dir,"/",var_type,"-",molphe_name,"/QTL_results/cis.txt")
        qtl = data.table(vroom(qtl_path, delim = "\t", show_col_types = FALSE, progress = FALSE))
        tmp = str_split(str_split(qtl$Phenotype_position, ":", simplify = TRUE)[,2],"-", simplify = TRUE)
        start = tmp[,1]
        end = tmp[,2]
        qtl_df = data.frame(
            Phenotype = qtl$Phenotype,
            chr = str_split(qtl$Phenotype_position, ":", simplify = TRUE)[,1],
            start = as.numeric(start), end = as.numeric(end)
        )
        rm(tmp, start, end)
        gc()

        qtl_df_dist = qtl_df %>% distinct() %>% arrange(Phenotype, chr, start, end)
        rm(qtl_df)
        gc()

        range_subdir = paste0(range_dir,"eQTL_",var_type)
        eqtl_flt_file = qtl_flt_files[j]
        eqtl_flt = data.table(vroom(eqtl_flt_file, delim = "\t", show_col_types = FALSE, progress = FALSE))
        eqtl_flt = eqtl_flt %>% filter(molphe == "eQTL")
        genes = str_split(unique(eqtl_flt$Phenotype), "\\.", simplify = TRUE)[, 1]
        print(paste0("effective gene numbers:", length(genes)))
        rm(eqtl_flt)
        gc()

        mclapply(1:length(genes), function(g) {
            tryCatch({
                overlap_phe_path = paste0(overlap_phe_dir,"/",genes[g],".txt")
                if (file.exists(overlap_phe_path)) { return(NULL) }
                gene = genes[g]
                sentinelVARs_range_path = paste0(range_subdir,"/",gene,"_sentinelVARs_1Mb_range.bed")
                sentinelVARs_range = read.table(sentinelVARs_range_path, sep = "\t", header = F)
                colnames(sentinelVARs_range) = c("chr", "start", "end")
                sentinelVARs_range$chr = paste0("chr", sentinelVARs_range$chr)
                
                overlap_phe = filter_molphe_by_range(qtl_df_dist, sentinelVARs_range)$Phenotype
                if (length(overlap_phe) > 0) {
                    overlap_phe = unique(overlap_phe)
                    write.table(overlap_phe, file = overlap_phe_path, quote = FALSE, row.names = FALSE, col.names = FALSE)
                }
            }, error = function(e) {
                cat("Error processing var_type:", var_type, "molphe_name:", molphe_name, "gene:", genes[g], "\n")
                cat("Error message:", conditionMessage(e), "\n")
            })
        }, mc.cores = mc)
    }
}