#!/usr/bin/env Rscript

# -------------
# FileName     : code.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Extract GWAS and eQTL sentinel variants within 1Mb ranges for colocalization analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(xQTLbiolinks))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/01_get_chr1-22_GWAS_and_eQTL_sentinelVARs_range")

mc <- getOption("mc.cores", 200)

gwas_flt_dir = "../data/GWAS_flt"
qtl_flt_dir = "../data/QTL_flt"

gwas_flt_files = list.files(gwas_flt_dir, pattern = "\\.tsv$", full.names = TRUE)
gwas_var_type_list = str_split(basename(gwas_flt_files), "\\.", simplify = TRUE)[, 1]

qtl_flt_files = list.files(qtl_flt_dir, pattern = "\\.tsv$", full.names = TRUE)
qtl_var_type_list = str_split(basename(qtl_flt_files), "_", simplify = TRUE)[, 3]

eff_qtl_var_type_list = intersect(gwas_var_type_list, qtl_var_type_list)

merge_sentinel_ranges <- function(df) {
    if (!is.data.table(df)) df <- as.data.table(df)
    df[, position := as.integer(position)]
    df <- df[order(chr, position)]
    df[, start := pmax(1, position - 1e6)]
    df[, end := position + 1e6]
    gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start, end = df$end))
    gr_merged <- reduce(gr)
    result <- data.frame(chr = as.character(seqnames(gr_merged)),
                        start = start(gr_merged),
                        end = end(gr_merged))
    result$chr <- sub("chr", "", result$chr)
    return(result)
}

print("GWAS")
for (j in 1:length(eff_qtl_var_type_list)) {
    var_type = eff_qtl_var_type_list[j]
    print(var_type)
    gwas_flt_file = gwas_flt_files[gwas_var_type_list == var_type]
    gwas_flt = fread(gwas_flt_file)
    print(paste0("effective trait numbers: ", length(unique(gwas_flt$Indicator))))
    tarits = unique(gwas_flt$Indicator)
    mclapply(1:length(tarits), function(t) {
        tryCatch({
            trait = tarits[t]
            outdir = paste0("GWAS_",var_type,"/")
            if (!dir.exists(outdir)) dir.create(outdir)
            outpath = paste0(outdir, trait, "_sentinelVARs_1Mb_range.bed")
            gwas_flt_trait = gwas_flt[gwas_flt$Indicator == trait, ]
            gwas_flt_trait$AF = 0.3
            gwas_format = select(gwas_flt_trait, c("SNP","CHR","BP","P","AF","BETA","SE"))
            colnames(gwas_format) <- c("rsid","chrom","position","pValue","AF","beta","se")
            if (nrow(gwas_format) == 1){
                sentinelVARs_bed = data.frame(gwas_format$chrom, max(1, as.numeric(gwas_format$position)-1e6), as.numeric(gwas_format$position)+1e6)
                colnames(sentinelVARs_bed) = c("chr","start","end")
                sentinelVARs_bed$chr = sub("chr", "", sentinelVARs_bed$chr)
            } else {
                sentinelVARs_df = suppressMessages(xQTLanalyze_getSentinelSnp(gwas_format, pValueThreshold = 1, mafThreshold = 0, centerRange = 1e6))
                sentinelVARs_bed = merge_sentinel_ranges(sentinelVARs_df)
            }
            write.table(sentinelVARs_bed, file = outpath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
        }, error = function(e) {
            message(paste0("GWAS - Error at index ", t, " (trait: ", tarits[t], "): ", e$message))
        })
    }, mc.cores = mc)
}

print("eQTL")
for (j in 1:length(eff_qtl_var_type_list)) {
    var_type = eff_qtl_var_type_list[j]
    print(var_type)
    qtl_flt_file = qtl_flt_files[qtl_var_type_list == var_type]
    qtl_flt = fread(qtl_flt_file)
    eqtl_flt = qtl_flt[qtl_flt$molphe == "eQTL", ]
    eqtl_flt = eqtl_flt[!grepl("chrX", eqtl_flt$Variant_position) & !grepl("chrX", eqtl_flt$Phenotype_position), ]
    print(paste0("effective gene numbers: ", length(unique(eqtl_flt$Gene))))
    genes = unique(eqtl_flt$Gene)
    mclapply(1:length(genes), function(t) {
        tryCatch({
            gene = genes[t]
            outdir = paste0("eQTL_",var_type,"/")
            if (!dir.exists(outdir)) dir.create(outdir)
            outpath = paste0(outdir, gene, "_sentinelVARs_1Mb_range.bed")
            eqtl_flt_gene = eqtl_flt[eqtl_flt$Gene == gene, ]
            eqtl_flt_gene$AF = 0.3
            tmp = str_split(eqtl_flt_gene$Variant_position, ":", simplify = TRUE)
            eqtl_flt_gene$CHR = tmp[, 1]
            eqtl_flt_gene$BP = str_split(tmp[, 2], "-", simplify = TRUE)[, 1]
            eqtl_format = select(eqtl_flt_gene, c("Variant","CHR","BP","P","AF","Beta","Se"))
            colnames(eqtl_format) <- c("rsid","chrom","position","pValue","AF","beta","se")
            if (nrow(eqtl_format) == 1){
                sentinelVARs_bed = data.frame(eqtl_format$chrom, max(1, as.numeric(eqtl_format$position)-1e6), as.numeric(eqtl_format$position)+1e6)
                colnames(sentinelVARs_bed) = c("chr","start","end")
                sentinelVARs_bed$chr = sub("chr", "", sentinelVARs_bed$chr)
            } else {
                sentinelVARs_df = suppressMessages(xQTLanalyze_getSentinelSnp(eqtl_format, pValueThreshold = 1, mafThreshold = 0, centerRange = 1e6))
                sentinelVARs_bed = merge_sentinel_ranges(sentinelVARs_df)
            }
            write.table(sentinelVARs_bed, file = outpath, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
        }, error = function(e) {
            message(paste0("eQTL - Error at index ", t, " (gene: ", genes[t], "): ", e$message))
        })
    }, mc.cores = mc)
}