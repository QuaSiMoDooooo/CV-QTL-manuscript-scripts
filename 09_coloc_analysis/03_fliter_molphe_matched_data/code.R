#!/usr/bin/env Rscript

# -------------
# FileName     : code.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Filter molecular phenotype matched data for colocalization analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))

options(stringsAsFactors = FALSE)
options(scipen = 999)
mc <- getOption('mc.cores', 36)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/03_fliter_molphe_matched_data")

flt_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/02_get_molphe_within_sentinelVARs_range"
flt_files = list.files(flt_dir, pattern = "*.tsv", full.names = TRUE)
var_type_list = sub(".tsv","",str_split(basename(flt_files), "_", simplify = TRUE)[, 4])

molphe_names = c("expression","splicing","APA","methylation")
qtl_names = c("eQTL","sQTL","apaQTL","meQTL")
matched_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA/"

cov_filename = "covariant.txt"
phe_pos = "phenotype.position.txt"
phe_qua = "phenotype.quantity.matched.txt"
var_pos = "variant.position.txt"
var_qua = "variant.genotype.matched.txt"

flt_molphe <- function(flt_file, qtl_name, input_dir, outdir) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    cov_file = paste0(input_dir, "/", cov_filename)
    phe_pos_file = paste0(input_dir, "/", phe_pos)
    phe_qua_file = paste0(input_dir, "/", phe_qua)
    var_pos_file = paste0(input_dir, "/", var_pos)
    var_qua_file = paste0(input_dir, "/", var_qua)
    system(paste0("cp ", cov_file, " ", outdir))
    system(paste0("cp ", var_pos_file, " ", outdir))
    system(paste0("cp ", var_qua_file, " ", outdir))

    flt_df = fread(flt_file)
    flt_molphe_df = flt_df %>% filter(molphe == qtl_name)
    phe_pos_df = fread(phe_pos_file)
    phe_qua_df = fread(phe_qua_file)

    if (qtl_name == "sQTL") {
        flt_molphe_df$Phenotype = paste(
            flt_molphe_df$chr,flt_molphe_df$start,flt_molphe_df$end,str_split(flt_molphe_df$Phenotype, ":", simplify = TRUE)[,2],
            sep = ":"
        )
    }
    
    phe_qua_df = phe_qua_df %>% filter(ID %in% flt_molphe_df$Phenotype)
    phe_pos_df = phe_pos_df %>% filter(ID %in% flt_molphe_df$Phenotype)
    fwrite(phe_pos_df, paste0(outdir, "/", phe_pos), sep = "\t")
    fwrite(phe_qua_df, paste0(outdir, "/", phe_qua), sep = "\t")
}

for (j in 1:length(var_type_list)) {
    var_type = var_type_list[j]
    flt_file = flt_files[j]
    for (m in 1:length(molphe_names)) {
        molphe_name = molphe_names[m]
        qtl_name = qtl_names[m]
        input_dir = paste0(matched_dir, var_type, "-", molphe_name, "/matched_data")
        outdir = paste0("results/",var_type, "-", molphe_name, "/matched_data")
        flt_molphe(flt_file, qtl_name, input_dir, outdir)
    }
}