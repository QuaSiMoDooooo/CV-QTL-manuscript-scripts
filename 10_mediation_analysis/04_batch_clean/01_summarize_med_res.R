#!/usr/bin/env Rscript

# -------------
# FileName     : 01_summarize_med_res.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Summarize mediation analysis results from batch processing
# -------------

suppressPackageStartupMessages(library(medflex))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/04_batch_clean")

mc = getOption("mc.cores", 100)

med_res_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/03_batch_mediation_analysis"
med_type_list = c("02_cis_QTL_overlap","02_QTL_GWAS_overlap")
outdir_list = paste0("01_",str_split_fixed(med_type_list, "_",2)[,2])
variant_type_list = c("SNP","SV")

for (i in 1:length(med_type_list)){
    med_type = med_type_list[i]
    outdir = outdir_list[i]
    if (!dir.exists(outdir)){dir.create(outdir)}
    for (j in 1:length(variant_type_list)){
        variant_type = variant_type_list[j]
        out_subdir = paste0(outdir,"/",variant_type)
        if (!dir.exists(out_subdir)){dir.create(out_subdir)}
        med_res_subdir = paste0(med_res_dir,"/",med_type,"/",variant_type)
        
        med_res_paths = list.files(med_res_subdir,pattern = "*.output",full.names = T)
        tmp = str_split(basename(med_res_paths),":", simplify = T)
        var_ids = tmp[,1]
        exp_ids = tmp[,2]
        med_ids = tmp[,3]
        
        s_index = grepl(":clu", med_res_paths)
        med_ids[s_index] = paste(tmp[s_index,3],tmp[s_index,4],tmp[s_index,5],tmp[s_index,6],sep = ":")
        
        results_list <- mclapply(1:length(med_res_paths), function(m) {
            tryCatch({
                med_res <- read.table(med_res_paths[m], header = TRUE, sep = "\t")
                med_res_long <- med_res %>%
                pivot_longer(
                    cols = colnames(med_res)[-1],
                    names_to = "stats",
                    values_to = "value"
                )
                list(
                TE = med_res_long$value[med_res_long$eff_names == "Total effect" & med_res_long$stats == "eff_est"],
                DE = med_res_long$value[med_res_long$eff_names == "Direct effect" & med_res_long$stats == "eff_est"],
                ME = med_res_long$value[med_res_long$eff_names == "Mediated effect" & med_res_long$stats == "eff_est"],
                TE_p = med_res_long$value[med_res_long$eff_names == "Total effect" & med_res_long$stats == "eff_pval"],
                DE_p = med_res_long$value[med_res_long$eff_names == "Direct effect" & med_res_long$stats == "eff_pval"],
                ME_p = med_res_long$value[med_res_long$eff_names == "Mediated effect" & med_res_long$stats == "eff_pval"],
                TE_lcl = med_res_long$value[med_res_long$eff_names == "Total effect" & med_res_long$stats == "eff_lcl"],
                DE_lcl = med_res_long$value[med_res_long$eff_names == "Direct effect" & med_res_long$stats == "eff_lcl"],
                ME_lcl = med_res_long$value[med_res_long$eff_names == "Mediated effect" & med_res_long$stats == "eff_lcl"],
                TE_hcl = med_res_long$value[med_res_long$eff_names == "Total effect" & med_res_long$stats == "eff_hcl"],
                DE_hcl = med_res_long$value[med_res_long$eff_names == "Direct effect" & med_res_long$stats == "eff_hcl"],
                ME_hcl = med_res_long$value[med_res_long$eff_names == "Mediated effect" & med_res_long$stats == "eff_hcl"]
                )
            }, error = function(e) {
                cat("Error processing file:", med_res_paths[m], "\n")
                cat("Error message:", e$message, "\n")
                return(NULL)
            })
        }, mc.cores = mc)
        
        TE <- unlist(lapply(results_list, function(x) x[["TE"]]))
        DE <- unlist(lapply(results_list, function(x) x[["DE"]]))
        ME <- unlist(lapply(results_list, function(x) x[["ME"]]))
        TE_p <- unlist(lapply(results_list, function(x) x[["TE_p"]]))
        DE_p <- unlist(lapply(results_list, function(x) x[["DE_p"]]))
        ME_p <- unlist(lapply(results_list, function(x) x[["ME_p"]]))
        TE_lcl <- unlist(lapply(results_list, function(x) x[["TE_lcl"]]))
        DE_lcl <- unlist(lapply(results_list, function(x) x[["DE_lcl"]]))
        ME_lcl <- unlist(lapply(results_list, function(x) x[["ME_lcl"]]))
        TE_hcl <- unlist(lapply(results_list, function(x) x[["TE_hcl"]]))
        DE_hcl <- unlist(lapply(results_list, function(x) x[["DE_hcl"]]))
        ME_hcl <- unlist(lapply(results_list, function(x) x[["ME_hcl"]]))
        
        df = data.frame(var = var_ids,med = med_ids,exp = exp_ids,
        TE = TE,DE = DE,ME = ME,
        TE_p = TE_p,DE_p = DE_p,ME_p = ME_p,
        TE_lcl = TE_lcl,DE_lcl = DE_lcl,ME_lcl = ME_lcl,
        TE_hcl = TE_hcl,DE_hcl = DE_hcl,ME_hcl = ME_hcl
        )
       fwrite(df,paste0(out_subdir,"/med_res_all.tsv"),quote = F,row.names = F,col.names = T, sep = "\t")
    }
}