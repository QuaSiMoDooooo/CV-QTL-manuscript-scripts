#!/usr/bin/env Rscript

# -------------
# FileName     : 02_QTL_GWAS_overlap_batch.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Batch mediation analysis for QTL and GWAS overlap data
# -------------

suppressPackageStartupMessages(library(medflex))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/03_batch_mediation_analysis")

data_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/02_format/02_QTL_GWAS_overlap"
outdir = paste0("02_",str_split_fixed(basename(data_dir), "_",2)[,2])
variant_type_list = c("SNP","SV")

mc = getOption("mc.cores",100)

med_analysis <- function(data_path, out_dir) {
    if (!dir.exists(out_dir)) {
        dir.create(out_dir)
    }
    resfile = paste0(out_dir, "/",sub("input","output",basename(data_path)))
    if (!file.exists(resfile)) {
        data = fread(data_path)
        data = na.omit(data)
        
        expData = neImpute(e_values~var_values+m_values+gender+age, data=data,family=gaussian)
        
        tryCatch({
            nMod = neModel(e_values ~ var_values0 + var_values1 + gender + age, expData, family = gaussian, se = "robust")
            summary(nMod)
            neEffdecomp(nMod)
            
            cf = as.data.frame(coef(summary(neEffdecomp(nMod))))
            ci = confint(neEffdecomp(nMod))
            eff_names = c("Direct effect", "Mediated effect", "Total effect")
            eff_est = cf$Estimate
            eff_lcl = unname(ci[, 1])
            eff_hcl = unname(ci[, 2])
            eff_pval = cf$`Pr(>|z|)`
            
            df = data.frame(
                eff_names = factor(eff_names, levels = c("Total effect", "Mediated effect", "Direct effect")),
                eff_est = eff_est,
                eff_lcl = eff_lcl,
                eff_hcl = eff_hcl,
                eff_pval = eff_pval
            )
            
            write.table(df, resfile, col.names = TRUE, sep = "\t")
            
            sig_res_plot = paste0(out_dir, "/sig_res_plot")
            if (!dir.exists(sig_res_plot)) { dir.create(sig_res_plot) }
            
            plot_file = paste0(sig_res_plot, "/", sub("input", "plot.pdf", basename(data_path)))
            
            if (any(df$eff_pval < 0.05) & !file.exists(plot_file)) {
                p = ggplot(df, aes(x = eff_est, y = eff_names)) +
                    geom_point(color = "#0872AF") +
                    geom_errorbarh(aes(xmin = eff_lcl, xmax = eff_hcl), height = 0.2, color = "#0872AF") +
                    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
                    geom_text(aes(x = eff_hcl + 0.02, label = sprintf("%.3e", eff_pval)), hjust = 0, size = 4, color = "black") + 
                    labs(x = "Effect estimate", y = "") +
                    xlim(min(df$eff_lcl) - 0.1, max(df$eff_hcl) + 0.3) +
                    theme_minimal() +
                    theme(
                        axis.text.y = element_text(size = 14, face = "bold"),
                        axis.title.x = element_text(size = 14, face = "bold"),
                        axis.title.y = element_blank(),
                        panel.grid = element_blank(),
                        axis.line = element_line(color = "black", size = 0.5),
                        axis.text.x = element_text(size = 14)
                    )
                ggsave(p, filename = plot_file, width = 6, height = 4, units = "in", dpi = 300)
            }
        }, error = function(e) {
            print(paste0("Failed to analyze ", basename(data_path)))
            cat("An error occurred: ", conditionMessage(e), "\n")
        })
    }
}

for (j in 1:length(variant_type_list)) {
    variant_type = variant_type_list[j]
    input_dir = paste0(data_dir, "/", variant_type)
    out_dir = paste0(outdir, "/", variant_type)
    if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = T)
    }
    data_paths = list.files(input_dir, pattern = "*.input", full.names = T)
    mclapply(data_paths, function(data_path) {
        med_analysis(data_path, out_dir)
    }, mc.cores = mc)
}