#!/usr/bin/env Rscript

# -------------
# FileName     : 02_med_res_all_stats.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Statistical analysis of mediation results including effect significance and mediation rates
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/04_batch_clean")

med_type_list = c("01_cis_QTL_overlap","01_QTL_GWAS_overlap")
variant_type_list = c("SNP","SV")

all_stats = data.frame()
for (m in 1:length(med_type_list)){
    med_type = med_type_list[m]
    for (j in 1:length(variant_type_list)){
        variant_type = variant_type_list[j]
        med_res_all_path = paste0(med_type,"/",variant_type,"/med_res_all.tsv")
        med_res_all = fread(med_res_all_path)
        
        all_events = dim(med_res_all)[1]
        sig_TE = sum(med_res_all$TE_p < 0.05)
        sig_DE = sum(med_res_all$DE_p < 0.05)
        sig_ME = sum(med_res_all$ME_p < 0.05)
        
        sig_TE_med_res = med_res_all[TE_p < 0.05]
        sig_TE_DE = sum(sig_TE_med_res$DE_p < 0.05)
        sig_TE_ME = sum(sig_TE_med_res$ME_p < 0.05)
        sig_TE_DE_ME = sum(sig_TE_med_res$DE_p < 0.05 & sig_TE_med_res$ME_p < 0.05)
        
        single_stats = data.frame(
            med_type = med_type, variant_type = variant_type,
            all_events = all_events, sig_TE = sig_TE, sig_DE = sig_DE, sig_ME = sig_ME,
            sig_TE_DE = sig_TE_DE, sig_TE_ME = sig_TE_ME, sig_TE_DE_ME = sig_TE_DE_ME
        )
        all_stats = rbind(all_stats,single_stats)
    }
}
write.table(all_stats,"02_all_stats.tsv",quote = F,row.names = F,sep = "\t")

check_all_positive_or_negative <- function(x) {
  all_positive <- all(x > 0)
  all_negative <- all(x < 0)
  return(all_positive | all_negative)
}

all_med_rates = data.frame()
for (m in 1:length(med_type_list)){
    med_type = med_type_list[m]
    for (j in 1:length(variant_type_list)){
        variant_type = variant_type_list[j]
        med_res_all_path = paste0(med_type,"/",variant_type,"/med_res_all.tsv")
        med_res_all = fread(med_res_all_path)
        sig_TE_ME = med_res_all[med_res_all$TE_p < 0.05 & med_res_all$ME_p < 0.05]
        
        judge_direct = unlist(lapply(1:dim(sig_TE_ME)[1], function(x) {
            tmp = sig_TE_ME[x,]
            check_all_positive_or_negative(c(tmp$TE,tmp$DE,tmp$ME))
        }))
        
        sig_TE_ME_same_direct = sig_TE_ME[judge_direct,]
        
        med_rate = unlist(lapply(1:dim(sig_TE_ME_same_direct)[1], function(x) {
            tmp = sig_TE_ME_same_direct[x,]
            abs(tmp$ME)/abs(tmp$TE)
        }))
        
        single_med_rate = data.frame(
            med_type = med_type, variant_type = variant_type,
            var = sig_TE_ME_same_direct$var, med = sig_TE_ME_same_direct$med,  exp = sig_TE_ME_same_direct$exp,
            med_rate = med_rate
        )
        all_med_rates = rbind(all_med_rates,single_med_rate)
    }
}

write.table(all_med_rates,"02_all_med_rates_under_sig_TE_ME.tsv",quote = F,row.names = F,sep = "\t")

all_med_rates_grouped <- all_med_rates %>%
  group_by(med_type, variant_type) %>%
  summarise(
    max_med_rate = max(med_rate),
    min_med_rate = min(med_rate)
  ) %>% as.data.frame()

write.table(all_med_rates_grouped,"02_all_med_rates_under_sig_TE_ME_grouped.tsv",quote = F,row.names = F,sep = "\t")

library(RColorBrewer)
mycolor <- brewer.pal(8, "Paired")[c(1,2,5,6)]
pdf("02_all_med_rates_under_sig_TE_ME_boxplot.pdf", width = 4, height = 5)
ggpubr::ggboxplot(tmp, x = "Type", y = "med_rate", color = "Type",  palette = mycolor) +
labs(x = "Type", y = "Proportion of mediating effect") +
theme(
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.title.y = element_text(size = 14, face = "bold"),
  axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.position = "none"
)
dev.off()