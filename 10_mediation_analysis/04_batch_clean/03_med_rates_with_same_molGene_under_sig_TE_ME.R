#! Rscript

# -------------
# FileName     : 03_med_rates_with_same_molGene_under_sig_TE_ME.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Calculate mediation rates for same molecular gene under significant total and mediated effects
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/09_mediation_analysis/04_batch_clean")

data = fread("02_all_med_rates_under_sig_TE_ME.tsv", header = T, sep = "\t")

data = data[!grepl("QTL_GWAS_overlap",data$med_type)]
data = data[!grepl("meth",data$med),]
data$exp_gene = str_split(data$exp,"\.",simplify = T)[,1] # nolint: error.

index_s = grepl("clu", data$med)
data$med_gene = str_split(data$med, "_", simplify = T)[,1]

clu2ensembl = fread("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/08_trans_QTL_hotspot/05_transQTL_cisQTL_Genes/02_Genes_enrich/01_get_sGenes/02_clu_ensembl.tsv")
data$med_gene[index_s] = paste0(str_split(data$med[index_s], ":", simplify = T)[,1], ":", str_split(data$med[index_s], ":", simplify = T)[,4])
data$med_gene[index_s] = clu2ensembl$ENSEMBL[match(data$med_gene[index_s], clu2ensembl$clu)]

data_same_gene = data[data$exp_gene == data$med_gene,]

fwrite(data_same_gene, "03_med_rates_with_same_molGene_under_sig_TE_ME.tsv", sep = "\t", quote = F, row.names = F)

data_same_gene = fread("03_med_rates_with_same_molGene_under_sig_TE_ME.tsv", header = T, sep = "\t")

tmp = data.frame(
    type = paste0(str_split_fixed(data_same_gene$med_type,"_",2)[,2],":",data_same_gene$variant_type),
    med_rate = data_same_gene$med_rate
)

library(RColorBrewer)
mycolor <- brewer.pal(8, "Paired")[c(1,2)]
pdf("03_med_rates_with_same_molGene_under_sig_TE_ME_boxplot.pdf", width = 2.6, height = 4.7)
ggpubr::ggboxplot(tmp, x = "type", y = "med_rate",color = "type", palette = mycolor) +
labs(x = "Type", y = "Proportion of mediating effect") +
theme(
  axis.title.x = element_text(size = 14, face = "bold"),
  axis.title.y = element_text(size = 14, face = "bold"),
  axis.text.x = element_text(size = 12, face = "bold", angle = 45, hjust = 1),
  axis.text.y = element_text(size = 12, face = "bold"),
  legend.title = element_text(size = 10),
  legend.text = element_text(size = 10),
  legend.position = "none")
dev.off()