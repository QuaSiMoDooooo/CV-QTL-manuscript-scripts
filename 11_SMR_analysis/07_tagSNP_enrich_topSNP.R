#!/usr/bin/env Rscript

# -------------
# FileName     : 07_tagSNP_enrich_topSNP.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Analyze tagSNP enrichment in top SNPs from SMR analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA/"
qtl_name = c("expression","splicing","APA","methylation")

qtl_path = paste0(qtl_dir,"SNP-",qtl_name,"/QTL_results/cis.filtered.txt.gz")
qtl = do.call(
  rbind,
  lapply(
    qtl_path,
    function(x) {
      fread(x,showProgress = FALSE)
    }
  )
)

tmp = str_split(unique(qtl$Variant), "_" , simplify = TRUE)
uniq_cis_qtl = paste0(tmp[,1],"_",tmp[,2])
uniq_cis_qtl_num = length(uniq_cis_qtl)

ld = fread("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/02_complexVar_SNP_LD/05_ld_process/01_ld_subsets/SNP_SV.ld", sep=",")
ld = ld[abs(ld$R2) > 0.75,]
tag_snp = paste0("chr",ld$CHR_A,"_",ld$BP_A)
uniq_tag_snp = unique(tag_snp)
uniq_tag_snp_num = length(uniq_tag_snp)

smr = fread("05_smr_all_cis_flt.tsv")
uniq_top_snp = unique(paste0("chr",smr$topSNP_chr,"_",smr$topSNP_bp))

library(ggplot2)
library(reshape2)

total_snp <- 1861465
tag_total <- 488049
top_total <- 220
tag_top <- 22

mat <- matrix(c(
  tag_top,
  top_total - tag_top,
  tag_total - tag_top,
  (total_snp - tag_total) - (top_total - tag_top)
), nrow = 2, byrow = TRUE)
dimnames(mat) <- list(TopSNP = c("Yes", "No"), TagSNP = c("Yes", "No"))

fisher_result <- fisher.test(mat)
p_text <- paste0("Fisher's exact p = ", signif(fisher_result$p.value, 3))

df <- data.frame(
  Group = c("All SNPs", "Top SNPs"),
  TagSNP = c(tag_total, tag_top),
  NonTagSNP = c(total_snp - tag_total, top_total - tag_top)
)

df_long <- melt(df, id.vars = "Group")
colnames(df_long) <- c("Group", "SNP_Type", "Count")
df_long$Proportion <- ave(df_long$Count, df_long$Group, FUN = function(x) x / sum(x))

ggplot(df_long, aes(x = Group, y = Proportion, fill = SNP_Type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("TagSNP" = "#E69F00", "NonTagSNP" = "#56B4E9")) +
  labs(
    title = "Proportion of tagSNPs in All SNPs vs Top SNPs",
    y = "Proportion", x = "", fill = "SNP Type"
  ) +
  annotate("text", x = 1.5, y = 1.05, label = p_text, size = 5, fontface = "italic") +
  theme_minimal(base_size = 14) +
  theme(plot.margin = margin(10, 20, 10, 10))