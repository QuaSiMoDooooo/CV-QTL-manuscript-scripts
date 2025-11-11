#!/usr/bin/env Rscript

# -------------
# FileName     : methbank_CRMs.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Compare methylation levels between HZAU cohort and MethBank CRMs database
# -------------

library(dplyr)
library(data.table)

setwd("/home/wtian/project/HZAU_cohort_meth/05_extra_database/methbank/01_comp_CRMs/")

crm_path <- "/home/wtian/data/meth/methbank4.0/whole_blood_850K.txt.gz"
crm_df <- fread(crm_path, data.table = FALSE, header = TRUE, sep = "\t")

crm_df$average_level <- rowMeans(crm_df[,2:ncol(crm_df)], na.rm = TRUE)
crm_df$level <- ifelse(crm_df$average_level > 0.7, "high", 
                      ifelse(crm_df$average_level < 0.3, "low", "medium"))
crm_level <- crm_df[,c("probe","average_level","level")]
colnames(crm_level) <- c("probe_id","crm_average_level","crm_level")

data <- fread("meth_cpg_148.tsv", data.table = FALSE, header = TRUE, sep = "\t")
data$average_level <- rowMeans(data[,2:ncol(data)], na.rm = TRUE)
data$level <- ifelse(data$average_level > 0.7, "high", 
                    ifelse(data$average_level < 0.3, "low", "medium"))
data <- data[,c("methid","average_level","level")]

methloc <- fread("methloc.tsv", data.table = FALSE, header = TRUE, sep = "\t")
methloc$meth <- paste(methloc$chr, methloc$right, sep = ":")
merge_df <- left_join(data, methloc[,c("methid","meth")], by = "methid")

meth2cpg <- fread("merge_manifest.hg38.tsv", data.table = FALSE, header = TRUE, sep = "\t")
meth2cpg$meth <- paste(meth2cpg$chr, meth2cpg$hg38_start, sep = ":")
meth2cpg <- meth2cpg[,c("meth","probe_id")]
merge_df <- merge_df %>% left_join(meth2cpg, by = "meth")

fwrite(crm_level, "01_crm_cpg_average_level.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
fwrite(merge_df, "01_hzau_cpg_average_level.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

crm_level <- fread("01_crm_cpg_average_level.tsv", data.table = FALSE, header = TRUE, sep = "\t")
merge_df <- fread("01_hzau_cpg_average_level.tsv", data.table = FALSE, header = TRUE, sep = "\t")

merge_df2 <- merge_df %>% left_join(crm_level, by = "probe_id")

source('https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r')
library(ggpubr)
library(gridExtra)
library(tidyverse)

p1 <- ggscatter(merge_df2, x = "average_level", y = "crm_average_level",
                color = "#82B0D2", size = 2, alpha = 0.03, shape = 16, stroke = 0,
                add = "reg.line", conf.int = FALSE,
                add.params = list(color = "#FA7F6F")) +
  stat_cor(size = 3.6, fontface = "bold") +
  wdy_theme() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 14, face = "bold")
  ) +
  labs(x = "Methylation level of the cohort", y = "Methylation level of the CRMs")

p2 <- ggplot(merge_df2, aes(average_level)) +
  geom_density(fill = "#FFBE7A") +
  theme_void()

p3 <- ggplot(merge_df2, aes(crm_average_level)) +
  geom_density(fill = "#8ECFC9") +
  coord_flip() +
  theme_void()

library(patchwork)
empty_plot <- plot_spacer()
layout_design <- c("AAAAD", "BBBBC", "BBBBC", "BBBBC", "BBBBC")
p2 + p1 + p3 + empty_plot + plot_layout(design = layout_design)

ggsave("01_cor_scatter.pdf", width = 10, height = 10, units = "cm")
