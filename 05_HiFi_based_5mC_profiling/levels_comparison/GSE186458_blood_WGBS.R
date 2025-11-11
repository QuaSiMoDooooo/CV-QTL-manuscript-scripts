#!/usr/bin/env Rscript

# -------------
# FileName     : GSE186458_blood_WGBS.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Compare methylation levels between HZAU cohort and GSE186458 blood WGBS data
# -------------

library(data.table)
library(parallel)
library(dplyr)
library(vroom)

setwd("/home/wtian/project/HZAU_cohort_meth/05_extra_database/GSE186458")

bed_dir <- "~/data/meth/a_dna_methylation_atlas_of_normal_human_cell_types/02_wgbstools_beta2bed/253samples_hg38_bed"
bed_files <- list.files(bed_dir, pattern = "Blood.*\\.bed$", full.names = TRUE)

read_bed_unique <- function(file) {
  dt <- vroom(file, col_names = c("chr", "start", "end", "level"), delim = "\t", show_col_types = FALSE)
  dt <- dt %>% select(chr, start, level) %>% distinct(chr, start, .keep_all = TRUE)
  return(dt)
}

ncores <- detectCores() - 1
bed_list <- mclapply(bed_files, read_bed_unique, mc.cores = ncores)

merged_dt <- Reduce(function(x, y) inner_join(x, y, by = c("chr", "start")), bed_list)
merged_dt[merged_dt == -1] <- NA

merged_dt$average_level <- rowMeans(merged_dt[, -c(1:2)], na.rm = TRUE)
colnames(merged_dt) <- c("chr", "start", basename(bed_files), "average_level")
merged_dt <- merged_dt[!is.na(merged_dt$average_level), ]

fwrite(merged_dt, file = file.path(bed_dir, "Blood_samples_average_methylation.bed"), sep = "\t")

merged_dt <- fread("~/data/meth/a_dna_methylation_atlas_of_normal_human_cell_types/02_wgbstools_beta2bed/253samples_hg38_bed/Blood_samples_average_methylation.bed")

merged_df <- data.frame(
    chr_pos = paste(merged_dt$chr, merged_dt$start+1, sep = "_"),
    average_level = merged_dt$average_level
)

wb_meth_df_loc2 <- fread("/home/wtian/project/HZAU_cohort_meth/05_extra_database/ENCODE/01_comp_WGBS/02_wb_meth_level.tsv")
wb_meth_df_loc3 <- data.frame(
    chr_pos = paste(wb_meth_df_loc2$chr, wb_meth_df_loc2$pos, sep = "_"),
    wb_average_level = wb_meth_df_loc2$wb_average_level
)

merge_df <- inner_join(wb_meth_df_loc3, merged_df, by = "chr_pos")

source('https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r')
library(ggpubr)
library(gridExtra)
library(tidyverse)

cor_val <- cor.test(merge_df$wb_average_level, merge_df$average_level, method = "p")

sample_df <- merge_df %>% sample_n(100000)

p1 <- ggscatter(sample_df, x = "wb_average_level", y = "average_level",
           color = "#d8f1ff", size = 2, alpha = 0.05, shape = 16, stroke = 0,
           add = "reg.line", conf.int = FALSE,
           add.params = list(color = "#FA7F6F")) +
  wdy_theme() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 14, face = "bold")
  ) + 
  labs(x = "Methylation level of the cohort", y = "Methylation level of the WGBS")

p2 <- ggplot(sample_df, aes(wb_average_level)) +
  geom_density(fill = "#FFBE7A") +
  theme_void()

p3 <- ggplot(sample_df, aes(average_level)) +
  geom_density(fill = "#8ECFC9") +
  coord_flip() +
  theme_void()

library(patchwork)
empty_plot <- plot_spacer()
layout_design <- c("AAAAD", "BBBBC", "BBBBC", "BBBBC", "BBBBC")
p2 + p1 + p3 + empty_plot + plot_layout(design = layout_design)

ggsave("02_cor_scatter.pdf", width = 10, height = 10, units = "cm")