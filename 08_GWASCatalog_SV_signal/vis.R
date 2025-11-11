#!/usr/bin/env Rscript

# -------------
# FileName     : vis.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Generate visualization for GWAS Catalog SV signal analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(pryr))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/15_GWASCatalog_SV_signal")

df_full <- fread("gwas_flt_snp_ld_sv_info_gene_region.txt") %>% arrange(V2)
df_sub  <- fread("gwas_flt_snp_ld_sv_info_qtl_gene_region.txt") %>% arrange(V2)

region_levels <- c("intergenic","upstream","downstream","upstream;downstream",
                   "UTR5","UTR3","ncRNA_intronic","ncRNA_exonic",
                   "intronic","exonic","splicing")

df1 <- data.frame(
  Category = "SVs",
  Type = factor(df_full$V1, levels = region_levels),
  Count = df_full$V2
) %>%
  group_by(Category) %>%
  mutate(Percent = Count / sum(Count)) %>%
  ungroup()

df_ratio <- data.frame(
  Type = factor(df_full$V1, levels = region_levels),
  Full = df_full$V2,
  Sub = df_sub$V2
) %>%
  mutate(Ratio = Sub / Full) %>%
  arrange(Ratio) %>%
  mutate(Type = factor(Type, levels = Type))

custom_colors <- c(
  "intronic" = "#FB7F72", "exonic" = "#cc649a", "splicing" = "#990010",
  "ncRNA_intronic" = "#FFBD80", "ncRNA_exonic" = "#d37c0d",
  "UTR5" = "#80C5BF", "UTR3" = "#4a8f8a",
  "upstream" = "#82B1D1", "downstream" = "#4d7c9a", "upstream;downstream" = "#384955",
  "intergenic" = "#E7DBD2"
)

source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")

p1 <- ggplot(df1, aes(x = Category, y = Percent, fill = Type)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(add = c(0, 0.02))
  ) +
  labs(x = NULL, y = "Percentage of genomic regions (%)", fill = "Region") +
  wdy_theme() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

p2 <- ggplot(df_ratio, aes(x = Ratio, y = Type, fill = Type)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = custom_colors, guide = "none") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "Percentage of SV-QTLs (%)", y = NULL) +
  wdy_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_combined <- p1 + p2 + plot_layout(widths = c(1.5, 2))

ggsave("GWAS-catalog_SV_region_and_QTL_ratio.pdf", p_combined, width = 8, height = 4.5, dpi = 300)