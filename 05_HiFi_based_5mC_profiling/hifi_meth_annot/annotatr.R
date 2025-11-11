#! Rscript
# -------------
# FileName     : annotatr
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : annotate HiFi-based 5mC sites with their genomic context
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
suppressPackageStartupMessages(library(annotatr))
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)

dm_file = "methloc.bed"

options(scipen = 999)
tmp = fread(dm_file, header = FALSE)
fwrite(tmp, file = "methloc.bed", sep = "\t", row.names = FALSE, col.names = FALSE)

dm_regions = read_regions(con = dm_file, genome = 'hg38', format = 'bed')
annots = c('hg38_cpgs','hg38_basicgenes','hg38_genes_intergenic')
annotations = build_annotations(genome = 'hg38', annotations = annots)

unique_types <- unique(annotations$type)
ind_reg = c("hg38_cpg_islands","hg38_cpg_shores","hg38_cpg_shelves","hg38_cpg_inter",
            "hg38_genes_1to5kb","hg38_genes_promoters","hg38_genes_5UTRs","hg38_genes_exons","hg38_genes_introns","hg38_genes_3UTRs","hg38_genes_intergenic")

type_length_sum <- integer(length = length(unique_types))
for (i in seq_along(unique_types)) {
  current_type <- unique_types[i]
  current_ranges <- annotations[annotations$type == current_type, ]
  type_length_sum[i] <- sum(width(reduce(current_ranges)@ranges))
}

type_length_sum_df <- data.frame(
  type = unique_types,
  total_length = type_length_sum
)

type_length_sum_df$category <- ifelse(grepl("cpg",type_length_sum_df$type),"cpg","gene")

p <- ggplot(type_length_sum_df, aes(x = type , y = total_length, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.8) +
  geom_text(aes(label = total_length),position = position_dodge(width = 0.5),vjust=-0.5, size = 4)+
  labs(y = "total length (bp)") +
  scale_y_continuous(limits = c(0,max(type_length_sum_df$total_length)+100000000),expand = c(0, 0)) +
  scale_fill_manual(values = c("cpg" = "#80C5BF", "gene" = "#F9B27C")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", size = 1, linetype = "solid"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, margin = margin(r = 10)),
    axis.text.x = element_text(size=10,angle = 45, vjust = 0.5, hjust = 1, margin = margin(t = -42)),
    axis.text.y = element_text(size=15),
    plot.margin = unit(c(1,1,3,1),"lines"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )
ggsave(p, filename = "annot_regions_length_uniq.pdf",width = 12,height = 8)

dm_annotated = annotate_regions(
    regions = dm_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
fwrite(df_dm_annotated,"dm_annotated.tsv",sep = "\t",quote = FALSE,row.names = FALSE)
df_dm_annotated = fread("dm_annotated.tsv",header = TRUE)

df_dm_annotated_symbol = df_dm_annotated[!is.na(df_dm_annotated$annot.symbol),]
fwrite(df_dm_annotated_symbol,"02_dm_annotated_symbol.tsv",sep = "\t",quote = FALSE,row.names = FALSE)

dm_annsum = summarize_annotations(
    annotated_regions = dm_annotated,
    quiet = TRUE)
dm_annsum$category <- ifelse(grepl("cpg",dm_annsum$annot.type),"cpg","gene")
dm_annsum$annot.type = factor(dm_annsum$annot.type, levels = ind_reg)

p1 <- ggplot(dm_annsum, aes(x = annot.type , y = n, fill = category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.8) +
  geom_text(aes(label = n),position = position_dodge(width = 0.5),vjust=-0.5, size = 5)+
  xlab("annotation") +
  ylab("counts") +
  scale_y_continuous(expand=expansion(add = c(0, 1)), limits = c(0,max(dm_annsum$n)+1000000)) +
  expand_limits(y = c(0)) +
  scale_fill_manual(values = c("cpg" = "#80C5BF", "gene" = "#F9B27C")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", size = 1, linetype = "solid"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, margin = margin(r = 10)),
    axis.text.x = element_text(size=10,angle = 45, vjust = 0.5, hjust = 1, margin = margin(t = -42)),
    axis.text.y = element_text(size=15),
    plot.margin = unit(c(1,1,3,1),"lines"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )
ggsave(p1, filename = "annot_cpg_gene_counts.pdf",width = 12,height = 8)

names(type_length_sum_df)[1] = "annot.type"
df_merge = left_join(type_length_sum_df, dm_annsum, by = "annot.type")
df_merge$percentage = df_merge$n / df_merge$total_length * 100
df_merge$percentage = format(round(df_merge$percentage, 2), nsmall = 2)
df_merge$percentage = as.numeric(df_merge$percentage)
df_merge$annot.type = factor(df_merge$annot.type, levels = ind_reg)

p2 <- ggplot(df_merge, aes(x = annot.type , y = percentage, fill = category.y)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.3), width = 0.8) +
  geom_text(aes(label = percentage),position = position_dodge(width = 0.5),vjust=-0.5, size = 5)+
  ylab("percent of methylation sites (%)") +
  scale_y_continuous(expand=expansion(add = c(0, 1)), limits = c(0,max(df_merge$percentage)+0.005)) +
  expand_limits(y = c(0)) +
  scale_fill_manual(values = c("cpg" = "#80C5BF", "gene" = "#F9B27C")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.border = element_rect(color = "black", size = 1, linetype = "solid"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 24, margin = margin(r = 10)),
    axis.text.x = element_text(size=10,angle = 45, vjust = 0.5, hjust = 1, margin = margin(t = -42)),
    axis.text.y = element_text(size=15),
    plot.margin = unit(c(1,1,3,1),"lines"),
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16)
  )+
  labs(fill = "category")
ggsave(p2, filename = "annot_cpg_gene_percentage.pdf",width = 6.5,height = 4)

ind_reg = c("hg38_cpg_islands","hg38_cpg_shores","hg38_cpg_shelves","hg38_cpg_inter",
            "hg38_genes_1to5kb","hg38_genes_promoters","hg38_genes_5UTRs","hg38_genes_exons","hg38_genes_introns","hg38_genes_3UTRs","hg38_genes_intergenic")
df_dm_annotated = fread("dm_annotated.tsv",header = TRUE)
methlevel = fread("meth_average.tsv")
methloc = fread("HZAU_meth_148_imputed_flted/methloc.tsv")
methloc_level = left_join(methloc, methlevel, by = "methid")

df_annot = data.frame(
  chr = df_dm_annotated$seqnames,
  left = df_dm_annotated$start-1,
  right = df_dm_annotated$end,
  annot.type = df_dm_annotated$annot.type
)

region_avg_meth = left_join(df_annot, methloc_level, by = c("chr", "left", "right"))
region_avg_meth$annot.type = factor(region_avg_meth$annot.type, levels = ind_reg)
region_avg_meth$category = ifelse(grepl("cpg",region_avg_meth$annot.type), "cpg", "gene")
region_avg_meth$mean = round(region_avg_meth$mean*100, 2)
fwrite(region_avg_meth, "region_avg_meth.tsv", sep = "\t", quote = FALSE)
region_avg_meth = fread("region_avg_meth.tsv")
region_avg_meth$annot.type = factor(region_avg_meth$annot.type, levels = ind_reg)
p = ggplot(region_avg_meth, aes(x = annot.type, y = mean, fill = category)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c("cpg" = "#7FC5BF", "gene" = "#F9B280")) +
  labs(y = "DNA Methylation level (%)") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
ggsave(p, filename = "annot_cpg_gene_methlevel.pdf",width = 6.5,height = 4)
