#!/usr/bin/env Rscript

# -------------
# FileName     : 05_med_TF.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Transcription factor enrichment analysis for mediation results
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(patchwork))

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

data_flt = data[!data$exp_gene == data$med_gene,]

write.table(unique(data_flt$med_gene), "05_med_no_same_molGene.txt", quote = F, row.names = F, col.names = F, sep = "\t")

gene_ensembls = unique(data_flt$med_gene)
library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_positions <- getBM(
  attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position"),
  filters = "ensembl_gene_id",
  values = gene_ensembls,
  mart = ensembl
)
gene_positions$start_position <- as.numeric(gene_positions$start_position)
gene_positions$end_position <- as.numeric(gene_positions$end_position)
filtered_genes <- subset(
  gene_positions,
  !(chromosome_name == "6" & start_position < 35000000 & end_position > 25000000)
)
filtered_gene_ids <- filtered_genes$ensembl_gene_id
write.table(filtered_gene_ids, "05_med_no_same_molGene_noMHC.txt", quote = F, row.names = F, col.names = F, sep = "\t")

df = fread("05_med_no_same_molGene_noMHC_gprofiler.csv")

plot_enrichment_dotplot <- function(df, output_pdf = "enrichment_dotplot.pdf") {
    set.seed(123)
    source_levels <- c("GO:MF", "GO:BP", "GO:CC", "KEGG", "REAC", "WP", "TF", "MIRNA", "HP")
    source_colors <- c(
      "GO:MF" = "#d73027", 
      "GO:BP" = "#fc8d59",
      "GO:CC" = "#1a9850", 
      "KEGG"  = "#df65b0", 
      "REAC"  = "#7570b3", 
      "WP"    = "#1c91c0", 
      "TF"    = "#4c72b0", 
      "MIRNA" = "#1b9e77", 
      "HP"    = "#9e14ea"
    )
    df <- df %>%
    dplyr::select(term_name, source, adjusted_p_value) %>%
    mutate(
        log10_padj = -log10(adjusted_p_value),
        source = factor(source, levels = source_levels)
    )
    label_df <- df %>%
    group_by(source) %>%
    group_modify(~ slice_head(.x[order(.x$adjusted_p_value), ], n = min(2, nrow(.x)))) %>%
    ungroup()

    y_range <- range(df$log10_padj, na.rm = TRUE)
    y_height <- diff(y_range) * 0.01
    p <- ggplot(df, aes(x = source, y = log10_padj)) +
    geom_tile(
      data = data.frame(source = factor(source_levels, levels = source_levels)),
      aes(x = source, y = 0),
      fill = source_colors[source_levels],
      height = y_height,
      inherit.aes = FALSE
    ) +
    geom_jitter(aes(color = source), width = 0.2, size = 2, alpha = 0.7, show.legend = FALSE) +
    geom_text_repel(data = label_df, aes(label = term_name, color = source,), 
                    max.overlaps = Inf, size = 4.5, box.padding = 0.6, point.padding = 0.2,force = 1, show.legend = FALSE) +
    scale_color_manual(values = source_colors) +
    theme_minimal(base_size = 13) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        legend.position = "none"
    ) +
    labs(x = NULL, y = expression(-log[10](adjusted~p~value)))
    return(p)
}

p = plot_enrichment_dotplot(as.data.frame(df),"05_med_no_same_molGene_noMHC_gprofiler.pdf")
ggsave("05_med_no_same_molGene_noMHC_gprofiler.pdf", p, width = 6, height = 5)