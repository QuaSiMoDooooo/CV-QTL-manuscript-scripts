#!/usr/bin/env Rscript

# -------------
# FileName     : 04_trait_or_gene_based_graph_vis.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Visualize colocalization results as network graphs for trait or gene-based analysis
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggraph))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(stringsAsFactors = FALSE)
mc <- getOption('mc.cores', 36)
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc")

gwas_qtl_coloc = fread("03_gwas_qtl_colocResultsig.tsv")
gwas_qtl_coloc$phe = gsub("\\.[0-9]+","",gwas_qtl_coloc$phe)
eqtl_qtl_coloc = fread("03_eqtl_qtl_colocResultsig.tsv")
colnames(gwas_qtl_coloc)[1] = "node"
colnames(eqtl_qtl_coloc)[1] = "node"

gwas_qtl_coloc$node1_type = "trait"
eqtl_qtl_coloc$node1_type = "expression"

all_coloc = rbind(gwas_qtl_coloc, eqtl_qtl_coloc)
colnames(all_coloc)[1] = "node1"
colnames(all_coloc)[4] = "node2"
all_coloc$node2_type = all_coloc$molphe_name
all_coloc = all_coloc %>% select(node1,node1_type,node2,node2_type,var_type)
write.table(all_coloc, file = "04_all_coloc_asso.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

all_coloc <- all_coloc %>%
  mutate(node1 = as.character(node1),
         node2 = as.character(node2))

edges <- all_coloc %>%
  group_by(node1, node2) %>%
  summarize(var_type_combined = paste(sort(unique(var_type)), collapse = ","), .groups = "drop") %>%
  mutate(var_type = case_when(
    var_type_combined == "SNP" ~ "Only SNV",
    var_type_combined == "SV" ~ "Only SV",
    var_type_combined == "SNP,SV" ~ "SV and SNV",
    TRUE ~ var_type_combined
  )) %>%
  select(node1, node2, var_type)
edges$var_type <- factor(edges$var_type, levels = c("Only SNV", "Only SV", "SV and SNV"))

nodes1 <- all_coloc %>%
  select(name = node1, node_type = node1_type)
nodes2 <- all_coloc %>%
  select(name = node2, node_type = node2_type)
nodes <- bind_rows(nodes1, nodes2) %>%
  distinct(name, .keep_all = TRUE)

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)

node_type_colors <- c(
  "trait" = "#1f77b4",
  "expression" = "#ff7f0e",
  "methylation" = "#B77554",
  "splicing" = "#d62728",
  "APA" = "#9467bd"
)

var_type_colors <- c(
  "Only SNV" = "#F37C72",
  "Only SV" = "#78A7C4",
  "SV and SNV" = "#cec0ff"
)

pdf("04_all_coloc_asso_component10.pdf", width = 12, height = 10)
components_info <- components(g)
V(g)$component_size <- components_info$csize[components_info$membership]
g_filtered <- induced_subgraph(g, vids = V(g)[component_size >= 10])
V(g_filtered)$is_trait <- ifelse(V(g_filtered)$node_type == "trait", TRUE, FALSE)
V(g_filtered)$is_node1 <- ifelse(V(g_filtered)$name %in% all_coloc$node1, TRUE, FALSE)

ggraph(g_filtered, layout = 'kk') + 
  geom_edge_link(aes(color = var_type), alpha = 0.8) +
  geom_node_point(aes(color = node_type), size = 2) +
  geom_node_text(aes(
      label = ifelse(is_node1, name, ""), 
      fontface = ifelse(is_trait, "bold", "plain")
    ),
    color = ifelse(V(g_filtered)$is_trait, "red", "black"), 
    repel = TRUE, size = 2, max.overlaps = 10000, show.legend = FALSE
  ) +
  scale_color_manual(values = node_type_colors) +
  scale_edge_color_manual(values = var_type_colors) +
  theme_void() +
  theme(legend.position = "right",
  legend.text = element_text(size = 12),
  legend.title = element_text(size = 14)
  )
dev.off()

components_info <- components(g)
largest_comp_id <- which.max(components_info$csize)
g_largest <- induced_subgraph(g, vids = V(g)[components_info$membership == largest_comp_id])
V(g_largest)$is_trait <- ifelse(V(g_largest)$node_type == "trait", TRUE, FALSE)
V(g_largest)$is_node1 <- ifelse(V(g_largest)$name %in% all_coloc$node1, TRUE, FALSE)
pdf("04_all_coloc_asso_component_the_most.pdf", width = 12, height = 10)
ggraph(g_largest, layout = 'kk') + 
  geom_edge_link(aes(color = var_type), alpha = 0.8) +
  geom_node_point(aes(color = node_type), size = 5) +
  geom_node_text(aes(
      label = ifelse(is_node1, name, ""), 
      fontface = ifelse(is_trait, "bold", "plain")
    ),
    color = ifelse(V(g_largest)$is_trait, "red", "black"), 
    repel = TRUE, size = 5, max.overlaps = 10000, show.legend = FALSE
  ) +
  scale_color_manual(values = node_type_colors) +
  scale_edge_color_manual(values = var_type_colors) +
  theme_void() +
  theme(legend.position = "right",
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16)
  )
dev.off()

components_info <- components(g)
is_expression_node <- V(g)$node_type == "expression"
expression_components <- unique(components_info$membership[is_expression_node])
expression_comp_sizes <- components_info$csize[expression_components]
largest_expression_comp_id <- expression_components[which.max(expression_comp_sizes)]
g_largest_expression <- induced_subgraph(g, vids = V(g)[components_info$membership == largest_expression_comp_id])
V(g_largest_expression)$is_expression <- ifelse(V(g_largest_expression)$node_type == "expression", TRUE, FALSE)
V(g_largest_expression)$is_node1 <- ifelse(V(g_largest_expression)$name %in% all_coloc$node1, TRUE, FALSE)
pdf("04_all_coloc_asso_component_the_most_expression.pdf", width = 12, height = 10)
ggraph(g_largest_expression, layout = 'kk') + 
  geom_edge_link(aes(color = var_type), alpha = 0.8) +
  geom_node_point(aes(color = node_type), size = 4) +
  geom_node_text(aes(
      label = ifelse(is_node1, name, ""), 
      fontface = ifelse(is_expression, "bold", "plain")
    ),
    color = ifelse(V(g_largest_expression)$is_expression, "red", "black"), 
    repel = TRUE, size = 5, max.overlaps = 10000, show.legend = FALSE
  ) +
  scale_color_manual(values = node_type_colors) +
  scale_edge_color_manual(values = var_type_colors) +
  theme_void() +
  theme(legend.position = "right",
  legend.text = element_text(size = 14),
  legend.title = element_text(size = 16)
  )
dev.off()

if (!dir.exists("04_trait_or_gene_based_graph")) { dir.create("04_trait_or_gene_based_graph")}
components_info <- components(g)
target_nodes <- V(g)[node_type %in% c("trait", "expression")]
for (target_node in target_nodes) {
  node_name <- target_nodes$name[target_node]
  comp_id <- components_info$membership[target_node]
  g_sub <- induced_subgraph(g, vids = V(g)[components_info$membership == comp_id])
  V(g_sub)$is_trait <- ifelse(V(g_sub)$node_type == "trait", TRUE, FALSE)
  V(g_sub)$is_node1 <- ifelse(V(g_sub)$name %in% all_coloc$node1, TRUE, FALSE)
  output_file <- file.path("04_trait_or_gene_based_graph", paste0(node_name, ".pdf"))
  pdf(output_file, width = 12, height = 10, family = "Helvetica")
  p = ggraph(g_sub, layout = 'kk') + 
    geom_edge_link(aes(color = var_type), alpha = 0.8) +
    geom_node_point(aes(color = node_type), size = 5) +
    geom_node_text(aes(
        label = name,
        fontface = ifelse(is_trait, "bold", "plain")
      ),
      color = ifelse(V(g_sub)$is_trait, "red", "black"), 
      repel = TRUE, size = 5, max.overlaps = 10000, show.legend = FALSE
    ) +
    scale_color_manual(values = node_type_colors) +
    scale_edge_color_manual(values = var_type_colors) +
    theme_void() +
    theme(legend.position = "right",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
    ) 
  print(p)
  dev.off()
}