#! Rscript
# -------------
# FileName     : plof_GO_enrich
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-03 23:55
# Last Modified: 2025-03-04 10:28
# Modified By  : EastsunW
# -------------
# Description  : 对被全SV覆盖的基因进行富集分析
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(clusterProfiler))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
setwd("/home/wangdy/Projects/Weibin/Downstreams/variant_annotation")

# 导入所有被SV覆盖的基因信息
gene_info <- fread("results/epi_enrichment/plof/whole_gene_results.txt")

gene_id <- bitr(
    gene_info$Gene_name,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)

go_enrich <- enrichGO(
    gene = gene_info$Gene_name,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
)
pdf("results/epi_enrichment/plof/GO_enrichment_dot.pdf", width = 5, height = 4)
dotplot(
    go_enrich,
    showCategory = 10
) + wdy_theme(10, 1)
dev.off()
