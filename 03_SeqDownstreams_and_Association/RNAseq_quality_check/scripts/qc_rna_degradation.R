# -------------
# FileName     : qc_rna_degradation
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-05-29 21:31
# Last Modified: 2024-06-06 16:10
# Modified By  : EastsunW
# -------------
# Description  : RNA降解检查的热图绘制和表达量检查，同时把覆盖度进行标准化，输出pdf和样本聚类信息
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")

script_root = "/home/wangdy/Projects/HZAU-Weibin/Downstreams/RNA_degradation_check"

sample_group_info_path = paste0(script_root, "/data/sample_group.txt")
sample_genebody_coverage_path = paste0(script_root, "/data/RNAseq_20240419.geneBodyCoverage.txt")
housekeeping_gene_path = paste0(script_root, "/data/hg38.HouseKeepingGenes.bed")
expression_log2tpm_path = paste0(script_root, "/data/RNAseq_20240419.log2tpm.txt")
refseq_id_table_path = paste0(script_root, "/data/refseq_name.txt")

get_cluster_info = function(order_list, raw_names) {
    level_list = lapply(
        order_list,
        function(item) {
            return(raw_names[item])
        }
    )
    level_df = data.frame(
        sample = NULL,
        cluster = NULL
    )
    for(item in names(level_list)) {
        temp_df = data.frame(
            sample = level_list[[item]],
            cluster = rep(item, length(level_list[[item]]))
        )
        level_df <- rbind(
            level_df,
            temp_df
        )
    }
    return(level_df)
}

# 导入数据

## 导入分组和批次信息
sample_group <- fread(
    sample_group_info_path,
    colClasses = "character",
    header = TRUE
) %>%
    mutate(
        quality = factor(quality, levels = c("A", "B", "C")),
        batch = factor(batch, levels = c("0119", "0429", "0531")),
    )
## 导入coverage数据
genebody_coverage <- fread(
    sample_genebody_coverage_path,
    header = TRUE,
    col.names = c("sample", paste0("R",1:100))
) %>%
    column_to_rownames("sample") %>%
    t() %>%
    as.data.frame() %>%
    mutate_all(~scale(., center = FALSE, scale = TRUE)) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    arrange(factor(sample, levels = sample_group$sample)) %>%
    column_to_rownames("sample") %>%
    as.matrix()

## 导入refseq转换表
refseq_id_table <- fread(refseq_id_table_path)
## 导入管家基因的列表（bed）
housekeeping_gene <- fread(
    housekeeping_gene_path,
    header = FALSE
) %>%
    select(4) %>%
    rename(refseq=1) %>%
    mutate(geneID = refseq_id_table$`Gene stable ID`[match(
        refseq,
        refseq_id_table$`RefSeq mRNA ID`)
    ]) %>%
    filter(!is.na(geneID))

## 导入表达量信息并筛选出housekeeping
expression_matrix <- fread(
    expression_log2tpm_path,
    header = TRUE
) %>%
    select(all_of(c("gene_ID", "gene_symbol", sample_group$sample))) %>%
    filter(gene_ID %in% housekeeping_gene$geneID) %>%
    mutate(
        mean = apply(., 1, function(row) {
            return(mean(as.numeric(row[-c(1:2)])))
        }),
        median = apply(., 1, function(row) {
            return(median(as.numeric(row[-c(1:2)])))
        }),
        CV = apply(., 1, function(row) {
            sd = sd(as.numeric(row[-c(1:2)]))
            mean = mean(as.numeric(row[-c(1:2)]))
            return(sd / mean)
        })
    ) %>%
    filter(mean > 1, median > 1, CV < 0.2) %>%
    select(-mean, -median, -CV)

# 可视化基因覆盖度（热图）
heatmap_no_group = draw(Heatmap(
    genebody_coverage,
    name = "Coverage",
    col = circlize::colorRamp2(
        seq(0, 1.2, 0.2),
        c("white", "#FCDE83", "#FB9D00", "#FF3A00", "red", "#AA1500", "black")
    ),
    rect_gp = gpar(col = NA, lwd = 0),
    border = "black",
    row_title = "Samples and their quality level",
    row_title_side = "left",
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    row_km = 3,
    row_names_gp = gpar(fontsize = 6),
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_title = "Gene body percentile (5' to 3')",
    column_title_side = "bottom"
))
## 获取聚类信息，加到groupinfo后面
sample_info_with_cluster = get_cluster_info(
    row_order(heatmap_no_group),
    rownames(genebody_coverage)
) %>%
    inner_join(
        y = sample_group,
        by = "sample"
    ) %>%
    arrange(factor(sample, levels = rownames(genebody_coverage)))

group_col = c("A" = "orange", "B" = "orangered", "C" = "orangered4")
batch_col = c("0119" = "orchid", "0429" = "orchid4", "0531" = "maroon2")
cluster_col = c("1" = "royalblue", "2" = "skyblue", "3" = "skyblue4")
sample_annotation_cluster = rowAnnotation(
    group = sample_info_with_cluster$quality,
    batch = sample_info_with_cluster$batch,
    cluster = sample_info_with_cluster$cluster,
    col = list(
        group = group_col,
        batch = batch_col,
        cluster = cluster_col
    )
)

## 按照质量分组可视化
heatmap_by_group = draw(Heatmap(
    genebody_coverage,
    name = "Coverage",
    col = circlize::colorRamp2(
        seq(0, 1.2, 0.2),
        c("white", "#FCDE83", "#FB9D00", "#FF3A00", "red", "#AA1500", "black")
    ),
    rect_gp = gpar(col = NA, lwd = 0),
    border = "black",
    row_title = "Samples and their group info",
    row_title_side = "left",
    left_annotation = sample_annotation_cluster,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    row_names_gp = gpar(fontsize = 8),
    row_split = sample_info_with_cluster$quality,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_title = "Gene body percentile (5' to 3')",
    column_title_side = "bottom"
))


## 按照cluster可视化
heatmap_by_cluster = draw(Heatmap(
    genebody_coverage,
    name = "Coverage",
    col = circlize::colorRamp2(
        seq(0, 1.2, 0.2),
        c("white", "#FCDE83", "#FB9D00", "#FF3A00", "red", "#AA1500", "black")
    ),
    rect_gp = gpar(col = NA, lwd = 0),
    border = "black",
    row_title = "Samples and their group info",
    row_title_side = "left",
    left_annotation = sample_annotation_cluster,
    cluster_rows = FALSE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    row_order = get_cluster_info(
        row_order(heatmap_no_group),
        rownames(genebody_coverage)
    )$sample,
    row_names_gp = gpar(fontsize = 8),
    row_split = sample_info_with_cluster$cluster,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_title = "Gene body percentile (5' to 3')",
    column_title_side = "bottom"
))

## 按照batch可视化
heatmap_by_batch = draw(Heatmap(
    genebody_coverage,
    name = "Coverage",
    col = circlize::colorRamp2(
        seq(0, 1.2, 0.2),
        c("white", "#FCDE83", "#FB9D00", "#FF3A00", "red", "#AA1500", "black")
    ),
    rect_gp = gpar(col = NA, lwd = 0),
    border = "black",
    row_title = "Samples and their group info",
    row_title_side = "left",
    left_annotation = sample_annotation_cluster,
    cluster_rows = TRUE,
    show_row_names = FALSE,
    show_row_dend = FALSE,
    row_order = rownames(genebody_coverage),
    row_names_gp = gpar(fontsize = 8),
    row_split = sample_info_with_cluster$batch,
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_title = "Gene body percentile (5' to 3')",
    column_title_side = "bottom"
))

# 管家基因的表达量差异测试
gene_expression_long = expression_matrix %>%
    as.data.frame() %>%
    select(-gene_symbol) %>%
    pivot_longer(-gene_ID, names_to = "sample", values_to = "expression") %>%
    left_join(
        y = sample_info_with_cluster,
        by = "sample"
    )
## 按照样本分组
exp_group <- gene_expression_long %>%
    select(sample, quality, expression) %>%
    group_by(sample, quality) %>%
    mutate(quality_mean = mean(expression)) %>%
    distinct(sample, quality, quality_mean) %>%
    ggplot(aes(
        x = quality, y = quality_mean, fill = quality
    )) +
    geom_boxplot(width = 0.6) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("A", "B"),
            c("B", "C"),
            c("A", "C")
        ),
        label = "p.signif"
    ) +
    ggsci::scale_fill_aaas() +
    labs(
        x = "Quality",
        y = "Average expression"
    ) +
    wdy_theme(base_size = 10, legend.position = "none")
## 按照聚类
exp_cluster <- gene_expression_long %>%
    select(sample, cluster, expression) %>%
    group_by(sample, cluster) %>%
    mutate(cluster_mean = mean(expression)) %>%
    distinct(sample, cluster, cluster_mean) %>%
    ggplot(aes(
        x = cluster, y = cluster_mean, fill = cluster
    )) +
    geom_boxplot(width = 0.6) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("1", "2"),
            c("2", "3"),
            c("1", "3")
        ),
        label = "p.signif"
    ) +
    ggsci::scale_fill_aaas() +
    labs(
        x = "Cluster",
        y = "Average expression"
    ) +
    wdy_theme(base_size = 10, legend.position = "none")
## 按照批次
exp_batch <- gene_expression_long %>%
    select(sample, batch, expression) %>%
    group_by(sample, batch) %>%
    mutate(batch_mean = mean(expression)) %>%
    distinct(sample, batch, batch_mean) %>%
    ggplot(aes(
        x = batch, y = batch_mean, fill = batch
    )) +
    geom_boxplot(width = 0.6) +
    ggpubr::stat_compare_means(
        comparisons = list(
            c("0119", "0429"),
            c("0429", "0531"),
            c("0119", "0531")
        ),
        label = "p.signif"
    ) +
    ggsci::scale_fill_aaas() +
    labs(
        x = "Batch",
        y = "Average expression"
    ) +
    wdy_theme(base_size = 10, legend.position = "none")

pdf(paste0(script_root, "/results/geneBody_coverage_heatmaps.pdf"), width = 4.5, height = 4.5, onefile = TRUE)
heatmap_by_group
heatmap_by_cluster
heatmap_by_batch
dev.off()

pdf(paste0(script_root, "/results/gene_expression.pdf"), width = 3, height = 3, onefile = TRUE)
exp_group
exp_cluster
exp_batch
dev.off()

genebody_coverage %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    fwrite(paste0(script_root, "/results/genebody_coverage.normalized.txt"), sep = "\t", quote = FALSE)

fwrite(
    sample_info_with_cluster,
    paste0(script_root, "/results/sample_group_with_cluster.txt"),
    sep = "\t",
    quote = FALSE
)
