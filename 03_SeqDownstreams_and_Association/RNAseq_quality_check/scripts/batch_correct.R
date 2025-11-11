# -------------
# FileName     : batch_correct
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-05-30 21:55
# Last Modified: 2024-06-02 20:30
# Modified By  : EastsunW
# -------------
# Description  : 使用不同的批次变量分别矫正RNAseq表达矩阵的批次效应，比较矫正的结果，输出各种矫正的PCA图和校正后的结果
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(factoextra))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")


if (.Platform$OS.type == "windows") {
    script_root = "D:/OneDrive/投稿管理/未病人群队列分析/批次效应"
} else {
    script_root = "/home/wangdy/Projects/HZAU-Weibin/Downstreams/expression_batch_correct"
}
setwd(script_root)

sample_group_info_path = paste0(script_root, "/data/sample_group_with_cluster.txt")
expression_log2tpm_path = paste0(script_root, "/data/RNAseq_20240419.log2tpm.txt")


# 读取表达矩阵数据
expression_log2tpm_raw <- fread(expression_log2tpm_path)
expression_gene_list <- expression_log2tpm_raw[, c(1,2)]
expression_log2tpm_before <- expression_log2tpm_raw %>%
    column_to_rownames("gene_ID") %>%
    select(-gene_symbol)
# 读取分组信息
sample_group <- fread(sample_group_info_path, colClasses="character") %>%
    column_to_rownames("sample")

#===============================================
known_batch_remove1 <- function(var, rawdata, group_data) {
    model0 <- model.matrix(~1, data = group_data)
    data_corrected <- ComBat(
        dat = rawdata,
        batch = sample_group[[var]],
        mod = model0
    )
    pca_before <- PCA(t(rawdata), graph = FALSE)
    pca_ind_plot_before1 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before ", var, " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["cluster"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = paste0(var, collapse = " and "),
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_before2 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before ", var, " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["batch"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = paste0(var, collapse = " and "),
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_before3 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before ", var, " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["quality"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = paste0(var, collapse = " and "),
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_after <- PCA(t(data_corrected), graph = FALSE)
    pca_ind_plot_after1 <- fviz_pca_ind(
        pca_after,
        title = paste0("After ", var, " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["cluster"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = paste0(var, collapse = " and "),
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_after2 <- fviz_pca_ind(
        pca_after,
        title = paste0("After ", var, " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["batch"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = paste0(var, collapse = " and "),
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_after3 <- fviz_pca_ind(
        pca_after,
        title = paste0("After ", var, " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["quality"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = paste0(var, collapse = " and "),
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )

    pca_compare_plot <- ggarrange(
        plotlist = list(
            ggarrange(
                pca_ind_plot_before1, pca_ind_plot_after1,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            ),
            ggarrange(
                pca_ind_plot_before2, pca_ind_plot_after2,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            ),
            ggarrange(
                pca_ind_plot_before3, pca_ind_plot_after3,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            )
        ),
        nrow = 1, ncol = 3
    ) + wdy_theme(10)
    return(list(
        corrected_data = data_corrected,
        plot = pca_compare_plot
    ))
}

pca_cluster = known_batch_remove1(
    var = "cluster",
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

pca_batch = known_batch_remove1(
    var = "batch",
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

pca_quality = known_batch_remove1(
    var = "quality",
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

fwrite(
    cbind(expression_gene_list, pca_quality$corrected_data),
    "results/quality_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
fwrite(
    cbind(expression_gene_list, pca_batch$corrected_data),
    "results/batch_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
fwrite(
    cbind(expression_gene_list, pca_cluster$corrected_data),
    "results/cluster_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
pdf("results/PCA_single_correction.pdf", width = 10, height = 6, onefile = TRUE)
pca_cluster$plot
pca_batch$plot
pca_quality$plot
dev.off()


#===============================================
known_batch_remove2 <- function(var, rawdata, group_data) {
    var_temp = var
    model0 <- model.matrix(~1, data = group_data)
    while(length(var) > 0) {
        data_corrected <- ComBat(
            dat = rawdata,
            batch = sample_group[[var[1]]],
            mod = model0
        )
        var = var[-1]
    }
    pca_before <- PCA(t(rawdata), graph = FALSE)
    pca_ind_plot_before1 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before ", paste0(var_temp, collapse = "/"), " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["cluster"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "cluster",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_before2 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before ", paste0(var_temp, collapse = "/"), " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["batch"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "batch",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_before3 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before ", paste0(var_temp, collapse = "/"), " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["quality"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "quality",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_after <- PCA(t(data_corrected), graph = FALSE)
    pca_ind_plot_after1 <- fviz_pca_ind(
        pca_after,
        title = paste0("After ", paste0(var_temp, collapse = "/"), " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["cluster"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "cluster",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_after2 <- fviz_pca_ind(
        pca_after,
        title = paste0("After ", paste0(var_temp, collapse = "/"), " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["batch"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "batch",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_after3 <- fviz_pca_ind(
        pca_after,
        title = paste0("After ", paste0(var_temp, collapse = "/"), " corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["quality"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "quality",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_compare_plot <- ggarrange(
        plotlist = list(
            ggarrange(
                pca_ind_plot_before1, pca_ind_plot_after1,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            ),
            ggarrange(
                pca_ind_plot_before2, pca_ind_plot_after2,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            ),
            ggarrange(
                pca_ind_plot_before3, pca_ind_plot_after3,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            )
        ),
        nrow = 1, ncol = 3
    ) + wdy_theme(10)
    return(list(
        corrected_data = data_corrected,
        plot = pca_compare_plot
    ))
}

pca_cluster_batch = known_batch_remove2(
    var = c("cluster", "batch"),
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

pca_cluster_quality = known_batch_remove2(
    var = c("cluster", "quality"),
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

pca_batch_quality = known_batch_remove2(
    var = c("batch", "quality"),
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

fwrite(
    cbind(expression_gene_list, pca_cluster_batch$corrected_data),
    "results/cluster_batch_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
fwrite(
    cbind(expression_gene_list, pca_batch_quality$corrected_data),
    "results/batch_quality_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
fwrite(
    cbind(expression_gene_list, pca_cluster_quality$corrected_data),
    "results/cluster_quality_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
pdf("results/PCA_double_correction.pdf", width = 10, height = 6, onefile = TRUE)
pca_cluster_batch$plot
pca_cluster_quality$plot
pca_batch_quality$plot
dev.off()

#=============================================
known_batch_remove3 <- function(var, rawdata, group_data) {
    model0 <- model.matrix(~1, data = group_data)
    while(length(var) > 0) {
        data_corrected <- ComBat(
            dat = rawdata,
            batch = sample_group[[var[1]]],
            mod = model0
        )
        var = var[-1]
    }
    pca_before <- PCA(t(rawdata), graph = FALSE)
    pca_ind_plot_before1 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before triple corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["cluster"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "cluster",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_before2 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before triple corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["batch"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "batch",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_before3 <- fviz_pca_ind(
        pca_before,
        title = paste0("Before triple corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["quality"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "quality",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_after <- PCA(t(data_corrected), graph = FALSE)
    pca_ind_plot_after1 <- fviz_pca_ind(
        pca_after,
        title = paste0("After triple corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["cluster"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "cluster",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_after2 <- fviz_pca_ind(
        pca_after,
        title = paste0("After triple corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["batch"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "batch",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_ind_plot_after3 <- fviz_pca_ind(
        pca_after,
        title = paste0("After triple corrected"),
        geom = c("point"),
        col.ind = factor(group_data[["quality"]]),
        pointshape = 19,
        pointsize = 1,
        alpha.ind = 0.8,
        addEllipses = TRUE,
        legend.title = "quality",
        palette = "aaas"
    ) + wdy_theme(
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank()
    )
    pca_compare_plot <- ggarrange(
        plotlist = list(
            ggarrange(
                pca_ind_plot_before1, pca_ind_plot_after1,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            ),
            ggarrange(
                pca_ind_plot_before2, pca_ind_plot_after2,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            ),
            ggarrange(
                pca_ind_plot_before3, pca_ind_plot_after3,
                ncol = 1,
                legend = "bottom",
                common.legend = TRUE
            )
        ),
        nrow = 1, ncol = 3
    ) + wdy_theme(10)
    return(list(
        corrected_data = data_corrected,
        plot = pca_compare_plot
    ))
}
pca_cluster_batch_quality = known_batch_remove3(
    var = c("cluster", "batch", "quality"),
    rawdata = expression_log2tpm_before,
    group_data = sample_group
)

fwrite(
    cbind(expression_gene_list, pca_cluster_batch_quality$corrected_data),
    "results/cluster_quality_batch_corrected.txt",
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
)
pdf("results/PCA_triple_correction.pdf", width = 10, height = 6, onefile = TRUE)
pca_cluster_batch_quality$plot
dev.off()

