#! Rscript
# -------------
# FileName     : plot_QTL_result
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-14 16:57
# Last Modified: 2024-12-14 21:00
# Modified By  : EastsunW
# -------------
# Description  : 输入QTLtype、要画的变异、要画的表型，输出一个QTL的箱线图
# -------------

suppressPackageStartupMessages(library(optparse))

get_options <- function() {
    usage <- "usage: Rscript %prog [options] <arguments>"
    description <- "description"
    option_list <- list(
        make_option(
            c("-q", "--qtl"),
            type    = "character",
            dest    = "qtl",
            metavar = "<SNP-eQTL|...>",
            default = FALSE,
            action  = "store",
            help    = "QTL type, variant-xQTL"
        ),
        make_option(
            c("-v", "--variant-id"),
            type    = "character",
            dest    = "variant",
            metavar = "<variant ID>",
            default = FALSE,
            action  = "store",
            help    = "Variant ID to plot, using ',' to split multiple IDs"
        ),
        make_option(
            c("-p", "--phenotype-id"),
            type    = "character",
            dest    = "phenotype",
            metavar = "<phenotype ID>",
            default = FALSE,
            action  = "store",
            help    = "Phenotype ID to plot, using ',' to split multiple IDs"
        ),
        make_option(
            c("-o", "--output"),
            type    = "character",
            dest    = "output",
            metavar = "<path>",
            default = "results/plot/QTL_plots",
            action  = "store",
            help    = "Output dir, the filename will be <dir>/<variant_type>_<qtl_type>_<variant_ID>_<phenotype_ID>_<cis|trans>.pdf"
        )
    )
    opt_obj <- OptionParser(
        usage           = usage,
        description     = description,
        add_help_option = TRUE,
        option_list     = option_list
    )
    opts <- optparse::parse_args(opt_obj)
    if (length(commandArgs(TRUE)) == 0) {
        optparse::print_help(opt_obj)
        if (interactive()) {
            stop("help requested")
        } else {
            quit(status = 0)
        }
    } else {
        return(opts)
    }
}

options <- get_options()
options[["variant_type"]] <- strsplit(options[["qtl"]], "-")[[1]][1]
options[["qtl_type"]] <- strsplit(options[["qtl"]], "-")[[1]][2]
options[["variant"]] <- strsplit(options[["variant"]], ",")[[1]]
options[["phenotype"]] <- strsplit(options[["phenotype"]], ",")[[1]]
if (any(unlist(options) == "")) {
    optparse::print_help(opt_obj)
    stop(
        paste0(
            "These argument are not specified:\n\t --",
            paste0(
                names(unlist(options)[unlist(options) == ""]),
                collapse = "\n\t --"
            )
        )
    )
}
if (!dir.exists(options[["output"]])) {
    dir.create(options[["output"]], recursive = TRUE)
}

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
source("https://gitee.com/eastsunw/personal_code_notebook/raw/master/plot_tools/wdy_theme.r")
custom_colors <- c(
    "#FA7F6F",
    "#8ECFC9",
    "#FFBE7A",
    "#82B0D2",
    "#BEB8DC",
    "#E7DAD2",
    "#999999"
)
setwd("/home/wangdy/Projects/Weibin/Downstreams/data_stat")

# 根据输入的QTL类型，读取对应的文件
qtl_result <- fread(
    file = paste0(
        "results/stat/QTL/",
        options[["variant_type"]], "_", options[["qtl_type"]],
        ".filtered.txt"
    ),
    sep = "\t"
) %>%
    filter(
        Variant %in% options[["variant"]],
        Phenotype %in% options[["phenotype"]] | Gene %in% options[["phenotype"]]
    )
genotype <- fread(
    file = paste0(
        "data/QTL_data/",
        options[["variant_type"]],
        "_",
        options[["qtl_type"]],
        ".genotype.txt"
    ),
    sep = "\t"
)
phenotype <- fread(
    file = paste0(
        "data/QTL_data/",
        options[["variant_type"]],
        "_",
        options[["qtl_type"]],
        ".phenotype.txt"
    ),
    sep = "\t"
)

plot_box_plot <- function(row) {
    # 处理数据来画图
    plot_data <- rbind(
        genotype %>% filter(ID == row[["Variant"]]),
        phenotype %>% filter(ID == row[["Phenotype"]])
    ) %>%
        mutate(ID = c("Genotype", "Quantity")) %>%
        column_to_rownames("ID") %>%
        t() %>%
        as.data.frame() %>%
        mutate_at(vars("Genotype"), ~ factor(., levels = as.character(0:2)))
    ggplot(
        data = plot_data,
        mapping = aes(
            x = Genotype,
            y = Quantity
        )
    ) +
        geom_boxplot(
            aes(
                group = Genotype,
                fill = Genotype
            ),
            color = "black"
        ) +
        geom_jitter(
            width = 0.1,
            height = 0,
            size = 1,
            alpha = 0.8,
            shape = 21,
            color = "gray20",
        ) +
        scale_fill_manual(
            values = c(
                "0" = custom_colors[2],
                "1" = custom_colors[3],
                "2" = custom_colors[1]
            )
        ) +
        labs(
            x = paste(
                row[["Variant"]],
                "genotype",
                sep = " "
            ),
            y = paste(
                "Normalized",
                ifelse(is.na(row[["Gene"]]), row[["Phenotype"]], row[["Gene"]]),
                "quantity",
                sep = " "
            )
        ) +
        # 在右上角添加QTL的样本量、P、Beta等信息
        annotate(
            "text",
            x = Inf,
            y = Inf,
            label = paste(
                "R = ", formatC(as.numeric(row[["R_pearson"]]), format = "f", digits = 2), "\n",
                "P = ", formatC(as.numeric(row[["P"]]), format = "e", digits = 2), "\n",
                "Beta = ", formatC(as.numeric(row[["Beta"]]), format = "f", digits = 2), "\n"
            ),
            hjust = 1.1,
            vjust = 1.1,
            size = 10 / .pt,
        ) +
        wdy_theme(
            legend.position = "none",
        )
}

plots <- apply(qtl_result, 1, plot_box_plot)
names(plots) <- apply(qtl_result, 1, function(row) {
    paste(
        paste(row[["qtl_type"]], options[["qtl_type"]], sep = "-"),
        row[["Variant"]],
        row[["Phenotype"]],
        sep = "_"
    )
})
lapply(
    names(plots),
    function(name) {
        ggsave(
            filename = paste0(options[["output"]], "/", name, ".pdf"),
            plot = plots[[name]],
            width = 3,
            height = 3
        )
    }
)
