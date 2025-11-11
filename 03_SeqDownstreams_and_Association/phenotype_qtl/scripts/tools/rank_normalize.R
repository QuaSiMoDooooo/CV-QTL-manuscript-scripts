# -------------
# FileName     : rank_normalize
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-19 20:48
# Last Modified: 2024-09-20 22:21
# Modified By  : EastsunW
# -------------
# Description  : 接受表型数据，对表型数据进行分位数归一化，输入的数据框每行一个特征，每列一个样本
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

rank_phenotype <- function(df) {
    df_makers <- df[[1]]
    df_headers <- colnames(df)
    df_values <- df[, -1]
    expr_rank <- t(
        apply(
            X = as.matrix(df_values),
            MARGIN = 1,
            FUN = function(x) {
                qnorm(rank(x, ties.method = "average") / (length(x) + 1))
            }
        )
    )
    expr_rank <- round(expr_rank, 4)
    expr_rank <- cbind(df_makers, expr_rank)
    expr_rank <- as.data.frame(expr_rank, stringsAsFactors = FALSE)
    colnames(expr_rank) <- df_headers
    return(expr_rank)
}

# 导入原始的定量信息
input_phenotype <- fread(snakemake@input[[1]])

ranked_phenotype <- rank_phenotype(input_phenotype)

# 输出rank后的定量信息
fwrite(
    x = ranked_phenotype,
    file = snakemake@output[[1]],
    sep = "\t",
    quote = FALSE
)
