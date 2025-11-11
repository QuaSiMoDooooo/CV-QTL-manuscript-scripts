# -------------
# FileName     : impute_biochemistry
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-20 19:19
# Last Modified: 2024-10-14 15:51
# Modified By  : EastsunW
# -------------
# Description  : 使用mice包进行缺失值填充，假设输入的数据框每行一个指标，每列一个样本，第一列和第一行作为行名和列名。
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mice))

# 插补中，需要列是指标，行是样本
impute_data <- function(df, method = "pmm") {
    col_names <- colnames(df)
    colnames(df) <- paste0("S", seq_len(ncol(df)))
    invisible({
        imputed_result <- mice(
            df,
            m = 5,
            maxit = 10,
            method = method,
            seed = 2024,
            printFlag = FALSE
        )
    })
    imputed_df <- complete(imputed_result)
    colnames(imputed_df) <- col_names
    return(imputed_df)
}

data_raw <- fread(
    snakemake@input[[1]],
    sep = "\t",
    na.strings = c("", "NA")
)

data_imputed <- data_raw %>%
    column_to_rownames(colnames(.)[1]) %>%
    t() %>%
    as.data.frame() %>%
    impute_data(method = snakemake@params[["method"]]) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(colnames(data_raw)[1])

fwrite(
    data_imputed,
    snakemake@output[[1]],
    sep = "\t",
    na = "NA",
    quote = FALSE
)
