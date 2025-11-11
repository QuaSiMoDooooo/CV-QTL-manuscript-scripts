# -------------
# FileName     : merge_count
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-02-25 17:18
# Last Modified: 2024-11-14 20:13
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

merged_count <- data.frame()
merged_tpm <- data.frame()

test <- suppressMessages({
    suppressWarnings({
        sapply(
            X = unlist(snakemake@input[["quantity"]]),
            FUN = function(file) {
                sample_name <- str_split(basename(file), "\\.", simplify = TRUE)[1]
                temp_exp <- fread(file)
                if (nrow(merged_count) == 0) {
                    merged_count <<- temp_exp %>%
                        select(Name, NumReads) %>%
                        column_to_rownames("Name") %>%
                        setnames(sample_name)
                } else {
                    merged_count <<- cbind(
                        merged_count,
                        temp_exp %>%
                            select(Name, NumReads) %>%
                            column_to_rownames("Name") %>%
                            setnames(sample_name)
                    )
                }
                if (nrow(merged_tpm) == 0) {
                    merged_tpm <<- temp_exp %>%
                        select(Name, TPM) %>%
                        column_to_rownames("Name") %>%
                        setnames(sample_name)
                } else {
                    merged_tpm <<- cbind(
                        merged_tpm,
                        temp_exp %>%
                            select(Name, TPM) %>%
                            column_to_rownames("Name") %>%
                            setnames(sample_name)
                    )
                }
            }
        )
    })
})

# id_mapping <- fread(snakemake@input[["idmapping"]])

merged_count_name <- merged_count %>%
    rownames_to_column("gene_ID") %>%
    # inner_join(
    #     y = id_mapping,
    #     by = "gene_ID",
    #     na_matches = "never",
    #     relationship = "many-to-many"
    # ) %>%
    # select(-gene_ID) %>%
    # group_by(gene_name) %>%
    # summarise_if(is.numeric, sum) %>%
    rename("ID" = "gene_ID") %>%
    filter(!is.na(ID), ID != "NA", ID != "")

merged_tpm_name <- merged_tpm %>%
    rownames_to_column("gene_ID") %>%
    # inner_join(
    #     y = id_mapping,
    #     by = "gene_ID",
    #     na_matches = "never",
    #     relationship = "many-to-many"
    # ) %>%
    # select(-gene_ID) %>%
    # group_by(gene_name) %>%
    # summarise_if(is.numeric, sum) %>%
    rename("ID" = "gene_ID") %>%
    filter(!is.na(ID), ID != "NA", ID != "")

log2tpm <- merged_tpm_name %>%
    column_to_rownames("ID") %>%
    mutate_all(~ log2(. + 1)) %>%
    rownames_to_column("ID")

fwrite(
    merged_count_name,
    snakemake@output[["rawcount"]],
    quote = FALSE,
    sep = "\t"
)
fwrite(
    merged_tpm_name,
    snakemake@output[["rawtpm"]],
    quote = FALSE,
    sep = "\t"
)
fwrite(
    log2tpm,
    snakemake@output[["log2tpm"]],
    quote = FALSE,
    sep = "\t"
)
