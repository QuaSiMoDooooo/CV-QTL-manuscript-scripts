# -------------
# FileName     : merge
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-11 00:52
# Last Modified: 2024-04-11 15:27
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

temp_df <- data.frame(
    chr = factor(character(length = 0), levels = c(paste0("chr", c(as.character(1:22), "X", "Y", "M")))),
    start = numeric(length = 0),
    end = numeric(length = 0)
)
for (file_path in unlist(snakemake@input)) {
    sample_name <- str_split(basename(file_path), "\\.", simplify = TRUE)[1]
    my_header <- c(
        "chr", "start", "end",
        sample_name, "hap",
        "reads_all", "reads_mod", "reads_unmod",
        "mod_frac"
    )
    temp_df <<- fread(file_path, col.names = my_header) %>%
        select(chr, start, end, !!sample_name) %>%
        full_join(
            y = temp_df,
            by = c("chr", "start", "end")
        ) %>%
        mutate_at("chr", ~factor(., levels = c(paste0("chr", c(as.character(1:22), "X", "Y", "M")))))
    temp = gc()
    rm(temp)
}
fwrite(temp_df, snakemake@output[[1]], sep = "\t", quote = FALSE, na = "NA")
