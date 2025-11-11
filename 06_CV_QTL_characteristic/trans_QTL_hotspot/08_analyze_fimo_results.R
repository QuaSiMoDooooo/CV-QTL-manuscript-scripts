#!/usr/bin/env Rscript

# -------------
# FileName     : 08_analyze_fimo_results.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Analyze FIMO motif scanning results for trans-QTL hotspot regions
# -------------

library(data.table)
library(dplyr)

fimo_file <- "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/01_quick_stats_for_paper/trans_hotspot_loci/fimo_out/fimo.tsv"

fimo_data <- fread(fimo_file, header = TRUE)

p_summary <- fimo_data %>%
  summarise(
    min_p = min(`p-value`, na.rm = TRUE),
    max_p = max(`p-value`, na.rm = TRUE),
    median_p = median(`p-value`, na.rm = TRUE),
    mean_p = mean(`p-value`, na.rm = TRUE),
    q25 = quantile(`p-value`, 0.25, na.rm = TRUE),
    q75 = quantile(`p-value`, 0.75, na.rm = TRUE)
  )

print(p_summary)

unique_sequences_q05 <- length(unique(fimo_data$sequence_name[fimo_data$`q-value` < 0.05]))
percentage_q05 <- unique_sequences_q05 / length(unique(fimo_data$sequence_name)) * 100

cat("Sequences with q-value < 0.05:", unique_sequences_q05, "(", round(percentage_q05, 2), "%)\n")