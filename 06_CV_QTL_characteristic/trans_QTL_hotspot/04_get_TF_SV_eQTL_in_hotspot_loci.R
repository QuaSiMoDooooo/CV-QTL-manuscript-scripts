#!/usr/bin/env Rscript

# -------------
# FileName     : 04_get_TF_SV_eQTL_in_hotspot_loci.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Identify transcription factor SV-eQTLs located in trans-QTL hotspot regions
# -------------

library(vroom)
library(dplyr)
library(stringr)

qtl_file <- "~/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA/SV-expression/QTL_results/cis.filtered.txt.gz"
tf_file  <- "/home/wtian/data/TF/humanTFdb/TFs_Ensembl_v_1.01.txt"
hotspot_dir <- "~/project/HZAU_cohort_meth/wdy_assist_rerun/01_quick_stats_for_paper/trans_hotspot_loci/trans_hotspot_locis"

tfs <- vroom(tf_file, col_names = FALSE, delim = " ") %>% pull(1)

qtl <- vroom(qtl_file) %>% as.data.frame()
qtl$Phenotype = str_split(qtl$Phenotype, "\\.", simplify = TRUE)[,1]
qtl = qtl[qtl$Phenotype %in% tfs,]

out_file <- file.path("TF_SV_eQTL.tsv")
write.table(qtl, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

coords <- str_split_fixed(qtl$Variant_position, "[:-]", 3)
qtl$chr <- str_remove(coords[,1], "chr")
qtl$start <- as.integer(coords[,2])
qtl$end   <- as.integer(coords[,3])

hotspot_files <- list(
  APA         = file.path(hotspot_dir, "trans-apaQTL.hotspot_area_merged.tsv"),
  expression  = file.path(hotspot_dir, "trans-eQTL.hotspot_area_merged.tsv"),
  methylation = file.path(hotspot_dir, "trans-meQTL.hotspot_area_merged.tsv"),
  splicing    = file.path(hotspot_dir, "trans-sQTL.hotspot_area_merged.tsv")
)

res_list <- list()

for (phe in names(hotspot_files)) {
  hs <- vroom(hotspot_files[[phe]], col_names = c("chr","start","end"), delim = "\t") %>% as.data.frame()
  hs$chr <- str_remove(hs$chr, "chr")
  
  for (i in seq_len(nrow(hs))) {
    h_chr   <- hs$chr[i]
    h_start <- hs$start[i]
    h_end   <- hs$end[i]
    
    matches <- qtl[qtl$chr == h_chr & qtl$start >= h_start & qtl$start <= h_end, ]
    
    if (nrow(matches) > 0) {
      matches$Hotspot_type <- phe
      matches$Hotspot_chr  <- h_chr
      matches$Hotspot_start <- h_start
      matches$Hotspot_end   <- h_end
      res_list[[length(res_list)+1]] <- matches
    }
  }
}

final_res <- do.call(rbind, res_list)

out_file <- file.path("TF_SV_eQTL_in_trans_hotspots.tsv")
write.table(final_res, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("Done! Results saved to: ", out_file)