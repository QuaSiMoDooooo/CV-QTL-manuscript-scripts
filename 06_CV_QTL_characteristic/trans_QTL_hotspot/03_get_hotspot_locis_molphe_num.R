#!/usr/bin/env Rscript

# -------------
# FileName     : 03_get_hotspot_locis_molphe_num.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Count unique molecular phenotypes in trans-QTL hotspot regions across variant types
# -------------

library(vroom)
library(dplyr)
library(stringr)

qtl_dir <- "~/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA"
hotspot_dir <- "~/project/HZAU_cohort_meth/wdy_assist_rerun/01_quick_stats_for_paper/trans_hotspot_loci/trans_hotspot_locis"

phenotypes <- c("APA", "expression", "methylation", "splicing")
variants <- c("SNP", "MNV", "InDel", "SV")

results <- list()

for (phe in phenotypes) {
    hotspot_file <- file.path(
        hotspot_dir,
        paste0(
            "trans-", 
            ifelse(phe=="expression", "eQTL",
                   ifelse(phe=="methylation","meQTL",
                          ifelse(phe=="splicing","sQTL","apaQTL"))),
            ".hotspot_area_merged.tsv"
        )
    )
    hotspots <- vroom(hotspot_file, col_names = c("chr","start","end"), delim = "\t") %>% 
        as.data.frame()
    hotspots$chr <- str_remove(hotspots$chr, "chr")
    
    phe_phenos <- c()
    
    for (var in variants) {
        qtl_file <- file.path(qtl_dir, paste0(var, "-", phe), "QTL_results", "trans.filtered.txt.gz")
        if (!file.exists(qtl_file)) next
        
        message("Processing: ", var, "-", phe)
        
        qtl <- vroom(qtl_file, col_select = c("Variant_position","Phenotype")) %>% 
            as.data.frame()
        
        coords <- str_split_fixed(qtl$Variant_position, "[:-]", 3)
        qtl$chr <- str_remove(coords[,1], "chr")
        qtl$pos <- as.integer(coords[,2])
        
        for (i in seq_len(nrow(hotspots))) {
            h_chr <- hotspots$chr[i]
            h_start <- hotspots$start[i]
            h_end <- hotspots$end[i]
            
            in_region <- unique(qtl$Phenotype[qtl$chr == h_chr & qtl$pos >= h_start & qtl$pos <= h_end])
            phe_phenos <- c(phe_phenos, in_region)
        }
    }
    
    phe_phenos <- unique(phe_phenos)
    
    results[[length(results)+1]] <- data.frame(
        Phenotype_type = phe,
        Unique_Phenotype_N = length(phe_phenos),
        stringsAsFactors = FALSE
    )
}

final_res <- do.call(rbind, results)
out_file <- file.path("trans_hotspot_qtl_counts_merged.tsv")
write.table(final_res, out_file, sep = "\t", row.names = FALSE, quote = FALSE)

message("Done! Results saved to: ", out_file)