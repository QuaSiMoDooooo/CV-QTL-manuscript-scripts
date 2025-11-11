#!/usr/bin/env Rscript

# -------------
# FileName     : 02_02_eqtl_qtl_each_file.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Perform colocalization analysis between eQTL and other QTL data for each molecular phenotype
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(xQTLbiolinks))
suppressPackageStartupMessages(library(parallel))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc")

mc <- getOption("mc.cores", 200)

range_dir = "../01_get_chr1-22_GWAS_and_eQTL_sentinelVARs_range/"

qtl_flt_dir = "../data/QTL_flt"
qtl_flt_files = list.files(qtl_flt_dir, pattern = "\\.tsv$", full.names = TRUE)
qtl_var_type_list = str_split(basename(qtl_flt_files), "_", simplify = TRUE)[, 3]
var_type_list = c("SNP","SV")
qtl_flt_files = qtl_flt_files[qtl_var_type_list %in% var_type_list]
qtl_var_type_list = qtl_var_type_list[qtl_var_type_list %in% var_type_list]

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/04_rerun_QTL/results"
molphe_names = c("splicing","APA","methylation")

overlap_phe_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc/02_01_overlap_phe"

filter_molphe_by_range <- function(df, range) {
    range_gr <- GRanges(
    seqnames = range$chr,
    ranges = IRanges(start = range$start, end = range$end)
    )
    qtl_gr <- GRanges(
    seqnames = df$chr,
    ranges = IRanges(start = df$start, end = df$end)
    )
    hits <- findOverlaps(qtl_gr, range_gr)
    df_overlap <- df[queryHits(hits), ]
    return(df_overlap)
}

for (j in 1:length(qtl_var_type_list)){
    var_type = qtl_var_type_list[j]
    print(var_type)

    freq_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_process_GWAS_bim2_and_bed/freq"
    freq_ref = fread(paste0(freq_dir,"/",var_type,".frq"))
    colnames(freq_ref) = c("rsid", "maf")
    setindex(freq_ref, rsid)

    range_subdir = paste0(range_dir,"eQTL_",var_type)
    eqtl_flt_file = qtl_flt_files[j]
    eqtl_flt = data.table(vroom(eqtl_flt_file, delim = "\t", show_col_types = FALSE, progress = FALSE))
    eqtl_flt = eqtl_flt %>% filter(molphe == "eQTL")
    genes = str_split(unique(eqtl_flt$Phenotype), "\\.", simplify = TRUE)[, 1]
    print(paste0("effective gene numbers:", length(genes)))
    rm(eqtl_flt)
    gc()

    eqtl_raw_dir = paste0(qtl_dir,"/",var_type,"-expression/QTL_results/cis_each")
    eqtl_raw_files = list.files(eqtl_raw_dir, pattern = "\\.txt$", full.names = TRUE)

    for (g in 1:length(genes)){
        gene = genes[g]
        colocResult_p_All_paths = paste0("02_02_eqtl_coloc_res/",var_type,"-",molphe_names,"/",gene,".tsv")
        if (all(file.exists(colocResult_p_All_paths))){
        print(paste0("gene ",gene," has been processed"))
        next
        } else {
        print(paste0("gene ",gene," is processing..."))
        }

        eqtl_raw_path = eqtl_raw_files[grepl(gene, basename(eqtl_raw_files))][1]
        if (is.na(eqtl_raw_path)){ 
            print(paste0("gene ",gene," not rerun"))
            next }
        eqtl_raw = data.table(vroom(eqtl_raw_path, delim = "\t", show_col_types = FALSE, progress = FALSE))
        eqtl_raw_format = select(eqtl_raw, c("Variant","Variant_position","P","Beta","Se"))
        tmp = str_split(eqtl_raw_format$Variant_position, ":", simplify = TRUE)
        eqtl_raw_format$chrom = tmp[, 1]
        eqtl_raw_format$position = as.numeric(str_split(tmp[, 2], "-", simplify = TRUE)[, 1])
        eqtl_raw_format = select(eqtl_raw_format, c("Variant","chrom","position","P","Beta","Se"))
        colnames(eqtl_raw_format) <- c("rsid","chrom","position","pValue","beta","se")
        setindex(eqtl_raw_format, rsid)
        rm(eqtl_raw, tmp)
        gc()

        sentinelVARs_range_path = paste0(range_subdir,"/",gene,"_sentinelVARs_1Mb_range.bed")
        if (!file.exists(sentinelVARs_range_path)){ next }
        sentinelVARs_range = read.table(sentinelVARs_range_path,sep = "\t",header = F)
        colnames(sentinelVARs_range) = c("chr","start","end")
        sentinelVARs_range$chr = paste0("chr",sentinelVARs_range$chr)

        for (m in 1:length(molphe_names)){
            molphe_name = molphe_names[m]

            colocResult_p_All_dir = paste0("02_02_eqtl_coloc_res/",var_type,"-",molphe_name,"/")
            colocResult_p_All_path = paste0(colocResult_p_All_dir,gene,".tsv")
            if (!dir.exists(colocResult_p_All_dir)){ dir.create(colocResult_p_All_dir, recursive = TRUE) }
            if (file.exists(colocResult_p_All_path)){ next } else{
                print(paste0("molphe type ",molphe_name," is processing..."))
            }
            qtl_subdir = paste0(qtl_dir,"/",var_type,"-",molphe_name,"/QTL_results/cis_each")

            overlap_phe_path = paste0(overlap_phe_dir,"/",var_type,"-",molphe_name,"/",gene,".txt")
            if (!file.exists(overlap_phe_path)){ next }
            overlap_phe = read.table(overlap_phe_path,sep = "\t",header = F)$V1

            colocResult_p_All = mclapply(1:length(overlap_phe), function(p) {
                tryCatch({
                    phe = overlap_phe[p]
                    qtl_phe_path = paste0(qtl_subdir, "/", phe, ".txt")
                    qtl_phe = data.table(vroom(qtl_phe_path, delim = "\t", show_col_types = FALSE, progress = FALSE))
                    qtl_flt_format_p = select(qtl_phe, c("Variant", "Phenotype", "P", "Beta", "Se"))
                    colnames(qtl_flt_format_p) <- c("rsid", "gene", "pValue", "beta", "se")

                    eqtl_raw_format_p = eqtl_raw_format[eqtl_raw_format$rsid %in% qtl_flt_format_p$rsid, ]
                    if (nrow(eqtl_raw_format_p) == 0) { return(NULL) }
                    eqtl_raw_format_p$chr = eqtl_raw_format_p$chrom
                    eqtl_raw_format_p$start = eqtl_raw_format_p$position
                    eqtl_raw_format_p$end = eqtl_raw_format_p$position + 1
                    eqtl_raw_format_p_range = filter_molphe_by_range(eqtl_raw_format_p, sentinelVARs_range)                
                    rm(eqtl_raw_format_p)
                    eqtl_raw_format_p_range[, c("chr", "start", "end")] = NULL
                    qtl_flt_format_p_range = qtl_flt_format_p[qtl_flt_format_p$rsid %in% eqtl_raw_format_p_range$rsid, ]
                    qtl_flt_format_p_range = merge(qtl_flt_format_p_range, eqtl_raw_format_p_range[, .(rsid, chrom, position)], by = "rsid")[, .(rsid, chrom, position, pValue, beta, se)]

                    qtl_flt_format_p_range = left_join(qtl_flt_format_p_range, freq_ref, by = "rsid")
                    eqtl_raw_format_p_range = left_join(eqtl_raw_format_p_range, freq_ref, by = "rsid")
                    qtl_flt_format_p_range = select(qtl_flt_format_p_range, c("rsid", "chrom", "position", "pValue", "maf", "beta", "se"))
                    eqtl_raw_format_p_range = select(eqtl_raw_format_p_range, c("rsid", "chrom", "position", "pValue", "maf", "beta", "se"))
                    qtl_flt_format_p_range[qtl_flt_format_p_range$pValue == 0] = 1e-20
                    eqtl_raw_format_p_range[eqtl_raw_format_p_range$pValue == 0] = 1e-20

                    colocResult_p <- suppressMessages(xQTLanalyze_coloc_diy(gwasDF = eqtl_raw_format_p_range, qtlDF = qtl_flt_format_p_range, method = "Both"))
                    rm(eqtl_raw_format_p_range, qtl_flt_format_p_range)
                    colocResult_p <- colocResult_p$coloc_Out_summary
                    colocResult_p = cbind(data.frame(gene = gene, var_type = var_type, molphe_name = molphe_name, phe = phe), colocResult_p)
                    return(colocResult_p)
                }, error = function(e) {
                    cat("Error processing var_type:", var_type, "molphe_name:", molphe_name, "gene:", gene, "phe:", phe, "\n")
                    cat("Error message:", conditionMessage(e), "\n")
                    return(NULL)
                })
            }, mc.cores = mc)
            colocResult_p_All = do.call(rbind, colocResult_p_All)
            colocResult_p_All = colocResult_p_All %>% arrange(gene, var_type, molphe_name, phe)
            fwrite(colocResult_p_All, colocResult_p_All_path,sep = "\t")
            rm(colocResult_p_All)
            gc()
        }
    }
}