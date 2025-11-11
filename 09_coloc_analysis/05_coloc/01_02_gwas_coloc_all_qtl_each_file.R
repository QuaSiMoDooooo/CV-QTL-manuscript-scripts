#!/usr/bin/env Rscript

# -------------
# FileName     : 01_02_gwas_coloc_all_qtl_each_file.R
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Perform colocalization analysis between GWAS and QTL data for each molecular phenotype
# -------------

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(xQTLbiolinks))
suppressPackageStartupMessages(library(parallel))
setwd("/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc")

mc <- getOption("mc.cores", 64)

gwas_flt_dir = "../data/GWAS_flt"
gwas_flt_files = list.files(gwas_flt_dir, pattern = "\\.tsv$", full.names = TRUE)
gwas_var_type_list = str_split(basename(gwas_flt_files), "\\.", simplify = TRUE)[, 1]
var_type_list = c("SNP","SV")
gwas_flt_files = gwas_flt_files[gwas_var_type_list %in% var_type_list]
gwas_var_type_list = gwas_var_type_list[gwas_var_type_list %in% var_type_list]

gwas_raw_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA"

range_dir = "../01_get_chr1-22_GWAS_and_eQTL_sentinelVARs_range/"

qtl_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/04_rerun_QTL/results"
molphe_names = c("expression","splicing","APA","methylation")

overlap_phe_dir = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/11_coloc/05_coloc/01_01_overlap_phe"

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

for (j in 1:length(gwas_var_type_list)){
  var_type = gwas_var_type_list[j]
  range_subdir = paste0(range_dir,"GWAS_",var_type)
  gwas_flt_file = gwas_flt_files[j]
  gwas_flt = data.table(vroom(gwas_flt_file, show_col_types = FALSE, progress = FALSE))
  gwas_flt = gwas_flt[gwas_flt$IndicatorType == "common", ]
  traits = unique(gwas_flt$Indicator)

  for(t in 1:length(traits)){
    trait = traits[t]

    colocResult_p_All_paths = paste0("01_02_gwas_coloc_res/",var_type,"-",molphe_names,"/",trait,".tsv")
    if (all(file.exists(colocResult_p_All_paths))){
      next
    }
    
    gwas_raw_path = paste0(gwas_raw_dir,"/",var_type,"-biochemistry/biochem_common/",trait,".PHENO1.glm.linear")
    gwas_raw = data.table(vroom(gwas_raw_path, delim = "\t", show_col_types = FALSE, progress = FALSE))
    gwas_raw_format = filter(gwas_raw, TEST=="ADD") %>%  select(c("ID","#CHROM","POS","P","A1_FREQ","BETA","SE")) 
    colnames(gwas_raw_format) <- c("rsid","chrom","position","pValue","AF","beta","se")
    setindex(gwas_raw_format, rsid)
    rm(gwas_raw)

    sentinelVARs_range_path = paste0(range_subdir,"/",trait,"_sentinelVARs_1Mb_range.bed")
    sentinelVARs_range = read.table(sentinelVARs_range_path,sep = "\t",header = F)
    colnames(sentinelVARs_range) = c("chr","start","end")
    sentinelVARs_range$chr = paste0("chr",sentinelVARs_range$chr)

    for (m in 1:length(molphe_names)){
      molphe_name = molphe_names[m]
      
      colocResult_p_All_dir = paste0("01_02_gwas_coloc_res/",var_type,"-",molphe_name,"/")
      colocResult_p_All_path = paste0(colocResult_p_All_dir,trait,".tsv")
      if (!dir.exists(colocResult_p_All_dir)){ dir.create(colocResult_p_All_dir, recursive = TRUE) }
      if (file.exists(colocResult_p_All_path)){ next }
      qtl_subdir = paste0(qtl_dir,"/",var_type,"-",molphe_name,"/QTL_results/cis_each")

      overlap_phe_path = paste0(overlap_phe_dir,"/",var_type,"-",molphe_name,"/",trait,".txt")
      if (!file.exists(overlap_phe_path)){ next }
      overlap_phe = read.table(overlap_phe_path,sep = "\t",header = F)$V1

      colocResult_p_All = mclapply(1:length(overlap_phe), function(p){
        phe = overlap_phe[p]
        qtl_phe_path = paste0(qtl_subdir,"/",phe,".txt")
        qtl_phe = data.table(vroom(qtl_phe_path, delim = "\t", show_col_types = FALSE, progress = FALSE))
        qtl_flt_format_p = select(qtl_phe, c("Variant","Phenotype","P","Beta","Se"))
        colnames(qtl_flt_format_p) <- c("rsid","gene","pValue","beta","se")
        gwas_raw_format_p = gwas_raw_format[gwas_raw_format$rsid %in% qtl_flt_format_p$rsid, ]
        if (nrow(gwas_raw_format_p) == 0){ return(NULL) }
        gwas_raw_format_p$chr = paste0("chr",gwas_raw_format_p$chrom)
        gwas_raw_format_p$start = gwas_raw_format_p$position
        gwas_raw_format_p$end = gwas_raw_format_p$position+1
        gwas_raw_format_p_range = filter_molphe_by_range(gwas_raw_format_p, sentinelVARs_range)
        rm(gwas_raw_format_p)
        gwas_raw_format_p_range[,c("chr","start","end")] = NULL
        qtl_flt_format_p_range = qtl_flt_format_p[qtl_flt_format_p$rsid %in% gwas_raw_format_p_range$rsid, ]
        qtl_flt_format_p_range = merge(qtl_flt_format_p_range, gwas_raw_format_p_range[, .(rsid, chrom, position)], by = "rsid")[, .(rsid, chrom, position, pValue, beta, se)]
        qtl_flt_format_p_range$maf = gwas_raw_format_p_range[match(qtl_flt_format_p_range$rsid, gwas_raw_format_p_range$rsid), "AF"]
        qtl_flt_format_p_range = select(qtl_flt_format_p_range, c("rsid","chrom","position","pValue","maf","beta","se"))
        qtl_flt_format_p_range[qtl_flt_format_p_range$pValue ==0 ] = 1e-20
        colocResult_p <- suppressMessages(xQTLanalyze_coloc_diy(gwasDF = gwas_raw_format_p_range, qtlDF = qtl_flt_format_p_range, method = "Both"))
        rm(gwas_raw_format_p_range, qtl_flt_format_p_range)
        colocResult_p <- colocResult_p$coloc_Out_summary
        colocResult_p = cbind(data.frame(trait = trait, var_type = var_type, molphe_name = molphe_name, phe = phe), colocResult_p)
        return(colocResult_p)
      }, mc.cores = mc)
      colocResult_p_All = do.call(rbind, colocResult_p_All)
      colocResult_p_All = colocResult_p_All %>% arrange(trait, var_type, molphe_name, phe)
      fwrite(colocResult_p_All, colocResult_p_All_path,sep = "\t")
      rm(colocResult_p_All)
      gc()
    }
  }
}