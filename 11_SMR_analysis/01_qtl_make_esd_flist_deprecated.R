suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(vroom))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(parallel))

cur_path <- "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/07_QTL_GWAS_SMR"
setwd(cur_path)

# 读取rsid注释信息
rsid <- fread("./simplify_rsid_ref_bed/all_qtl_filtered.sorted.uniq.rsid.redundant.bed")
judge_vec <- !duplicated(paste(rsid$V1, rsid$V2, sep = "_"))
rsid_ref <- data.frame(chr = rsid$V1, pos = rsid$V2, rsid = rsid$V4)[judge_vec, ]

# 读取frq信息
frq <- fread("/home/wtian/project/HZAU_cohort_meth/wdy_result/SNP_plink_genotype_process_bim3/SNP.biochem_common.frq")
tmp <- str_split(frq$SNP, "_", simplify = TRUE)
frq$chr_pos <- paste(tmp[, 1], tmp[, 2], sep = "_")
rsid_ref$chr_pos <- paste(rsid_ref$chr, rsid_ref$pos, sep = "_")
frq_ref <- frq[frq$chr_pos %in% rsid_ref$chr_pos, ]
rsid_ref$chr_pos <- NULL
tmp <- select(frq_ref, chr_pos, MAF)
tmp$chr <- str_split(tmp$chr_pos, "_", simplify = TRUE)[, 1]
tmp$pos <- str_split(tmp$chr_pos, "_", simplify = TRUE)[, 2]
frq_ref <- select(tmp, chr, pos, MAF) %>% filter(MAF > 0)
frq_ref$chr_pos <- paste(frq_ref$chr, frq_ref$pos, sep = "_")
rsid_ref$chr_pos <- paste(rsid_ref$chr, rsid_ref$pos, sep = "_")
rsid_ref <- rsid_ref %>% select(chr_pos, rsid)
frq_ref <- frq_ref %>% select(chr_pos, MAF)

# 任务函数
process_qtl_type <- function(qtl_info) {
  qtl <- qtl_info$qtl
  qtl_name <- qtl_info$qtl_name
  type <- qtl_info$type
  qtl_dir <- "/home/wangdy/Projects/Weibin/Downstreams/phenotype_qtl/results"

  file_path <- file.path(qtl_dir, paste0("SNP-", qtl_name), "QTL_results", paste0(type, ".filtered.txt.gz"))
  message("Processing: ", file_path)

  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }

  data <- tryCatch({
    fread(file_path, showProgress = TRUE, sep = "\t")
  }, error = function(e) {
    warning("Error reading file: ", file_path)
    return(NULL)
  })

  tmp <- str_split(data$Variant, "_", simplify = TRUE)
  data$chr_pos <- paste(tmp[, 1], tmp[, 2], sep = "_")
  data$chr <- tmp[, 1]
  data$pos <- tmp[, 2]

  data <- left_join(data, rsid_ref, by = "chr_pos") %>% filter(!is.na(rsid))
  data <- left_join(data, frq_ref, by = "chr_pos") %>% filter(!is.na(MAF))

  data_slt <- data %>%
    select(chr, pos, rsid, Phenotype, Alt, Ref, MAF, Beta, P, Se) %>%
    setNames(c("Chr", "BP", "SNP", "phe", "A1", "A2", "Freq", "Beta", "p", "se")) %>%
    select(phe, Chr, SNP, BP, A1, A2, Freq, Beta, se, p) %>%
    arrange(phe, Chr, BP) %>%
    mutate(Chr = sub("chr", "", Chr))

  # 保存esd
  esd_dir <- file.path(cur_path, "01_esd", paste0(type, "_", qtl))
  if (!dir.exists(esd_dir)) dir.create(esd_dir, recursive = TRUE)

  phe_list <- unique(data_slt$phe)
  for (phe in phe_list) {
    file_name <- file.path(esd_dir, paste0(phe, ".esd"))
    write.table(data_slt[data_slt$phe == phe, -1], file_name, quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  }

  # 构造flist
  flist_dir <- file.path(cur_path, "01_flist")
  if (!dir.exists(flist_dir)) dir.create(flist_dir)

  flist_name <- file.path(flist_dir, paste0(type, "_", qtl, ".flist"))
  tmp2 <- data[!duplicated(data$Phenotype), ]
  tmp2$ProbeBp <- as.numeric(tmp2$pos) - as.numeric(tmp2$distance)
  tmp2$phechr <- sub("chr", "", str_split(tmp2$Phenotype_position, ":", simplify = TRUE)[, 1])
  flist <- data.frame(
    ProbeID = tmp2$Phenotype,
    Chr = tmp2$phechr,
    ProbeBp = tmp2$ProbeBp,
    Gene = tmp2$Gene,
    GeneticDistance = 0,
    Orientation = NA
  )
  flist$PathOfEsd <- file.path(esd_dir, paste0(flist$ProbeID, ".esd"))
  flist <- flist %>% select(Chr, ProbeID, GeneticDistance, ProbeBp, Gene, Orientation, PathOfEsd)

  if (all(is.na(flist$Gene))) flist$Gene <- NA
  write.table(flist, flist_name, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  message("Finished: ", type, "_", qtl)
  return(flist_name)
}

# 构造参数列表
qtl_list <- c("eQTL", "sQTL", "apaQTL", "meQTL")
qtl_name_list <- c("expression", "splicing", "APA", "methylation")
type_list <- c("cis", "trans")

param_list <- expand.grid(qtl = qtl_list, type = type_list, stringsAsFactors = FALSE)
param_list$qtl_name <- qtl_name_list[match(param_list$qtl, qtl_list)]
param_list <- split(param_list, seq(nrow(param_list)))

# 并行运行
result_list <- mclapply(param_list, process_qtl_type, mc.cores = 32)
