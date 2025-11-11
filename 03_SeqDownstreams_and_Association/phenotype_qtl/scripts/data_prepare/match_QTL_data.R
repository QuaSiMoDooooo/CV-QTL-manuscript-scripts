# -------------
# FileName     : matrixEQTL
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-10-21 21:20
# Last Modified: 2024-11-17 21:38
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

# 加载必要的包
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))

# 读取变异数据
variant_genotype <- fread(snakemake@input[["variant_genotype"]]) %>%
  separate(
    col = "ID",
    into = c("chr", "range", "ID", "ref", "alt"),
    sep = ":"
  ) %>%
  select(-c("chr", "range", "ref", "alt"))
variant_position <- fread(snakemake@input[["variant_position"]])

# 读取表达数据
phenotype_quantity <- fread(snakemake@input[["phenotype_quantity"]])
phenotype_position <- fread(snakemake@input[["phenotype_position"]])

# 读取协变量数据
clinic_cov <- fread(
  snakemake@input[["clinic_covariant"]]
) %>%
  mutate_at(vars("gender"), ~ sapply(., function(x) {
    if (x %in% c(0, 1, 2)) {
      return(x)
    }
    if (x == "男") {
      return(1)
    } else if (x == "女") {
      return(2)
    } else {
      return(0)
    }
  }))
peer_cov <- fread(snakemake@input[["peer_covariant"]])
pop_cov <- fread(
  snakemake@input[["pop_covariant"]],
  sep = " ",
  col.names = c("FID", "IID", paste0("PC", 1:5))
)


# 用基因型的样本顺序作为最终的样本顺序
sample_order <- intersect(
  colnames(variant_genotype)[-1],
  colnames(phenotype_quantity)[-1]
)

# match 基因型
variant_genotype_matched <- variant_genotype %>%
  select(1, all_of(sample_order))

# match 表达量
phenotype_quantity_matched <- phenotype_quantity %>%
  select(1, all_of(sample_order))

# match协变量
clinic_cov_matched <- clinic_cov %>%
  column_to_rownames("ID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cov") %>%
  select(1, all_of(sample_order))
peer_cov_matched <- peer_cov %>%
  select(1, all_of(sample_order))
pop_cov_matched <- pop_cov %>%
  select(-FID) %>%
  column_to_rownames("IID") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("cov") %>%
  select(1, all_of(sample_order))
cov_merged_matched <- rbind(
  clinic_cov_matched,
  pop_cov_matched,
  peer_cov_matched
)

print(paste0("基因型数量：", nrow(variant_genotype_matched)))
print(paste0("基因型位置：", nrow(variant_position)))
print(paste0(
  "基因型重叠数量：",
  length(intersect(
    variant_genotype_matched[[1]],
    variant_position[[1]]
  ))
))
print(paste0("表型数量：", nrow(phenotype_quantity_matched)))
print(paste0("表型位置：", nrow(phenotype_position)))
print(paste0(
  "表型重叠数量：",
  length(intersect(
    phenotype_quantity_matched[[1]],
    phenotype_position[[1]]
  ))
))
# 确定所有的基因型和表型都有位置
if (!identical(
  variant_genotype_matched[[1]],
  variant_position[[1]]
)) {
  message("ID of variant genotype not match with variant position, remapping")
  if (!all(variant_position[[1]] %in% variant_genotype_matched[[1]])) {
    stop("some variant genotype not in variant genotype")
  } else {
    variant_genotype_matched <- variant_genotype_matched %>%
      filter(ID %in% variant_position[[1]])
    variant_genotype_matched <- variant_genotype_matched[match(
      variant_position[[1]],
      variant_genotype_matched[[1]]
    ), ]
  }
}
if (!identical(
  phenotype_quantity_matched[[1]],
  phenotype_position[[1]]
)) {
  message(
    "ID of phenotype quantity not match with phenotype position, remapping"
  )
  overlaped_phenotype <- intersect(
    phenotype_quantity_matched[[1]],
    phenotype_position[[1]]
  )
  phenotype_quantity_matched <- phenotype_quantity_matched[match(
    overlaped_phenotype,
    phenotype_quantity_matched[[1]]
  ), ]
}

# 输出文件
fwrite(
  variant_genotype_matched,
  snakemake@output[["genotype"]],
  sep = "\t",
  quote = FALSE
)
fwrite(
  phenotype_quantity_matched,
  snakemake@output[["phenotype"]],
  sep = "\t",
  quote = FALSE
)
fwrite(
  cov_merged_matched,
  snakemake@output[["covariant"]],
  sep = "\t",
  quote = FALSE
)
