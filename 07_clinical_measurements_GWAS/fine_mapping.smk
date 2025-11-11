#! python

# -------------
# FileName     : fine_mapping.smk
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Perform fine mapping for clinical measurements GWAS using SuSiE
# -------------

import os
from glob import glob
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection

PHE_DIR = "/home/wtian/project/HZAU_cohort_meth/wdy_assist_rerun/10_clinic_GWAS/05_conditional_analysis/02_batch/results"
GWAS_FILES = glob(os.path.join(PHE_DIR, "**", "all.ADD.sorted.*.PHENO1.glm.linear"), recursive=True)
PHENOS = [os.path.basename(p).replace(".PHENO1.glm.linear", "").replace("all.ADD.sorted.", "") for p in GWAS_FILES]
PHENO_MAP = dict(zip(PHENOS, GWAS_FILES))
CACHE_FILE = "sig_sv_cache.tsv"
RESULT_DIR = "results"

def load_cache():
    return pd.read_csv(CACHE_FILE, sep="\t")

sv_cache_df = load_cache()

all_targets = [
    os.path.join(RESULT_DIR, row.pheno, f"{row.sv_id}_pip_table.sorted.tsv")
    for idx, row in sv_cache_df.iterrows()
]

rule all:
    input:
        all_targets

rule prepare_region:
    input:
        sumstats = lambda wildcards: PHENO_MAP[wildcards.pheno]
    output:
        snplist = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_snplist.txt")),
        bimlist = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_newsnplist.txt")),
        ld_matrix = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_ld.ld.gz")),
        bim = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_locus.bim")),
        bed = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_locus.bed")),
        fam = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_locus.fam")),
        logf = temp(os.path.join(RESULT_DIR, "{pheno}", "{svid}_locus.log"))
    log:
        os.path.join(RESULT_DIR, "{pheno}", "{svid}_prepare.log")
    shell:
        r"""
        mkdir -p {RESULT_DIR}/{wildcards.pheno}
        Rscript -e '
        library(vroom); library(dplyr)
        sumstats <- vroom::vroom("{input.sumstats}") %>% filter(OBS_CT > 140)
        sv_row <- sumstats %>% filter(ID == "{wildcards.svid}")
        chr <- sv_row$`#CHROM`[1]
        pos <- sv_row$POS[1]
        locus_from <- max(pos - 1e6, 1)
        locus_to <- pos + 1e6
        reg_sumstats <- sumstats %>% filter(`#CHROM` == chr & POS > locus_from & POS < locus_to & !is.na(P))
        write.table(unique(reg_sumstats$ID), file="{output.snplist}", quote=F, row.names=F, col.names=F)
        '
        plink --bfile 03_plink_merge --keep-allele-order --extract {output.snplist} --make-bed --out {RESULT_DIR}/{wildcards.pheno}/{wildcards.svid}_locus >> {log} 2>&1
        plink --bfile {RESULT_DIR}/{wildcards.pheno}/{wildcards.svid}_locus --keep-allele-order -r square gz --out {RESULT_DIR}/{wildcards.pheno}/{wildcards.svid}_ld >> {log} 2>&1
        awk '{{print $2}}' {output.bim} > {output.bimlist}
        """

rule susie_analysis:
    input:
        sumstats = lambda wildcards: PHENO_MAP[wildcards.pheno],
        snplist = os.path.join(RESULT_DIR, "{pheno}", "{svid}_snplist.txt"),
        newsnplist = os.path.join(RESULT_DIR, "{pheno}", "{svid}_newsnplist.txt"),
        ld_matrix = os.path.join(RESULT_DIR, "{pheno}", "{svid}_ld.ld.gz")
    output:
        reg_sumstats = os.path.join(RESULT_DIR, "{pheno}", "{svid}_reg_sumstats.tsv"),
        table = os.path.join(RESULT_DIR, "{pheno}", "{svid}_pip_table.sorted.tsv")
    log:
        os.path.join(RESULT_DIR, "{pheno}", "{svid}_susie.log")
    shell:
        r"""
        Rscript -e '
        library(vroom); library(dplyr); library(data.table); library(susieR)
        sumstats <- vroom("{input.sumstats}") %>% filter(OBS_CT > 140)
        ld_mat <- as.matrix(fread("{input.ld_matrix}", header=F))
        snplist <- read.table("{input.newsnplist}", header=F)$V1
        reg_sumstats <- sumstats[match(snplist, sumstats$ID), ]

        valid_idx <- which(rowSums(is.na(ld_mat)) == 0)
        ld_mat_clean <- ld_mat[valid_idx, valid_idx]
        snplist_clean <- snplist[valid_idx]
        reg_sumstats_clean <- reg_sumstats[match(snplist_clean, reg_sumstats$ID), ]
        fwrite(reg_sumstats_clean, file="{output.reg_sumstats}", sep="\t", quote=F, row.names=F)
        susie_res <- tryCatch({{
            susie_rss(
                bhat = reg_sumstats_clean$BETA,
                shat = reg_sumstats_clean$SE,
                R = ld_mat_clean,
                n = 148,
                L = 10
            )
        }}, error = function(e1) {{
            message("❌ susie_rss() failed with default settings.")
            message("Reason: ", conditionMessage(e1))
            tryCatch({{
                susie_rss(
                    bhat = reg_sumstats_clean$BETA,
                    shat = reg_sumstats_clean$SE,
                    R = ld_mat_clean,
                    n = 148,
                    L = 10,
                    check_prior = FALSE
                )
            }}, error = function(e2) {{
                message("❌ susie_rss() failed again with check_prior=FALSE.")
                message("Reason: ", conditionMessage(e2))
                return(NULL)
            }})
        }})

        if (!is.null(susie_res)) {{
            credible_sets <- susie_res$sets$cs
            pip_table <- data.frame(index = 1:length(snplist_clean), SNP = snplist_clean, PIP = susie_res$pip)
            pip_table <- pip_table[order(-pip_table$PIP), ]
            write.table(pip_table, file="{output.table}", quote=F, row.names=F, sep="\t")
        }} else {{
            write.table(data.frame(), file="{output.table}", quote=F, row.names=F, sep="\t")
        }}
        '
        """