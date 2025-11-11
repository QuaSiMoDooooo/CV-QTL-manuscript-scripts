#! python

# -------------
# FileName     : conditional_analysis.smk
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Perform conditional analysis for clinical measurements GWAS using GCTA
# -------------

import os
from glob import glob

SNP_DIR = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA/SNP-biochemistry/biochem_common"
SV_DIR = "/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA/SV-biochemistry/biochem_common"

def extract_pheno(path):
    return os.path.basename(path).replace(".PHENO1.glm.linear", "")

SNP_FILES = glob(os.path.join(SNP_DIR, "*.PHENO1.glm.linear"))
SV_FILES = glob(os.path.join(SV_DIR, "*.PHENO1.glm.linear"))

PHENOS = sorted(set(map(extract_pheno, SNP_FILES)) & set(map(extract_pheno, SV_FILES)))

rule all:
    input:
        expand("results/{pheno}/05_{pheno}.sig.sv.gwas.cond.tsv", pheno=PHENOS)

rule merge_gwas:
    input:
        snp = lambda wc: os.path.join(SNP_DIR, f"{wc.pheno}.PHENO1.glm.linear"),
        sv = lambda wc: os.path.join(SV_DIR, f"{wc.pheno}.PHENO1.glm.linear")
    output:
        merged = "results/{pheno}/all.ADD.sorted.{pheno}.PHENO1.glm.linear",
        snpadd = "results/{pheno}/SNP.ADD.{pheno}.PHENO1.glm.linear",
        svadd = "results/{pheno}/SV.ADD.{pheno}.PHENO1.glm.linear"
    log:
        "results/{pheno}/log.merge_gwas.txt"
    shell:
        """
        mkdir -p results/{wildcards.pheno}
        awk 'NR==1 || /ADD/ {{print > "{output.svadd}"}}' {input.sv}
        awk 'NR==1 || /ADD/ {{print > "{output.snpadd}"}}' {input.snp}
        cat {output.snpadd} > tmp_all.ADD.{wildcards.pheno}.PHENO1.glm.linear
        awk -v OFS="\t" 'NR>1{{print $0}}' {output.svadd} >> tmp_all.ADD.{wildcards.pheno}.PHENO1.glm.linear
        {{ head -n 1 tmp_all.ADD.{wildcards.pheno}.PHENO1.glm.linear; sed '1d' tmp_all.ADD.{wildcards.pheno}.PHENO1.glm.linear | sort -k1,1n -k2,2n; }} > {output.merged}
        rm tmp_all.ADD.{wildcards.pheno}.PHENO1.glm.linear
        &> {log}
        """

rule format_ma:
    input:
        "results/{pheno}/all.ADD.sorted.{pheno}.PHENO1.glm.linear"
    output:
        "results/{pheno}/02_all.ADD.sorted.{pheno}.PHENO1.glm.linear.OBSflt141.ma"
    log:
        "results/{pheno}/log.format_ma.txt"
    run:
        import pandas as pd, sys
        sys.stdout = open(log[0], "w")
        sys.stderr = sys.stdout
        df = pd.read_csv(input[0], sep="\t")
        ma = pd.DataFrame({
            "SNP": df["ID"], "A1": df["ALT"], "A2": df["REF"],
            "freq": df["A1_FREQ"], "b": df["BETA"],
            "se": df["SE"], "p": df["P"], "N": df["OBS_CT"]
        }).dropna(subset=["b"])
        ma = ma[ma["N"] > 140]
        ma.to_csv(output[0], sep="\t", index=False)

rule plink_merge:
    input:
        "03_plink_merge_list.txt"
    output:
        "03_plink_merge.bed",
        "03_plink_merge.bim",
        "03_plink_merge.fam"
    log:
        "log.plink_merge.txt"
    shell:
        "plink --merge-list {input} --keep-allele-order --make-bed --out 03_plink_merge &> {log}"

rule gcta_cojo:
    input:
        ma = "results/{pheno}/02_all.ADD.sorted.{pheno}.PHENO1.glm.linear.OBSflt141.ma",
        bed = "03_plink_merge.bed",
        bim = "03_plink_merge.bim",
        fam = "03_plink_merge.fam"
    output:
        cma = "results/{pheno}/04_gcta_cojo_slct.cma.cojo",
        jma = "results/{pheno}/04_gcta_cojo_slct.jma.cojo",
        ldr = "results/{pheno}/04_gcta_cojo_slct.ldr.cojo"
    log:
        "results/{pheno}/log.gcta_cojo.txt"
    run:
        import subprocess, sys
        sys.stdout = open(log[0], "w")
        sys.stderr = sys.stdout
        out_prefix = f"results/{wildcards.pheno}/04_gcta_cojo_slct"
        cmd_basic = f"gcta64 --bfile 03_plink_merge --cojo-file {input.ma} --cojo-slct --cojo-p 1e-04 --out {out_prefix}"
        cmd_fallback = f"gcta64 --bfile 03_plink_merge --cojo-file {input.ma} --cojo-slct --cojo-top-SNPs 140 --out {out_prefix}"
        try:
            subprocess.run(cmd_basic, shell=True, check=True, capture_output=True)
        except subprocess.CalledProcessError:
            subprocess.run(cmd_fallback, shell=True, check=True, capture_output=True)

rule sv_gwas_cond:
    input:
        sv = "results/{pheno}/SV.ADD.{pheno}.PHENO1.glm.linear",
        cond = "results/{pheno}/04_gcta_cojo_slct.cma.cojo"
    output:
        "results/{pheno}/05_{pheno}.sig.sv.gwas.cond.tsv"
    log:
        "results/{pheno}/log.sv_gwas_cond.txt"
    run:
        import pandas as pd, sys
        from statsmodels.stats.multitest import multipletests
        sys.stdout = open(log[0], "w")
        sys.stderr = sys.stdout
        sv = pd.read_csv(input.sv, sep="\t")
        sv.dropna(subset=["P"], inplace=True)
        sv["fdr"] = multipletests(sv["P"], method="fdr_bh")[1]
        sig = sv[sv["fdr"] < 0.05]
        cond = pd.read_csv(input.cond, sep="\t")
        result = cond[cond["SNP"].isin(sig["ID"])]
        result.to_csv(output[0], sep="\t", index=False)