# -*- coding: UTF-8 -*-
#
# FileName     : process_GWAS
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-19 23:28
# Last Modified: 2024-11-20 19:04
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule process_GWAS:
    input:
        genotype=multiext(
            "results/{variant}-biochemistry/{variant}.biochem_{biochem_type}",
            ".bed",
            ".bim",
            ".fam",
        ),
        phenotype="results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.pheno",
        marker="results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.marker",
        covariant="results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.covariant",
    output:
        directory("results/{variant}-biochemistry/biochem_{biochem_type}"),
    log:
        "logs/GWAS/{variant}-{biochem_type}/.log",
    threads: 10
    run:
        import subprocess

        def run_plink(index, marker):
            cmd = (
                f"mkdir -p results/{wildcards.variant}-biochemistry/biochem_{wildcards.biochem_type} && "
                f"scripts/tools/plink-2 "
                f"--silent --neg9-pheno-really-missing --chr 1-22 "
                f"--threads {threads} "
                f"--bfile results/{wildcards.variant}-biochemistry/{wildcards.variant}.biochem_{wildcards.biochem_type} "
                f"--pheno {input.phenotype} "
                f"--covar {input.covariant} "
                f"--allow-no-sex --glm omit-ref --real-ref-alleles "
                f"--pheno-col-nums {str(index)} "
                f"--out results/{wildcards.variant}-biochemistry/biochem_{wildcards.biochem_type}/{marker} "
            )
            with open(f"logs/GWAS/{wildcards.variant}-{wildcards.biochem_type}.log", "w") as log:
                subprocess.run(cmd, shell=True, check=True, stdout=log, stderr=log)

        with open(input["marker"], "rt") as MARKER:
            markers = [line.strip().split("\t") for line in MARKER]

        for index, marker in markers:
            run_plink(int(index)+2, marker)
