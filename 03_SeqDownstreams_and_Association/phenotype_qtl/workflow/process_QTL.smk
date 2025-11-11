# -*- coding: UTF-8 -*-
#
# FileName     : multithread_QTL copy
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-15 17:04
# Last Modified: 2024-11-15 23:44
# Modified By  : EastsunW
# -------------
# Description  : 进行QTL分析
# -------------


# 匹配QTL的数据
rule match_QTL_data:
    input:
        variant_genotype="results/{variant}-{phenotype}/filtered_data/variant.genotype.filtered.txt",
        variant_position="results/{variant}-{phenotype}/matched_data/variant.position.txt",
        phenotype_quantity="results/{variant}-{phenotype}/filtered_data/phenotype.quantity.qnorm.txt",
        phenotype_position="results/{variant}-{phenotype}/matched_data/phenotype.position.txt",
        clinic_covariant="results/{variant}-{phenotype}/filtered_data/clinic_covariant.txt",
        pop_covariant="results/{variant}-{phenotype}/filtered_data/pop_covariant.txt",
        peer_covariant="results/{variant}-{phenotype}/filtered_data/peer_covariant.txt",
    output:
        genotype="results/{variant}-{phenotype}/matched_data/variant.genotype.matched.txt",
        phenotype="results/{variant}-{phenotype}/matched_data/phenotype.quantity.matched.txt",
        covariant="results/{variant}-{phenotype}/matched_data/covariant.txt",
    script:
        "../scripts/data_prepare/match_QTL_data.R"


rule QTL_analysis:
    input:
        variant_genotype="results/{variant}-{phenotype}/matched_data/variant.genotype.matched.txt",
        variant_position="results/{variant}-{phenotype}/matched_data/variant.position.txt",
        phenotype_quantity="results/{variant}-{phenotype}/matched_data/phenotype.quantity.matched.txt",
        phenotype_position="results/{variant}-{phenotype}/matched_data/phenotype.position.txt",
        covariant="results/{variant}-{phenotype}/matched_data/covariant.txt",
    output:
        [
            "results/{variant}-{phenotype}/QTL_results/sliced/raw.txt_" + str(i)
            for i in range(1, int(config["QTL"]["threads"]) + 1)
        ],
    params:
        pvalue=config["QTL"]["pvalue"],
        cis_window=config["QTL"]["cis_window"],
    threads: config["QTL"]["threads"]
    shell:
        """
        Rscript scripts/tools/matrixEQTL_slice.R \
            --thread {threads} \
            --variant_genotype {input.variant_genotype} \
            --variant_position {input.variant_position} \
            --phenotype_quantity {input.phenotype_quantity} \
            --phenotype_position {input.phenotype_position} \
            --covariate {input.covariant} \
            --pvalue {params.pvalue} \
            --output results/{wildcards.variant}-{wildcards.phenotype}/QTL_results/sliced/raw.txt
        """


rule merge_QTL_results:
    input:
        expand(
            "results/{{variant}}-{{phenotype}}/QTL_results/sliced/raw.txt_{i}",
            i=list(range(1, int(config["QTL"]["threads"]) + 1)),
        ),
    output:
        "results/{variant}-{phenotype}/QTL_results/raw.merged.txt",
    shell:
        """
        head -n 1 {input[0]} > {output}
        tail -n +2 -q {input} >> {output}
        """


rule format_qtl_results:
    input:
        qtl_results="results/{variant}-{phenotype}/QTL_results/raw.merged.txt",
        variant_position="results/{variant}-{phenotype}/matched_data/variant.position.txt",
        phenotype_position="results/{variant}-{phenotype}/matched_data/phenotype.position.txt",
    output:
        cis="results/{variant}-{phenotype}/QTL_results/cis.txt",
        trans="results/{variant}-{phenotype}/QTL_results/trans.txt",
    shell:
        """
        Rscript scripts/tools/format_QTL_results.R \
            --qtl {input.qtl_results} \
            --variant_position {input.variant_position} \
            --phenotype_position {input.phenotype_position} \
            --o_cis {output.cis} \
            --o_trans {output.trans}
        """


rule zip_qtl_results:
    input:
        "results/{variant}-{phenotype}/QTL_results/{qtl_type}.txt",
    output:
        "results/{variant}-{phenotype}/QTL_results/{qtl_type}.txt.gz",
    shell:
        """
        gzip -c {input} > {output}
        """


rule filter_qtl_results:
    input:
        "results/{variant}-{phenotype}/QTL_results/{qtl_type}.txt",
    output:
        "results/{variant}-{phenotype}/QTL_results/{qtl_type}.filtered.txt.gz",
    params:
        fdr=config["QTL"]["fdr"],
    shell:
        """
        awk -F '\t' \
        -v fdr={params.fdr} \
        'BEGIN {{OFS="\\t"}} NR==1 || ($11<fdr)' \
        {input} \
        | gzip -c > \
        {output}
        """
