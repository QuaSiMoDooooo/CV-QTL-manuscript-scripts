# -*- coding: UTF-8 -*-
#
# FileName     : prepare_phenotype
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-20 10:51
# Last Modified: 2024-09-20 15:23
# Modified By  : EastsunW
# -------------
# Description  : 准备用于QTL的表型文件，包括表型的质量控制、样本ID转换和分位数归一化
# -------------


rule phenotype_quantity_filter_expression:
    input:
        "resources/phenotypes/expression.quantity.txt",
    output:
        "results/common_data/phenotype/expression.quantity.filtered.txt",
    params:
        phenotype="expression",
        average_expression=config["data_filter"]["expression"]["average_expression"],
        median_expression=config["data_filter"]["expression"]["median_expression"],
    shell:
        """
        python scripts/data_prepare/filter_phenotype.py {params.phenotype} \
            -i {input} \
            -o {output} \
            --average {params.average_expression} \
            --median {params.median_expression}
        """


rule phenotype_quantity_filter_splicing:
    input:
        "resources/phenotypes/splicing.quantity.txt",
    output:
        "results/common_data/phenotype/splicing.quantity.filtered.txt",
    params:
        phenotype="splicing",
        max_na=config["data_filter"]["splicing"]["max_na_rate"],
    shell:
        """
        python scripts/data_prepare/filter_phenotype.py {params.phenotype} \
            -i {input} \
            -o {output} \
            --na {params.max_na}
        """


rule phenotype_quantity_filter_APA:
    input:
        "resources/phenotypes/APA.quantity.txt",
    output:
        "results/common_data/phenotype/APA.quantity.filtered.txt",
    params:
        phenotype="APA",
        max_na=config["data_filter"]["APA"]["max_na_rate"],
    shell:
        """
        python scripts/data_prepare/filter_phenotype.py {params.phenotype} \
            -i {input} \
            -o {output} \
            --na {params.max_na}
        """


rule phenotype_quantity_filter_methylation:
    input:
        "resources/phenotypes/methylation.quantity.txt"
    output:
        "results/common_data/phenotype/methylation.quantity.filtered.txt",
    shell:
        """
        cp {input} {output}
        """



rule phenotype_quantity_qnorm:
    input:
        "results/common_data/phenotype/{phenotype}.quantity.filtered.txt",
    output:
        "results/common_data/phenotype/{phenotype}.quantity.qnorm.txt",
    script:
        "../scripts/tools/rank_normalize.R"


rule phenotype_position:
    input:
        "resources/phenotypes/{phenotype}.position.raw.txt",
    output:
        "results/common_data/phenotype/{phenotype}.position.txt",
    shell:
        """
        python scripts/data_prepare/extract_phenotype_position.py {wildcards.phenotype} \
            -i {input} \
            -o {output} \
            --exclude_YM
        """


# 归一化加格式转换
rule biochemistry_to_plink:
    input:
        "resources/phenotypes/biochemistry.{biochem_type}.txt",
    output:
        quantity="results/common_data/phenotype/biochem_{biochem_type}.norm.txt",
        marker="results/common_data/phenotype/biochem_{biochem_type}.marker.txt",
    script:
        "../scripts/data_prepare/normalize_biochemistry.R"


for variant in ["SNP", "InDel", "MNV", "SV"]:

    rule:
        name:
            f"copy_phenotype_files_{variant}"
        input:
            quantity=rules.phenotype_quantity_qnorm.output,
            position=rules.phenotype_position.output,
        output:
            target_quantity=f"results/{variant}-{{phenotype}}/filtered_data/phenotype.quantity.qnorm.txt",
            target_position=f"results/{variant}-{{phenotype}}/matched_data/phenotype.position.txt",
        shell:
            """
            cp {input.quantity} {output.target_quantity}
            cp {input.position} {output.target_position}
            """

    rule:
        name:
            f"copy_biochemistry_files_{variant}"
        input:
            quantity="results/common_data/phenotype/biochem_{biochem_type}.norm.txt",
            marker="results/common_data/phenotype/biochem_{biochem_type}.marker.txt",
        output:
            quantity=f"results/{variant}-biochemistry/{variant}.biochem_{{biochem_type}}.pheno",
            marker=f"results/{variant}-biochemistry/{variant}.biochem_{{biochem_type}}.marker",
        shell:
            """
            cp {input.quantity} {output.quantity}
            cp {input.marker} {output.marker}
            """
