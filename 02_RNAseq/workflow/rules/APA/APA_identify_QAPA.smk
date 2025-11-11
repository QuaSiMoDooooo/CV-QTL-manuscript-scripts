# -*- coding: UTF-8 -*-
#
# FileName     : APA_identify_QAOA
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-11 21:37
# Last Modified: 2024-04-11 21:37
# Modified By  : EastsunW
# -------------
# Description  : 使用QAPA完成APA的鉴定
# Citation     : https://github.com/morrislab/qapa
# -------------


rule QAPA_UTR_quant:
    input:
        seq_1=lambda wildcards: f"{cohort.input_dir}/{cohort.list_files(wildcards.sample)[0]}",
        seq_2=lambda wildcards: f"{cohort.input_dir}/{cohort.list_files(wildcards.sample)[1]}",
    output:
        temp("results/{cohort}/APA/QAPA/UTR_quant/{sample}/quant.sf"),
    threads: config["APA"]["QAPA"]["threads_salmon"]
    conda: "../../envs/expression.yaml"
    log:
        "logs/{cohort}/APA/QAPA/UTR_quant/{sample}.log",
    params:
        utr_index=cohort.ref_config["APA"]["QAPA"]["utr_index"],
    shell:
        """
        (
            salmon quant \
                --gcBias -l A \
                --index {params.utr_index} \
                -1 {input.seq_1} \
                -2 {input.seq_2} \
                --threads {threads} \
                --output results/{wildcards.cohort}/APA/QAPA/UTR_quant/{wildcards.sample}
        ) &> {log}
        """


rule QAPA_identify:
    input:
        lambda wildcards: expand(
            f"results/{wildcards.cohort}/APA/QAPA/UTR_quant/{{sample}}/quant.sf",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/APA/{cohort}.QAPA.txt",
    log:
        "logs/{cohort}/APA/QAPA/QAPA_identify/{cohort}.APA_identify.log",
    benchmark:
        "benchmarks/{cohort}/APA/QAPA/QAPA_identify/{cohort}.APA_identify.benchmark"
    threads: config["APA"]["QAPA"]["threads_QAPA"]
    container:
        "docker://dockerproxy.net/sambrycesmith/qapa_fork:408d4ae"
    params:
        identifiers=cohort.ref_config["APA"]["QAPA"]["identifiers"],
        extra=config["APA"]["QAPA"]["extra_params"],
    shell:
        """
        (
            qapa quant \
                --db {params.identifiers} \
                {input} > \
                {output}
        ) &> {log}
        """
