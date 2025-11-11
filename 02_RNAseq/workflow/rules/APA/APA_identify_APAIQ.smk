# -*- coding: UTF-8 -*-
#
# FileName     : APA_identify_APAIQ
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-11 21:08
# Last Modified: 2024-04-11 21:08
# Modified By  : EastsunW
# -------------
# Description  : 使用APAIQ鉴定APA，需要合并
# Citation     : https://github.com/christear/APAIQ_release
# -------------


rule APAIQ_bam2bedgraph:
    input:
        rules.alignment.output.sorted,
    output:
        temp(f"results/{{cohort}}/APA/APAIQ/{{sample}}.bedgraph"),
    conda:
        "apaiq"
    shell:
        f"""
            genomeCoverageBed -split -bg -ibam {{input}} > {{output}}
        """


rule APAIQ_identify:
    input:
        rules.APAIQ_bam2bedgraph.output,
    output:
        f"results/{{cohort}}/APA/APAIQ/sample/{{sample}}.predicted.txt",
    log:
        f"logs/{{cohort}}/APA/{{sample}}.APAIQ_apapiq.log",
    benchmark:
        f"benchmarks/{{cohort}}/APA/{{sample}}.APAIQ_apapiq.benchmark"
    conda:
        "apaiq"
    params:
        reference=cohort.ref_config["alignment"]["fasta"],
        model=cohort.ref_config["APA"]["APAIQ"]["model"],
        pas_db=cohort.ref_config["APA"]["APAIQ"]["pas_db"],
        extra=config["APA"]["APAIQ"]["extra_params"],
    threads: config["APA"]["APAIQ"]["threads_APAIQ"]
    shell:
        f"""
        (
            apaiq \
                --t {{threads}} \
                --input_file {{input}} \
                --out_dir results/{{wildcards.cohort}}/APA/APAIQ/sample \
                --name {{wildcards.sample}} \
                --fa_file {{params.reference}} \
                --DB_file {{params.pas_db}} \
                --model {{params.model}} \
                {{params.extra}}
        ) &> {{log}}
        """


rule APAIQ_merge:
    input:
        expand(
            "results/{{cohort}}/APA/APAIQ/sample/{sample}.predicted.txt",
            sample=cohort.list_samples(),
        ),
    output:
        f"results/{{cohort}}/APA/{{cohort}}.APAIQ.txt",
    script:
        "APAIQ_merge.R"
