# -*- coding: UTF-8 -*-
#
# FileName     : SV_identify_SVision copy
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-19 09:21
# Last Modified: 2024-11-19 09:21
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule SV_SVisionpro_call:
    input:
        bam=rules.alignment_merge.output.bam,
        genome="resources/reference/alignment/GRCh38.ncbi.fasta",
    output:
        allsv=temp(
            "results/{cohort}/SV/sample_SVisionpro/{sample}.svision_pro_v2.1.s3.vcf"
        ),
        log=temp(
            "results/{cohort}/SV/sample_SVisionpro/{sample}.svision_pro_v2.1.s3.log"
        ),
    log:
        "logs/{cohort}/SV/sample_SVisionpro/{sample}.SVision.log",
    benchmark:
        "benchmarks/{cohort}/SV/sample_SVisionpro/{sample}.SVision.benchmark"
    container:
        f"docker://dockerproxy.net/caoboss/svisionpro:{config['SV']['SVisionpro']['version']}"
    threads: config["SV"]["SVisionpro"]["threads"]
    params:
        model_path=f'{config["SV"]["SVisionpro"]["model_dir"]}/model_liteunet_256_8_16_32_32_32.pth',
        access_path="resources/reference/SV/svisionpro.hg38.access.10M.bed",
        extra=config["SV"]["SVisionpro"]["extra_params"],
    shell:
        """
        (
            SVision-pro \
                --target_path {input.bam} \
                --genome_path {input.genome} \
                --model_path {params.model_path} \
                --process_num {threads} \
                --out_path results/{wildcards.cohort}/SV/sample_SVisionpro \
                --sample_name {wildcards.sample} \
                --access_path {params.access_path} \
                {params.extra}
        ) &> {log}
        """


rule separate_complex_SV:
    input:
        rules.SV_SVisionpro_call.output.allsv,
    output:
        simpleSV=temp("results/{cohort}/SV/sample_SVisionpro/{sample}.simple.vcf"),
        complexSV=temp("results/{cohort}/SV/sample_SVisionpro/{sample}.complex.vcf"),
    script:
        "separate_SVision_CSV.py"


rule gzip_SVisionpro:
    input:
        simpleSV=rules.separate_complex_SV.output.simpleSV,
        complexSV=rules.separate_complex_SV.output.complexSV,
    output:
        simpleSV="results/{cohort}/SV/sample_SVisionpro/{sample}.simple.vcf.gz",
        complexSV="results/{cohort}/SV/sample_SVisionpro/{sample}.complex.vcf.gz",
    conda:
        "../../envs/SV_process.yaml"
    shell:
        """
            bgzip -c {input.simpleSV} > {output.simpleSV}
            bgzip -c {input.complexSV} > {output.complexSV}
        """
