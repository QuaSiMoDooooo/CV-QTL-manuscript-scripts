# -*- coding: UTF-8 -*-
#
# FileName     : SV_identify_SVision
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-10-08 21:05
# Last Modified: 2024-10-09 11:24
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule SV_SVision_call:
    input:
        bam=rules.alignment_merge.output.bam,
        reference="resources/reference/alignment/GRCh38.ncbi.fasta",
    output:
        vcf=temp("results/{cohort}/SV/sample_SVision/{sample}/{sample}.svision.s3.graph.vcf"),
        graph_exactly=temp("results/{cohort}/SV/sample_SVision/{sample}/{sample}.graph_exactly_match.txt"),
        graph_symmetry=temp("results/{cohort}/SV/sample_SVision/{sample}/{sample}.graph_symmetry_match.txt"),
    log:
        "logs/{cohort}/SV/sample_SVision/{sample}.SVision_call.log",
    benchmark:
        "benchmarks/{cohort}/SV/sample_SVision/{sample}.SVision_call.benchmark"
    container:
        f"docker://dockerproxy.net/jiadongxjtu/svision:{config['SV']['SVision']['SVision_version']}"
    threads: config["SV"]["SVision"]["threads"]
    params:
        model_dir=f'{config["SV"]["SVision"]["model_dir"]}/svision-cnn-model.ckpt',
        extra=config["SV"]["SVision"]["extra_params"],
    shell:
        """
        (
            SVision \
                -b {input.bam} \
                -g {input.reference} \
                -o results/{wildcards.cohort}/SV/sample_SVision/{wildcards.sample} \
                -n {wildcards.sample} \
                -m {params.model_dir} \
                {params.extra}
        ) &> {log}
        """


rule separate_complex_SV:
    input:
        rules.SV_SVision_call.output.vcf,
    output:
        formalSV=temp("results/{cohort}/SV/sample_SVision/{sample}/{sample}.svision.formal.vcf"),
        complexSV=temp("results/{cohort}/SV/sample_SVision/{sample}/{sample}.svision.complex.vcf"),
    script:
        "separate_SVision_CSV.py"


rule gzip_SVision:
    input:
        formalSV=rules.separate_complex_SV.output.formalSV,
        complexSV=rules.separate_complex_SV.output.complexSV,
    output:
        formalSV="results/{cohort}/SV/sample_SVision/{sample}/{sample}.formal.vcf.gz",
        complexSV="results/{cohort}/SV/sample_SVision/{sample}/{sample}.complex.vcf.gz",
    conda:
        "SV"
    shell:
        """
            bgzip -c {input.formalSV} > {output.formalSV}
            bgzip -c {input.complexSV} > {output.complexSV}
        """
