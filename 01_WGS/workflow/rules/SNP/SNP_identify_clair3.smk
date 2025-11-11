# -*- coding: UTF-8 -*-
#
# FileName     : SNP_identify_clair3
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-23 21:35
# Last Modified: 2024-04-23 21:35
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule SNP_identify_clair3:
    input:
        bam=rules.alignment_merge.output.bam,
        reference="resources/reference/alignment/GRCh38.ncbi.fasta",
    output:
        vcf="results/{cohort}/SNP/clair3_sample/{sample}.vcf.gz",
        vcf_index="results/{cohort}/SNP/clair3_sample/{sample}.vcf.gz.tbi",
        gvcf="results/{cohort}/SNP/clair3_sample/{sample}.gvcf.gz",
        gvcf_index="results/{cohort}/SNP/clair3_sample/{sample}.gvcf.gz.tbi",
        temp_dir=temp(directory("results/{cohort}/SNP/clair3_sample/{sample}")),
    retries: 3
    log:
        "logs/{cohort}/SNP/clair3_sample/{sample}.clair3_identify.log",
    benchmark:
        "benchmarks/{cohort}/SNP/clair3_sample/{sample}.clair3_identify.benchmark"
    container:
        f"docker://dockerproxy.net/hkubal/clair3:{config['SNP']['clair3']['clair3_version']}"
    threads: config["SNP"]["clair3"]["clair3_threads"]
    params:
        model=config["SNP"]["clair3"]["clair3_model"],
        extra=config["SNP"]["clair3"]["clair3_params"],
    shell:
        """
        (
            /opt/bin/run_clair3.sh \
                --bam_fn={input.bam} \
                --ref_fn={input.reference} \
                --threads={threads} \
                --model_path="/opt/models/{params.model}" \
                --sample_name={wildcards.sample} \
                --output=results/{wildcards.cohort}/SNP/clair3_sample/{wildcards.sample} \
                {params.extra} && \
            mv results/{wildcards.cohort}/SNP/clair3_sample/{wildcards.sample}/merge_output.vcf.gz {output.vcf} && \
            mv results/{wildcards.cohort}/SNP/clair3_sample/{wildcards.sample}/merge_output.vcf.gz.tbi {output.vcf_index} && \
            mv results/{wildcards.cohort}/SNP/clair3_sample/{wildcards.sample}/merge_output.gvcf.gz {output.gvcf} && \
            mv results/{wildcards.cohort}/SNP/clair3_sample/{wildcards.sample}/merge_output.gvcf.gz.tbi {output.gvcf_index}
        ) &> {log}
        """


rule SNP_merge_clair3:
    input:
        expand(
            "results/{{cohort}}/SNP/clair3_sample/{sample}.gvcf.gz",
            sample=cohort.list_samples(),
            allow_missing=False,
        ),
    output:
        scratch_dir=temp(directory("results/{cohort}/SNP/clair3_cohort/GLnexus.DB")),
        merged_BCF=temp("results/{cohort}/SNP/clair3_cohort/{cohort}.cohort.bcf"),
        merged_VCF="results/{cohort}/SNP/clair3_cohort/{cohort}.cohort.vcf.gz",
        merged_VCF_index="results/{cohort}/SNP/clair3_cohort/{cohort}.cohort.vcf.gz.csi",
    log:
        "logs/{cohort}/SNP/clair3_cohort/{cohort}.merge.log",
    benchmark:
        "benchmarks/{cohort}/SNP/clair3_cohort/{cohort}.merge.benchmark"
    container:
        f"docker://ghcr.dockerproxy.net/dnanexus-rnd/glnexus:{config['SNP']['clair3']['glnexus_version']}"
    threads: config["SNP"]["clair3"]["glnexus_threads"]
    params:
        add_tags=(
            f"-- -t {config['SNP']['clair3']['glnexus_add_tags']}"
            if config["SNP"]["clair3"]["glnexus_add_tags"]
            else ""
        ),
        config=cohort.ref_config["SNP"]["clair3"]["glnexus_config"],
        extra=config["SNP"]["clair3"]["glnexus_params"],
    shell:
        """
        (
            glnexus_cli \
                --threads {threads} \
                --config {params.config} \
                --dir {output.scratch_dir} \
                {params.extra} \
                {input} \
            > {output.merged_BCF}
            bcftools view {output.merged_BCF} | \
                bcftools +fill-tags -Oz5 {params.add_tags} | \
                bcftools view -Oz5 > \
                {output.merged_VCF}
            bcftools index {output.merged_VCF}
        ) &> {log}
        """


rule SNP_merge_clair3_step_stat:
    input:
        expand(
            "results/{{cohort}}/SNP/clair3_sample/{sample}.vcf.gz",
            sample=cohort.list_samples(),
        ),
    output:
        "results/{cohort}/SNP/SNP_sample_increase.stat.txt",
    script:
        "stat_variant_sample_increase.py"
