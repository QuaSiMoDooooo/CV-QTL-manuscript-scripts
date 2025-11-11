# -*- coding: UTF-8 -*-
#
# FileName     : pepper_deepvariant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-02-24 17:20
# Last Modified: 2024-02-24 17:20
# Modified By  : EastsunW
# -------------
# Description  : 使用pepper-Deepvariant进行SNP/Indel的鉴定
# Citation     ：https://doi.org/10.1038/s41592-021-01299-w
# -------------


rule SNP_identify_pepper:
    input:
        bam=rules.alignment_merge.output.bam,
        reference=rules.alignment_movie.input.reference,
    output:
        vcf="results/{cohort}/SNP/pepper_sample/{sample}/{sample}.vcf.gz",
        vcf_index="results/{cohort}/SNP/pepper_sample/{sample}/{sample}.vcf.gz.tbi",
        gvcf="results/{cohort}/SNP/pepper_sample/{sample}/{sample}.g.vcf.gz",
        gvcf_index="results/{cohort}/SNP/pepper_sample/{sample}/{sample}.g.vcf.gz.tbi",
    retries: 5
    log:
        "logs/{cohort}/SNP/pepper_sample/{sample}.snp_identify.log",
    benchmark:
        "benchmarks/{cohort}/SNP/pepper_sample/{sample}.snp_identify.benchmark"
    container:
        f"docker://kishwars/pepper_deepvariant:{config['SNP']['pepper']['pepper_version']}"
    threads: config["SNP"]["pepper"]["pepper_threads"]
    params:
        extra=config["SNP"]["pepper"]["pepper_params"],
    shell:
        """
        (
            run_pepper_margin_deepvariant call_variant \
                --bam {input.bam} \
                --fasta {input.reference} \
                --output_dir results/{wildcards.cohort}/SNP/pepper_sample/{wildcards.sample} \
                --output_prefix {wildcards.sample} \
                --threads {threads} \
                --sample_name {wildcards.sample} \
                {params.extra}
        ) &> {log}

            bcftools view \
                --threads {threads} \
                -i 'FILTER!="NoCall"' \
                {output.gvcf} \
                -Oz5 -o results/{wildcards.cohort}/SNP/pepper_sample/{wildcards.sample}/{wildcards.sample}_filtered.g.vcf.gz
            bcftools index -t --threads {threads} results/{wildcards.cohort}/SNP/pepper_sample/{wildcards.sample}/{wildcards.sample}_filtered.g.vcf.gz
            rm {output.gvcf} {output.gvcf_index}
            mv results/{wildcards.cohort}/SNP/pepper_sample/{wildcards.sample}/{wildcards.sample}_filtered.g.vcf.gz {output.gvcf}
            mv results/{wildcards.cohort}/SNP/pepper_sample/{wildcards.sample}/{wildcards.sample}_filtered.g.vcf.gz.tbi {output.gvcf_index}
        """


rule SNP_merge_pepper:
    input:
        expand(
            "results/{{cohort}}/SNP/pepper_sample/{sample}/{sample}.g.vcf.gz",
            sample=cohort.list_samples(),
            allow_missing=False,
        ),
    output:
        scratch_dir=temp(directory("results/{cohort}/SNP/pepper_cohort/GLnexus.DB")),
        merged_BCF=temp("results/{cohort}/SNP/pepper_cohort/{cohort}.cohort.bcf"),
        merged_VCF="results/{cohort}/SNP/pepper_cohort/{cohort}.cohort.vcf.gz",
        merged_VCF_index="results/{cohort}/SNP/pepper_cohort/{cohort}.cohort.vcf.gz.csi",
    log:
        "logs/{cohort}/SNP/pepper_cohort/{cohort}.merge.log",
    benchmark:
        "benchmarks/{cohort}/SNP/pepper_cohort/{cohort}.merge.benchmark"
    container:
        f"docker://ghcr.io/dnanexus-rnd/glnexus:{config['SNP']['pepper']['glnexus_version']}"
    threads: config["SNP"]["pepper"]["glnexus_threads"]
    params:
        add_tags=f"-- -t {config['SNP']['pepper']['glnexus_add_tags']}"
        if config["SNP"]["pepper"]["glnexus_add_tags"]
        else "",
        extra=config["SNP"]["pepper"]["glnexus_params"],
    shell:
        """
        (
            glnexus_cli \
                --threads {threads} \
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
