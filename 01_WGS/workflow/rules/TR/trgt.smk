# -*- coding: UTF-8 -*-
#
# FileName     : trgt
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-18 09:44
# Last Modified: 2024-11-18 10:38
# Modified By  : EastsunW
# -------------
# Description  : https://www.nature.com/articles/s41587-023-02057-3
# -------------


rule trgt_identify_normal:
    input:
        genome="resources/reference/alignment/GRCh38.ncbi.fasta",
        repeats=cohort.ref_config["TR"]["normal"],
        bam=rules.alignment_merge.output.bam,
    output:
        trft_bam="results/{cohort}/TR/sample_trgt/{sample}/{sample}.normal.spanning.bam",
        trft_bam_index="results/{cohort}/TR/sample_trgt/{sample}/{sample}.normal.spanning.bam.bai",
        trgt_vcf="results/{cohort}/TR/sample_trgt/{sample}/{sample}.normal.vcf.gz",
        trgt_vcf_index="results/{cohort}/TR/sample_trgt/{sample}/{sample}.normal.vcf.gz.tbi",
    log:
        "logs/{cohort}/TR/sample_trgt/{sample}.normal.trgt.log",
    benchmark:
        "benchmarks/{cohort}/TR/sample_trgt/{sample}.normal.identity.benchmark"
    threads: config["TR"]["trgt"]["threads"]
    container:
        f"docker://quay.dockerproxy.net/pacbio/trgt:{config['TR']['trgt']['version']}"
    params:
        karyotype=lambda wildcards: "XX" if wildcards.sample[2] == "2" else "XY",
        extra=config["TR"]["trgt"]["extra_params"],
    shell:
        """
        (
            trgt genotype \
                --threads {threads} \
                --repeats {input.repeats} \
                --genome {input.genome} \
                --reads {input.bam} \
                --karyotype {params.karyotype} \
                --sample-name {wildcards.sample} \
                --output-prefix results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.raw \
                {params.extra}
            
            samtools sort \
                --threads {threads} \
                -o results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.spanning.bam \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.raw.spanning.bam
            
            samtools index \
                --threads {threads} \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.spanning.bam
            
            bcftools sort \
                --output-type z \
                --output results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.vcf.gz \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.raw.vcf.gz
            
            bcftools index \
                --threads {threads} \
                --tbi \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.vcf.gz
            rm results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.normal.raw.*
        ) &> {log}
        """


rule trgt_identify_disease:
    input:
        genome="resources/reference/alignment/GRCh38.ncbi.fasta",
        repeats=cohort.ref_config["TR"]["disease"],
        bam=rules.alignment_merge.output.bam,
    output:
        trft_bam="results/{cohort}/TR/sample_trgt/{sample}/{sample}.disease.spanning.bam",
        trft_bam_index="results/{cohort}/TR/sample_trgt/{sample}/{sample}.disease.spanning.bam.bai",
        trgt_vcf="results/{cohort}/TR/sample_trgt/{sample}/{sample}.disease.vcf.gz",
        trgt_vcf_index="results/{cohort}/TR/sample_trgt/{sample}/{sample}.disease.vcf.gz.tbi",
    log:
        "logs/{cohort}/TR/sample_trgt/{sample}.disease.trgt.log",
    benchmark:
        "benchmarks/{cohort}/TR/sample_trgt/{sample}.disease.identity.benchmark"
    threads: config["TR"]["trgt"]["threads"]
    container:
        f"docker://quay.dockerproxy.net/pacbio/trgt:{config['TR']['trgt']['version']}"
    params:
        karyotype=lambda wildcards: "XX" if wildcards.sample[2] == "2" else "XY",
        extra=config["TR"]["trgt"]["extra_params"],
    shell:
        """
        (
            trgt genotype \
                --threads {threads} \
                --repeats {input.repeats} \
                --genome {input.genome} \
                --reads {input.bam} \
                --karyotype {params.karyotype} \
                --sample-name {wildcards.sample} \
                --output-prefix results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.raw \
                {params.extra}
            
            samtools sort \
                --threads {threads} \
                -o results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.spanning.bam \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.raw.spanning.bam
            
            samtools index \
                --threads {threads} \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.spanning.bam
            
            bcftools sort \
                --output-type z \
                --output results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.vcf.gz \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.raw.vcf.gz
            
            bcftools index \
                --threads {threads} \
                --tbi \
                results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.vcf.gz
            rm results/{wildcards.cohort}/TR/sample_trgt/{wildcards.sample}/{wildcards.sample}.disease.raw.*
        ) &> {log}
        """



rule trgt_merge:
    input:
        genome="resources/reference/alignment/GRCh38.ncbi.fasta",
        vcf=lambda wildcards: expand(
            "results/{{cohort}}/TR/sample_trgt/{sample}/{sample}.{{tr_type}}.vcf.gz",
            sample=cohort.list_samples(),
            allow_missing=False,
        ),
    output:
        vcf="results/{cohort}/TR/cohort_trgt/{cohort}.{tr_type}.vcf.gz",
        index="results/{cohort}/TR/cohort_trgt/{cohort}.{tr_type}.vcf.gz.tbi",
    log:
        "logs/{cohort}/TR/cohort_trgt/{cohort}.{tr_type}.merge.log",
    benchmark:
        "benchmarks/{cohort}/TR/cohort_trgt/{cohort}.{tr_type}.merge.benchmark"
    container:
        f"docker://quay.dockerproxy.net/pacbio/trgt:{config['TR']['trgt']['version']}"
    shell:
        """
        trgt merge \
            --vcf {input.vcf} \
            --genome {input.genome} \
            --output {output.vcf}
        bcftools index --tbi {output.vcf}
        """
