# -*- coding: UTF-8 -*-
#
# FileName     : MNV_identify
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-04-08 22:23
# Last Modified: 2024-04-08 22:23
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule MNV_clair3_phase_split_chr:
    input:
        rules.SNP_merge_clair3.output.merged_VCF,
    output:
        vcf=temp("results/{cohort}/MNV/clair3/{cohort}.vcf.chr{chr}.gz"),
        index=temp("results/{cohort}/MNV/clair3/{cohort}.vcf.chr{chr}.gz.csi"),
    threads: 1
    shell:
        """
        bcftools view \
            -r chr{wildcards.chr} \
            -O z5 \
            -o {output.vcf} \
            --write-index \
            {input}
        """


rule MNV_clair3_phase_chr:
    input:
        vcf=rules.MNV_clair3_phase_split_chr.output.vcf,
        index=rules.MNV_clair3_phase_split_chr.output.index,
    output:
        bcf=temp("results/{cohort}/MNV/clair3/{cohort}.phased.chr{chr}.bcf"),
        bcf_index=temp("results/{cohort}/MNV/clair3/{cohort}.phased.chr{chr}.bcf.csi"),
        vcf=temp("results/{cohort}/MNV/clair3/{cohort}.phased.chr{chr}.vcf.gz"),
        vcf_index=temp(
            "results/{cohort}/MNV/clair3/{cohort}.phased.chr{chr}.vcf.gz.csi"
        ),
    log:
        "logs/{cohort}/MNV/clair3/phase_chr/{cohort}.phase_chr{chr}.log",
    benchmark:
        "benchmarks/{cohort}/MNV/clair3/phase_chr/{cohort}.phase_chr{chr}.benchmark"
    conda:
        "../../envs/MNV.yaml"
    threads: config["MNV"]["shapeit"]["threads"]
    params:
        ref_panel=lambda wildcards: f"{cohort.ref_config['MNV']['shapeit']['reference_panel']}".replace(
            "%chr%", f"{wildcards.chr}"
        ),
        extra=config["MNV"]["shapeit"]["extra_params"],
    shell:
        """
        (
            SHAPEIT5_phase_common \
                --thread {threads} \
                --reference {params.ref_panel} \
                --input {input.vcf} \
                --region chr{wildcards.chr} \
                {params.extra} \
                --output {output.bcf}
            bcftools view \
                --write-index \
                -O z5 \
                -o {output.vcf} \
                {output.bcf}
        ) &> {log}
        """


rule MNV_clair3_phase_merge:
    input:
        lambda wildcards: expand(
            "results/{{cohort}}/MNV/clair3/{{cohort}}.phased.chr{chr}.vcf.gz",
            chr=[str(i) for i in range(1, 23)] + ["X"],
            allow_missing=False,
        ),
    output:
        vcf="results/{cohort}/MNV/clair3/{cohort}.phased.vcf.gz",
        index="results/{cohort}/MNV/clair3/{cohort}.phased.vcf.gz.csi",
    conda:
        "../../envs/MNV.yaml"
    threads: 10
    shell:
        """
            bcftools concat \
                --threads {threads} \
                -O z5 \
                --write-index \
                -o {output.vcf} \
                {input}
        """


rule MNV_clair3_identify_process_GT:
    input:
        rules.MNV_clair3_phase_merge.output.vcf,
    output:
        temp("results/{cohort}/MNV/clair3/{cohort}.phased.gt.vcf"),
    threads: config["MNV"]["process_GT"]["threads"]
    shell:
        """
            bcftools view -h {input} > \
                {output}
            bcftools query \
                -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT[\t%GT]\n' \
                {input} >> \
                {output}
        """


rule MNV_clair3_identify:
    input:
        rules.MNV_clair3_identify_process_GT.output,
    output:
        mnv="results/{cohort}/MNV/clair3/{cohort}.mnv.txt",
        mnv_multi="results/{cohort}/MNV/clair3/{cohort}.mnv_multi.txt",
    log:
        "logs/{cohort}/MNV/clair3/{cohort}.MNV_identify.log",
    benchmark:
        "benchmarks/{cohort}/MNV/clair3/{cohort}.MNV_identify.benchmark"
    conda:
        "../../envs/MNV.yaml"
    params:
        extra=config["MNV"]["MNVidentify"]["extra_params"],
    shell:
        """
        (
            bash workflow/rules/MNV/MNV_identify_tools/MNVIdentify.sh \
                -i {input} \
                {params.extra}
            if [ -f results/{wildcards.cohort}/MNV/clair3/{wildcards.cohort}.phased.gt/mnv.txt ]; then
                mv results/{wildcards.cohort}/MNV/clair3/{wildcards.cohort}.phased.gt/mnv.txt {output.mnv}
            else
                touch {output.mnv}
            fi

            if [ -f results/{wildcards.cohort}/MNV/clair3/{wildcards.cohort}.phased.gt/mnv_multi.txt ]; then
                mv results/{wildcards.cohort}/MNV/clair3/{wildcards.cohort}.phased.gt/mnv_multi.txt {output.mnv_multi}
            else
                touch {output.mnv_multi}
            fi
            rm -rf results/{wildcards.cohort}/MNV/clair3/{wildcards.cohort}.phased.gt
        ) &> {log}
        """


rule give_MNV_ID:
    input:
        rules.MNV_clair3_identify.output.mnv,
    output:
        "results/{cohort}/MNV/clair3/{cohort}.mnv.withID.txt",
    script:
        "generate_MNV_ID.py"
