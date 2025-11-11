# -*- coding: UTF-8 -*-
#
# FileName     : prepare_genotype
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-09-20 09:50
# Last Modified: 2024-09-20 10:20
# Modified By  : EastsunW
# -------------
# Description  : 从原始的变异中提取出基因型信息,基因型中的变异是合并后的ID
# -------------


rule separate_snp_indel_and_filter:
    input:
        "resources/variants/snp_indel.raw.vcf.gz",
    output:
        snp="resources/variants/SNP.filtered.vcf.gz",
        indel="resources/variants/InDel.filtered.vcf.gz",
    params:
        hwe=config["data_filter"]["SNP"]["min_hwe"],
        maf=config["data_filter"]["SNP"]["min_maf"],
        na=config["data_filter"]["SNP"]["max_na_rate"],
        exclude_region=config["data_filter"]["blacklist"],
    conda:
        "../envs/variant.yaml"
    shell:
        """
        vcftools --gzvcf {input} \
            --not-chr chrY --not-chr chrM \
            --hwe {params.hwe} --maf {params.maf} --max-missing {params.na} \
            --exclude-bed {params.exclude_region} \
            --min-alleles 2 --max-alleles 2 \
            --remove-indels \
            --recode --recode-INFO-all \
            --stdout | gzip -c > {output.snp}
        vcftools --gzvcf {input} \
            --not-chr chrY --not-chr chrM \
            --hwe {params.hwe} --maf {params.maf} --max-missing {params.na} \
            --exclude-bed {params.exclude_region} \
            --min-alleles 2 --max-alleles 2 \
            --keep-only-indels \
            --recode --recode-INFO-all \
            --stdout | gzip -c > {output.indel}
        """


rule variant_genotype_SNP:
    input:
        rules.separate_snp_indel_and_filter.output.snp,
    output:
        "results/common_data/variant/SNP.genotype.filtered.txt",
    shell:
        """
        python scripts/data_prepare/extract_gt_from_vcf.py SNP -i {input} -o {output}
        """


rule variant_genotype_InDel:
    input:
        rules.separate_snp_indel_and_filter.output.indel,
    output:
        "results/common_data/variant/InDel.genotype.filtered.txt",
    shell:
        """
        python scripts/data_prepare/extract_gt_from_vcf.py SNP -i {input} -o {output}
        """


rule variant_genotype_MNV_filter:
    input:
        mnv="resources/variants/MNV.raw.vcf.gz",
        snpvcf=rules.variant_genotype_SNP.input,
    output:
        "results/common_data/variant/MNV.genotype.filtered.txt",
    params:
        maf=(
            f'--maf {config["data_filter"]["MNV"]["min_maf"]}'
            if config["data_filter"]["MNV"]["min_maf"]
            else ""
        ),
        exclude_YM="--exclude_YM",
        exclude_region=f'--exclude_region {config["data_filter"]["blacklist"]}',
    shell:
        """
        python scripts/data_prepare/extract_gt_from_vcf.py MNV \
            -i {input.mnv} \
            -v {input.snpvcf} \
            -o {output} \
            {params.maf} \
            {params.exclude_YM} \
            {params.exclude_region}
        """


rule variant_genotype_SV_filter:
    input:
        "resources/variants/SV.raw.vcf.gz",
    output:
        "results/common_data/variant/SV.genotype.filtered.txt",
    params:
        maf=(
            f'--maf {config["data_filter"]["SV"]["min_maf"]}'
            if config["data_filter"]["SV"]["min_maf"]
            else ""
        ),
        na=(
            f'--na {config["data_filter"]["SV"]["max_na_rate"]}'
            if config["data_filter"]["SV"]["max_na_rate"]
            else ""
        ),
        exclude_YM="--exclude_YM",
        exclude_region=f'--exclude_region {config["data_filter"]["blacklist"]}',
    shell:
        """
        python scripts/data_prepare/extract_gt_from_vcf.py SV \
            -i {input} \
            -o {output} \
            {params.maf} \
            {params.na} \
            {params.exclude_YM} \
            {params.exclude_region}
        """


rule variant_position:
    input:
        "results/common_data/variant/{variant}.genotype.filtered.txt",
    output:
        "results/common_data/variant/{variant}.position.txt",
    shell:
        """
        python scripts/data_prepare/extract_variant_position.py {wildcards.variant} \
            -i {input} \
            -o {output}
        """


for phenotype in ["expression", "APA", "splicing", "methylation"]:

    rule:
        name:
            f"copy_variant_files_{phenotype}"
        input:
            genotype="results/common_data/variant/{variant}.genotype.filtered.txt",
            position=rules.variant_position.output,
        output:
            target_genotype=f"results/{{variant}}-{phenotype}/filtered_data/variant.genotype.filtered.txt",
            target_position=f"results/{{variant}}-{phenotype}/matched_data/variant.position.txt",
        shell:
            """
            cp {input.genotype} {output.target_genotype}
            cp {input.position} {output.target_position}
            """
