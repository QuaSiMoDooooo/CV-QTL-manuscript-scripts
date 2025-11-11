#! /bin/bash
# FileName     : download_SV_data
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-12-17 21:21
# Last Modified: 2025-02-27 20:35
# Modified By  : EastsunW
# -------------
# Description  :
# -------------

set -euo pipefail

download_1KG() {
    local target_dir="data/SV_compare/1KG"
    mkdir -p "$target_dir"
    local chromosomes=($(seq 1 22))
    for chr in "${chromosomes[@]}"; do
        axel \
            -c \
            -n 10 \
            --quiet \
            "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz" \
            -o "${target_dir}/1KG_chr${chr}.alltype.vcf.gz"
        axel \
            -c \
            -n 10 \
            --quiet \
            "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi" \
            -o "${target_dir}/1KG_chr${chr}.alltype.vcf.gz.tbi"
    done
    axel \
        -c \
            -n 10 \
            --quiet \
            "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz" \
            -o "${target_dir}/1KG_chrX.alltype.vcf.gz"
    axel \
        -c \
            -n 10 \
            --quiet \
            "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi" \
            -o "${target_dir}/1KG_chrX.alltype.vcf.gz.tbi"
}

download_gnomAD() {
    local target_dir="data/SV_compare/gnomAD"
    mkdir -p "$target_dir"
    axel \
        -c \
        -n 10 \
        --quiet \
        "https://datasetgnomad.blob.core.windows.net/dataset/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz" \
        -o "${target_dir}/gnomAD/gnomAD.sv.vcf.gz"
}



download_Iceland2021() {
    local target_dir="data/SV_compare"
    mkdir -p "$target_dir"
    axel \
        -c \
        -n 10 \
        "https://github.com/DecodeGenetics/LRS_SV_sets/raw/refs/heads/master/ont_sv_high_confidence_SVs.sorted.vcf.gz" \
        -o "$target_dir/Iceland2021.vcf.gz"
    axel \
        -c \
        -n 10 \
        "https://github.com/DecodeGenetics/LRS_SV_sets/raw/refs/heads/master/ont_sv_high_confidence_SVs.sorted.vcf.gz.tbi" \
        -o "$target_dir/Iceland2021.vcf.gz.tbi"
    axel \
        -c \
        -n 10 \
        "https://github.com/DecodeGenetics/LRS_SV_sets/raw/refs/heads/master/ont_sv_high_confidence_tandemdup.csv" \
        -o "$target_dir/Iceland2021.tandemdup.csv"
}

download_HGDP() {
    local target_dir="data/SV_compare"
    mkdir -p "$target_dir/HGDP"
    axel \
        -c \
        -n 10 \
        --quiet \
        "https://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/del.manta.vcf.gz" \
        -o "${target_dir}/HGDP/HGDP.DEL.vcf.gz"
    axel \
        -c \
        -n 10 \
        --quiet \
        "https://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/dup.manta.vcf.gz" \
        -o "${target_dir}/HGDP/HGDP.DUP.vcf.gz"
    axel \
        -c \
        -n 10 \
        --quiet \
        "https://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/ins.manta.vcf.gz" \
        -o "${target_dir}/HGDP/HGDP.INS.vcf.gz"
    axel \
        -c \
        -n 10 \
        --quiet \
        "https://ngs.sanger.ac.uk/production/hgdp/hgdp_structural_variation/SV_CALLSET/inv.manta.vcf.gz" \
        -o "${target_dir}/HGDP/HGDP.INV.vcf.gz"
}



# 入口
# download_1KG
# download_gnomAD
# download_Iceland2021
download_HGDP
