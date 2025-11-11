#! /bin/bash
# FileName     : stat_sample_variant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2025-03-03 19:19
# Last Modified: 2025-03-03 19:44
# Modified By  : EastsunW
# -------------
# Description  : 
# -------------

set -euo pipefail

# 定义目录路径
vcf_directory="/home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SNP/clair3_sample"  # 替换为你的 VCF 文件所在目录
out_directory="/home/wangdy/Projects/Weibin/Downstreams/data_stat/results/stat/Variant_stat"  # 替换为你的输出目录

# 检查目录是否存在
if [ ! -d "$vcf_directory" ]; then
    echo "Error: Directory $vcf_directory does not exist."
    exit 1
fi

# 检查输出目录是否存在，不存在则创建
if [ ! -d "$out_directory" ]; then
    mkdir -p "$out_directory"
fi

# 创建一个输出文件用于保存统计结果
output_file=${out_directory}/vcf_variant_counts.txt
echo -e "Sample\tSNP\tInDel" > "$output_file"

# 遍历目录中的所有 VCF 文件
find "$vcf_directory" -name "*.vcf.gz" | while read -r vcf_file; do
    if [ -f "$vcf_file" ]; then
        echo "Processing $vcf_file..."
        # 提取文件名
        file_name=$(basename "$vcf_file")
        # 文件名去掉.vcf.gz就是样本名
        sample_name=${file_name%.vcf.gz}
        # 统计 SNP 和 InDel 数量
        snp_count=$(vcftools --gzvcf "$vcf_file" --remove-indels --recode --recode-INFO-all --stdout | grep -v "^#" | wc -l)
        indel_count=$(vcftools --gzvcf "$vcf_file" --keep-only-indels --recode --recode-INFO-all --stdout | grep -v "^#" | wc -l)
        echo "$snp_count"
        # 将结果写入输出文件
        echo -e "$sample_name\t$snp_count\t$indel_count" >> "$output_file"
    fi
done

echo "All files processed. Results saved to $output_file."
