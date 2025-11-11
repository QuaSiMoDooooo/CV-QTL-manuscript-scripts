#!/usr/bin/env bash
set -euo pipefail

# -------------
# FileName     : 03_all_vartype_vcf_merge.sh
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Merge multiple VCF files (SNP, InDel, MNV, SV) into a single VCF for LD pattern analysis
# -------------

INPUT_DIR="02_vcf_format_simplified"
OUTPUT_VCF="03_all_vartype_merged.vcf"
TMP_DIR=".tmp_merge_vcf_order"
FILES=(
  "SNP.biochem_common.vcf"
  "InDel.biochem_common.vcf"
  "MNV.biochem_common.vcf"
  "SV.biochem_common.vcf"
)

mkdir -p "$TMP_DIR"

for f in "${FILES[@]}"; do
  bcftools query -l "$INPUT_DIR/$f" | sort > "$TMP_DIR/$f.samples"
done

cp "$TMP_DIR/${FILES[0]}.samples" "$TMP_DIR/common.samples"
for f in "${FILES[@]:1}"; do
  comm -12 "$TMP_DIR/common.samples" "$TMP_DIR/$f.samples" > "$TMP_DIR/_tmp.samples"
  mv "$TMP_DIR/_tmp.samples" "$TMP_DIR/common.samples"
done

COMMON="$TMP_DIR/common.samples"
N_COMMON=$(wc -l < "$COMMON")
if (( N_COMMON == 0 )); then
  echo "Error: No common samples found in all VCF files"
  exit 1
fi

GZ_FILES=()
for f in "${FILES[@]}"; do
  base=${f%.vcf}
  src="$INPUT_DIR/$f"
  reorder="$TMP_DIR/${base}.reorder.vcf"
  gz="$reorder.gz"

  bcftools view -S "$COMMON" "$src" -Ov -o "$reorder"
  bgzip -c "$reorder" > "$gz"
  tabix -p vcf "$gz"

  GZ_FILES+=("$gz")
done

bcftools concat -a "${GZ_FILES[@]}" -Ou | bcftools sort -Ov -o "$OUTPUT_VCF"

rm -rf "$TMP_DIR"

echo "Completed: $OUTPUT_VCF"