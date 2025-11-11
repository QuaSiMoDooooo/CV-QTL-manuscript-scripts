#! Rscript
# -------------
# FileName     : 01_plink2vcf.sh
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Convert PLINK binary format files (.bed/.bim/.fam) to VCF format for LD pattern analysis
# -------------

# Input directory path prefix
INPUT_DIR="/home/wtian/project/HZAU_cohort_meth/wdy_result/phenotype_qtl_results_bim3_removaNA_cp"
# Output directory
OUTPUT_DIR="./01_flt_vcf"
mkdir -p "$OUTPUT_DIR"

# Process all biochem_common.bim files
for BIM_FILE in $INPUT_DIR/*biochem*/*.biochem_common.bim; do
    PREFIX="${BIM_FILE%.bim}"
    BED_FILE="${PREFIX}.bed"
    FAM_FILE="${PREFIX}.fam"
    BASENAME=$(basename "$PREFIX")
    OUT_VCF="${OUTPUT_DIR}/${BASENAME}.vcf"

    if [[ -f "$BED_FILE" && -f "$FAM_FILE" ]]; then
        echo "Processing: $BASENAME"
        plink --bfile "$PREFIX" --recode vcf --keep-allele-order --out "${OUTPUT_DIR}/${BASENAME}"
    else
        echo "Missing required files: $BED_FILE or $FAM_FILE, skipping $PREFIX"
    fi
done

find "$OUTPUT_DIR" -type f ! -name "*.vcf" -delete

for vcf in "$OUTPUT_DIR"/*.vcf; do
    tr -d '\000' < "$vcf" > "${vcf}.clean"
    mv "${vcf}.clean" "$vcf"
done