# -------------
# FileName     : 01_merge_bed_files.sh
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Merge multiple BED files and extract genomic regions for colocalization analysis
# -------------

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 file1.bed file2.bed ..."
  exit 1
fi

first_file=$(basename "$1")
tag=$(echo "$first_file" | awk -F'_' '{print $2}')

if [ -z "$tag" ]; then
  echo "Error: Failed to extract tag from filename '$first_file'"
  exit 1
fi

OUTPUT_FILE="merged_beds_${tag}.bed"
TMP_FILE=$(mktemp)

cat "$@" > "$TMP_FILE"
sort -k1,1n -k2,2n "$TMP_FILE" | bedtools merge > "$OUTPUT_FILE"
rm "$TMP_FILE"

echo "Merged BED saved to $OUTPUT_FILE"