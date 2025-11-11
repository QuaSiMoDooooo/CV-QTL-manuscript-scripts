#!/bin/bash

# -------------
# FileName     : 06_get_promoter_fa.sh
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Extract promoter region sequences from reference genome using bedtools
# -------------

BED=$1
FASTA=$2
OUTFILE=$3

> "$OUTFILE"

while read -r chr start end; do
    seq_name="${chr}_${start}_${end}"
    
    echo -e "${chr}\t${start}\t${end}" | \
        bedtools getfasta -fi "$FASTA" -bed - -fo - -name | \
        awk -v name="$seq_name" 'BEGIN{OFS="\n"} {if(NR%2==1){print ">"name}else{print $0}}' >> "$OUTFILE"
done < "$BED"
