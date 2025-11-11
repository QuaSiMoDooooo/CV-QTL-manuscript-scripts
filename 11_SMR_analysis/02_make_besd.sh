# FileName     : 02_make_besd.sh
# Author       : Tian Wu fabien012390wt@163.com
# Modified By  : Tian Wu
# -------------
# Description  : Generate BESD files from flist files using SMR tool
# -------------

input_dir="01_flist"
output_base_dir="02_besd"

for input_flist in "$input_dir"/*.flist; do
    flist_name=$(basename "$input_flist" .flist)
    output_dir="$output_base_dir/$flist_name"
    
    if [ ! -d "$output_dir" ]; then
        mkdir -p "$output_dir"
    fi
    
    output_prefix="$output_dir/$flist_name"
    
    if [ ! -f "$output_prefix.besd" ]; then
        echo "Processing $input_flist"
        echo "Output prefix: $output_prefix"
        smr --eqtl-flist "$input_flist" --make-besd --cis-wind 2000 --trans-wind 1000 \
            --peqtl-trans 5.0e-8 --peqtl-other 1.0e-5 --out "$output_prefix" --thread-num 36
    fi
done