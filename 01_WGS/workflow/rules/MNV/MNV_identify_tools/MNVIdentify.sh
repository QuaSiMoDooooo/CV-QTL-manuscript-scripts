#!/bin/bash

################  Description  ################
##         Name: MNVAnno - Systematic Identification and Annotation of MNVs
##      Authors: Weiwei Jin, Bioinfomatics Laboratory, HZAU, 2576720071@qq.com
##      Created: Feb 1, 2023 for creation
###############################################


################  Requirements  ################
##  (2) Python >= 3 
################################################

################  Quick Start  ################
## Test example: bash [scriptPath]/MNVIdentify.sh -i [scriptPath]/example/demo.vcf
## Parameters: bash [scriptPath]/MNVIdentify.sh -i [file url] -d [identification distance] -a [cutoff of AC] -m [cutoff of MAF] -f [switch] -c [chr list]
###############################################


############  Default Parameters  #################
    #Parameter of script path, used for both relative path and absolute path
    MNVIdentifyfolder=$0
    #Current_path record
    Current_path=($(pwd))
    #Auto detection and change the relative path into absolute path
    tmp=${MNVIdentifyfolder:0:${#MNVIdentifyfolder}-14}
    if [ x${tmp} != "x" ]
    then
    MNVIdentifyfolder=($(cd ${tmp} ; pwd)/)
    else
    MNVIdentifyfolder=($(pwd)/)
    fi
    #The absolute path of the VCF file
    file_url="0";
    mnv_distance=10;
    mnv_ac=0;
    mnv_maf=0;
    mnv_filter_switch='T';
    mnv_sample_index='F';
    multi_cutoff=6;
    chr_list='0';
    identifyMulti='F';
    delCache='T';
###################################################


###############  Reading user defined parameters  ###################
##Read in
    while getopts 'i:d:a:m:f:s:j:M:D:c:' optname
    do
        case $optname in
        i)
            file_url=$OPTARG;;
        d)
            mnv_distance=$OPTARG;;
        a)
            mnv_ac=$OPTARG;;
        m)
            mnv_maf=$OPTARG;;
        f)
            mnv_filter_switch=$OPTARG;;
        s)
            mnv_sample_index=$OPTARG;;
        j)
            multi_cutoff=$OPTARG;;
        M)
            identifyMulti=$OPTARG;;
        D)
            delCache=$OPTARG;;
        c)
            chr_list=$OPTARG;;
        esac
    done

mnv_pnumber=`expr $mnv_distance + 1`

#Showing the package version
echo -e "
====================================================
- MNVANNO [version 1.0]
====================================================\n"

#Showing help if there is not any input
tmp=$1
if [[ x${tmp} = "x" || ${tmp} = "-h" || ${tmp} = "-help" || ${tmp} = "--help" ]]
then
echo -e "Quick Start:

Test example: bash [scriptPath]/MNVIdentify.sh -i [scriptPath]/example/demo.vcf
Parameters: bash [scriptPath]/MNVIdentify.sh -i [file url] -d [identification distance] -a [cutoff of AC] -m [cutoff of MAF] -f [switch] -c [chr list]

\nRequired parameters:
-i <str> [file url] The absolute path of the VCF file;
-d <int> [identification distance] Max identification distance of MNV; Range of the parameter is >=2; Default is 10;
-a <int> [AC/adjust AC] Min AC/adjust_AC of MNV; Default is 0 (means >=1);
-m <int> [adjust MAF] Min adjust MAF of MNV; Default is 0;
-f <int> [filter switch] If 'T', using adjust AC and MAF to filter data; if 'F', using AC to filter data; Default 'T';
-s <int> [sample index] If 'T', indexs of sample containing MNV are output; Default 'F';
-M <int> [identify multiple MNV] If 'T', identify multi-MNV; Default 'F';
-j <int> [joint number for identifyMulti] Max joint number for identifyMulti; Range of the parameter is >=2; Default is 6;
-D <int> [delete cache] Default 'T';
-c <str> [chr list] Chromosomes list of vcf file; Multiple alignments must be in a comma separated list, such as 1,2,3,X,Y; [for uncompressed VCF file, no need to provide]

\n====================================================\n
"
exit
fi
###################################################

###################################################
echo -e "[$(date +%R:%S)] Start"
echo -e "===================================================="
echo -e "[$(date +%R:%S)] Check parameters"
# File url
if [ $file_url == "0" ]
then
echo -e "Please input the absolute path of the VCF file"
exit
fi
# Distance
if [[ $mnv_distance -lt 2 ]]
then 
echo -e "Range of the parameter is >= 2"
exit
fi
# AC
if [[ $mnv_ac -lt 0 ]]
then 
echo -e "mnv_ac must be >= 0"
exit
fi
# MAF
if [[ `echo "$mnv_maf <= 1.0" |bc` -eq 0 ]] || [[ `echo "$mnv_maf >= 0" |bc` -eq 0 ]]
then
echo $mnv_maf is out
exit
fi
# Joint number for identifyMulti
if [[ $multi_cutoff -lt 2 ]]
then 
echo -e "Range of the parameter is >= 2"
exit
fi
# Chr list
if [[ ${file_url:0-2} = 'gz' ]] && [[ $chr_list == '0' ]]
then
echo -e "For compressed VCF file, user need to fill in this parameter."
exit
elif [[ ${file_url:0-2} != 'gz' ]]
then
echo -e "[$(date +%R:%S)] Extract all chrID from mnv.vcf"
chr_list=($(awk '{print $1}' ${file_url} |grep -v '#'|sort -u|paste -s -d ","|sed 's/chr//g')) 
fi
# Create files
file_name=${file_url/.gz/}
file_name=${file_name/.vcf/}
mkdir $file_name
# Filter data
echo -e "===================================================="
python "${MNVIdentifyfolder}"VCFCheck.py -i $file_url -d $mnv_distance -c $chr_list
echo -e "===================================================="
# Determine whether filter.vcf file has data
line_count=$(cat ${file_name}'/filter.vcf' | wc -l)
if [ $line_count -gt 1 ]; then
# 在这里运行需要执行的代码

# Identify MNV
python "${MNVIdentifyfolder}"MNVIdentify.py -i ${file_name}'/filter.vcf' -d $mnv_distance -a $mnv_ac -m $mnv_maf -f $mnv_filter_switch -s $mnv_sample_index -c $chr_list
# Identify multiple MNV
if [ $identifyMulti == "T" ]
then
python "${MNVIdentifyfolder}"MNVIdentifyMulti.py -i ${file_name}'/filter.vcf' -d $mnv_distance -a $mnv_ac -m $mnv_maf -f $mnv_filter_switch -s $mnv_sample_index -j $multi_cutoff -c $chr_list
fi
grep -v '*' ${file_name}'/mnv.txt' > ${file_name}'/mnv2.txt'
rm  ${file_name}'/mnv.txt'
mv  ${file_name}'/mnv2.txt' ${file_name}'/mnv.txt'
# Delete cache
if [ $delCache == "T" ]
then
rm ${file_name}'/filter.vcf'
rm ${file_name}'/multiple_alleles.list'
fi
echo -e "===================================================="
echo -e "[$(date +%R:%S)] End\n"

else
rm ${file_name}'/filter.vcf'
rm ${file_name}'/multiple_alleles.list'
echo "no variants at ${mnv_distance} bp distance"
echo -e "===================================================="
echo -e "[$(date +%R:%S)] End\n"
fi
