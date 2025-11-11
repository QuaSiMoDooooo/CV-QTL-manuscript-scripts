# -*- coding: UTF-8 -*-
#
# FileName     : prepare_covariant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-02 10:45
# Last Modified: 2024-11-02 10:45
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


rule variant_to_plink_common:
    input:
        sample_info="resources/sample_info/sample_info.txt",
        genotype="results/common_data/variant/{variant}.genotype.filtered.txt",
        position="results/common_data/variant/{variant}.position.txt",
    output:
        ped_file="results/{variant}-biochemistry/{variant}.biochem_common.ped",
        map_file="results/{variant}-biochemistry/{variant}.biochem_common.map",
    params:
        variant=lambda wildcards: wildcards.variant,
        id_colname="ID",
        sample_prefix="HN-",
        age_colname="age",
        gender_colname="gender",
        gender_map="男=1,女=2",
    script:
        "../scripts/data_prepare/variant2plink.py"


for biochem_type in ["male", "female"]:

    rule:
        name:
            f"variant_to_plink_{biochem_type}"
        input:
            common_ped="results/{variant}-biochemistry/{variant}.biochem_common.ped",
            common_map="results/{variant}-biochemistry/{variant}.biochem_common.map",
            gender_pheno=f"results/common_data/phenotype/biochem_{biochem_type}.norm.txt",
        output:
            ped_file=f"results/{{variant}}-biochemistry/{{variant}}.biochem_{biochem_type}.ped",
            map_file=f"results/{{variant}}-biochemistry/{{variant}}.biochem_{biochem_type}.map",
        run:
            import shutil

            gender_samples = []
            with open(input["gender_pheno"], "rt") as INPUT:
                for line in INPUT:
                    line_splited = line.strip().split("\t")
                    gender_samples.append(line_splited[1])
            with open(output["ped_file"], "wt") as OUTPUT_PED:
                with open(input["common_ped"], "rt") as INPUT_PED:
                    for line in INPUT_PED:
                        line_splited = line.strip().split("\t")
                        if line_splited[1] in gender_samples:
                            OUTPUT_PED.write(line)
                            # map文件原样输出
            shutil.copy(input["common_map"], output["map_file"])



rule make_variant_ref:
    input:
        "results/common_data/variant/{variant}.position.txt",
    output:
        "results/common_data/variant/{variant}.refallele",
    shell:
        """
        awk -F'\t' 'NR>1 {{gsub(",", "_", $5); if ($5 == ".") $5="REF"; print $1 "\t" $5}}' {input} > {output}
        """


rule plink_text_to_plink_binary:
    input:
        ped_file="results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.ped",
        map_file="results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.map",
        ref_allele=rules.make_variant_ref.output,
    output:
        multiext(
            "results/{variant}-biochemistry/{variant}.biochem_{biochem_type}",
            ".bed",
            ".bim",
            ".fam",
        ),
    log:
        "logs/biochemistry/{variant}.biochem_{biochem_type}.to_plink.log",
    shell:
        """
        (scripts/tools/plink-1.9/plink \
            --pedmap results/{wildcards.variant}-biochemistry/{wildcards.variant}.biochem_{wildcards.biochem_type} \
            --a2-allele {input.ref_allele} 2 \
            --make-bed \
            --out results/{wildcards.variant}-biochemistry/{wildcards.variant}.biochem_{wildcards.biochem_type}) > {log} 2>&1
        """


# rule plink_binary_postprocess:
#     input:
#         position="results/common_data/variant/{variant}.position.txt",
#         bim="results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.bim.tmp",
#     output:
#         "results/{variant}-biochemistry/{variant}.biochem_{biochem_type}.bim",
#     script:
#         "../scripts/data_prepare/plink_bim_postprocess.py"
