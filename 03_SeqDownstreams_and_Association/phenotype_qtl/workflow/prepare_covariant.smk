# -*- coding: UTF-8 -*-
#
# FileName     : prepare_covariant
# Author       : EastsunW eastsunw@foxmail.com
# Create at    : 2024-11-02 12:42
# Last Modified: 2024-11-02 12:42
# Modified By  : EastsunW
# -------------
# Description  :
# -------------


# 所有的QTL都用SNP的PCA作为人群结构的协变量
rule covariant_population:
    input:
        multiext(
            "results/SNP-biochemistry/SNP.biochem_common",
            ".bed",
            ".bim",
            ".fam",
        ),
    output:
        val="results/common_data/covariant/population.covariant.eigenval",
        vec="results/common_data/covariant/population.covariant.txt",
    log:
        "logs/covariant/population_PCA.log",
    threads: 10
    shell:
        """
        scripts/tools/plink-1.9/plink \
            --silent \
            --bfile results/SNP-biochemistry/SNP.biochem_common \
            --pca 5 \
            --out results/common_data/covariant/population.covariant
        mv results/common_data/covariant/population.covariant.eigenvec {output.vec}
        """


rule covariant_clinic:
    input:
        "resources/sample_info/sample_info.txt",
    output:
        "results/common_data/covariant/clinic.covariant.txt",
    params:
        id_colname="ID",
        age_colname="age",
        gender_colname="gender",
        gender_map="男=1,女=2",
    run:
        with open(input[0], "rt") as INPUT:
            header = INPUT.readline().strip().split("\t")
            id_idx = header.index(params.id_colname)
            age_idx = header.index(params.age_colname)
            gender_idx = header.index(params.gender_colname)
            gender_map = dict([x.split("=") for x in params.gender_map.split(",")])
            with open(output[0], "wt") as OUTPUT:
                OUTPUT.write("ID\tage\tgender\n")
                for line in INPUT:
                    line_splited = line.strip().split("\t")
                    gender_code = (
                        gender_map[line_splited[gender_idx]]
                        if line_splited[gender_idx] in gender_map
                        else "0"
                    )
                    OUTPUT.write(
                        "\t".join(
                                [line_splited[id_idx], line_splited[age_idx], gender_code]
                            )
                        + "\n"
                    )


rule covariant_peer:
    input:
        "results/common_data/phenotype/{phenotype}.quantity.filtered.txt",
    output:
        "results/common_data/phenotype/{phenotype}.covariant.peer.txt",
    conda:
        "../envs/peer.yaml"
    shell:
        "Rscript scripts/tools/compute_peer.R -i {input} -o {output} -n 15"


rule biochem_covariant_common:
    input:
        clinic=rules.covariant_clinic.output,
        population=rules.covariant_population.output.vec,
    output:
        "results/common_data/covariant/biochem_common.covariant.txt",
    run:
        import pandas as pd
        clinic = pd.read_csv(str(input["clinic"]), sep="\t").rename(
            columns={"ID": "IID"}, inplace=False
        )
        population = pd.read_csv(
            input["population"],
            sep=" ",
            header=None,
            names=["FID", "IID"] + [f"PC{i}" for i in range(1, 6)],
        )
        merged_cov = pd.merge(clinic, population, on="IID", how="right")
        new_order = ["FID"] + [col for col in merged_cov.columns if col != "FID"]
        merged_cov = merged_cov.reindex(columns=new_order)
        merged_cov.to_csv(output[0], sep="\t", index=False, header=True)


for biochem_type in ["male", "female"]:

    rule:
        name:
            f"biochem_covariant_{biochem_type}"
        input:
            common_cov="results/common_data/covariant/biochem_common.covariant.txt",
            gender_pheno=f"results/common_data/phenotype/biochem_{biochem_type}.norm.txt",
        output:
            f"results/common_data/covariant/biochem_{biochem_type}.covariant.txt",
        run:
            import pandas as pd

            common_cov = pd.read_csv(
                input["common_cov"],
                sep="\t",
                header=None,
                names=["FID", "IID"]
                + ["age", "gender"]
                + [f"PC{i}" for i in range(1, 6)],
            )
            common_cov.drop(columns=["gender"], inplace=True)
            gender_pheno = pd.read_csv(
                input["gender_pheno"],
                sep="\t",
                header=None,
                names=["FID", "IID"] + [f"phenotype_{i}" for i in range(1, common_cov.shape[1] - 1)],
            )
            # 找到gender_pheno中的FID
            gender_cov = common_cov[common_cov["IID"].isin(gender_pheno["IID"])]
            gender_cov.to_csv(output[0], sep="\t", index=False, header=True)


for variant in ["SNP", "InDel", "MNV", "SV"]:

    rule:
        name:
            f"copy_covariant_{variant}"
        input:
            clinic=rules.covariant_clinic.output,
            population=rules.covariant_population.output.vec,
            peer=rules.covariant_peer.output,
        output:
            clinic=f"results/{variant}-{{phenotype}}/filtered_data/clinic_covariant.txt",
            population=f"results/{variant}-{{phenotype}}/filtered_data/pop_covariant.txt",
            peer=f"results/{variant}-{{phenotype}}/filtered_data/peer_covariant.txt",
        shell:
            """
            cp {input.clinic} {output.clinic}
            cp {input.population} {output.population}
            cp {input.peer} {output.peer}
            """

    rule:
        name:
            f"copy_covariant_biochem_{variant}"
        input:
            "results/common_data/covariant/biochem_{biochem_type}.covariant.txt",
        output:
            f"results/{variant}-biochemistry/{variant}.biochem_{{biochem_type}}.covariant",
        shell:
            "cp {input} {output}"
