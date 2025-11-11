# Variant-Phenotype Association

Identifying genetic variants associated with phenotypes, using QTL methods for transcriptomic molecular phenotypes and GWAS methods for clinical indicators.

# Data Sources

## Genetic Variants

- **SNP_InDel**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SNP/clair3_cohort/GRCh38_241031.cohort.vcf.gz resources/variants/snp_indel.raw.vcf.gz`
- **MNV**: `gzip -c /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/MNV/clair3/GRCh38_241031.mnv.withID.txt > resources/variants/MNV.raw.vcf.gz`
- **SV**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/SV/cohort_merged/GRCh38_241031.SV.vcf.gz resources/variants/SV.raw.vcf.gz`

## Phenotypes

- **APA Quantification**: `ln -s /home/wangdy/Projects/Weibin/RNAseq/results/GRCh38_241031/APA/GRCh38_241031.QAPA.txt resources/phenotypes/APA.position.raw.txt`
- **APA Positions**: `ln -s /home/wangdy/Projects/Weibin/RNAseq/results/GRCh38_241031/APA/GRCh38_241031.QAPA.txt resources/phenotypes/APA.quantity.txt`
- **Gene Positions**: `ln -s /home/data_admin/Reference/Annotation/gencode_47.full.gtf resources/phenotypes/expression.position.raw.txt`
- **Gene Quantification**: `ln -s /home/wangdy/Projects/Weibin/RNAseq/results/GRCh38_241031/expression/expression_cohort/GRCh38_241031.rawtpm.txt resources/phenotypes/expression.quantity.txt`
- **AS Quantification**: `ln -s /home/wangdy/Projects/Weibin/RNAseq/results/GRCh38_241031/AS/GRCh38_241031.AS.txt resources/phenotypes/splicing.position.raw.txt`
- **AS Positions**: `ln -s /home/wangdy/Projects/Weibin/RNAseq/results/GRCh38_241031/AS/GRCh38_241031.AS.txt resources/phenotypes/splicing.quantity.txt`
- **Methylation Quantification**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/methylation/methylation_cohort/GRCh38_241031.merged.txt resources/phenotypes/methylation.quantity.txt`
- **Methylation Positions**: `ln -s /home/wangdy/Projects/Weibin/WGS/results/GRCh38_241031/methylation/methylation_cohort/GRCh38_241031.merged.txt resources/phenotypes/methylation.position.raw.txt`

## Clinical Indicators

- **Common Indicators**: `ln -s /home/data_admin/Clinic/Biochemistry/Biochemistry.common.txt resources/phenotypes/biochemistry.common.txt`
- **Androgen Indicators**: `ln -s /home/data_admin/Clinic/Biochemistry/Biochemistry.male.txt resources/phenotypes/biochemistry.male.txt`
- **Estrogen Indicators**: `ln -s /home/data_admin/Clinic/Biochemistry/Biochemistry.female.txt resources/phenotypes/biochemistry.female.txt`
        