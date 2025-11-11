with open('../01.data/glnexus.output.vcf', 'r') as infile, open('../01.data/glnexus.output.GT.vcf', 'w') as outfile:
    for line in infile:
        if line.startswith('#'):
            # 保留 VCF 文件头部信息
            outfile.write(line)
        else:
            parts = line.rstrip().split('\t')
            format_field = parts[8]  # FORMAT 字段
            genotype_fields = parts[9:]  # 基因型字段

            # 只保留 GT 字段
            format_list = format_field.split(':')
            gt_index = format_list.index('GT')

            new_genotype_fields = [genotype.split(':')[gt_index] for genotype in genotype_fields]

            # 更新基因型字段信息并写入新文件
            parts[8] = 'GT'
            parts[9:] = new_genotype_fields
            outfile.write('\t'.join(parts) + '\n')
