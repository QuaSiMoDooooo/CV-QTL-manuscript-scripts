import argparse,time,gzip

########Help############################################################
parser = argparse.ArgumentParser(description='This is a description')
parser.add_argument('-i',"--infile",required=True,dest="infile",help="The absolute path of the VCF file;")
parser.add_argument('-d',"--distance",dest="distance",help="Max identification distance of MNV; Default is 10;",default=10)
parser.add_argument('-c',"--chr",required=True,dest="chr",help="Chromosomes list of vcf file;")
args=parser.parse_args()
####################################################################

########Functions############################################################
### Read compressed files
def openfile(filename, mode="r"):
    if filename[-3:] == '.gz':
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)
### Data Summation
def getPoint(all_point):
    sum_list = [0] * len(all_point[0])
    for i in range(len(all_point)):
        for ii,v in enumerate(all_point[i]):
            v_add = int(v) + sum_list[ii]
            sum_list[ii] = v_add
    return ''.join(str(e) for e in sum_list)
####################################################################

print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Run VCFCheck.py')
########Load parameters############################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load parameters')
file_url = args.infile
distance = int(args.distance)
chr_list = args.chr
file_prefix=file_url.replace('.gz','').replace('.vcf','')+'/'
output=file_prefix + 'filter.vcf'
output2=file_prefix + 'multiple_alleles.list'
####################################################################

########Data loading, processing and output################################################
# Create dicts
chr_list=chr_list.split(',')
chr_dict={}
for i in chr_list:
    chr_dict[i]=[]
# Open files
result=open(output,'w')
mresult=open(output2,'w')
# Load data
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load data')
with openfile(file_url) as f:
    if file_url[-3:] == ".gz":
        for i in f:
            i = i.decode('utf-8')
            if i[0]!='#':
                chr_dict[i[0:100].split("\t")[0].replace("chr","")].append(i.replace("chr",""))
    else:
        for i in f:
            if i[0]!='#':
                chr_dict[i[0:100].split("\t")[0].replace("chr","")].append(i.replace("chr",""))
f.close()
# Cycle for chr
for single_chr in chr_list:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Process chr' + single_chr)
    my_data_dict={}
    multi_data_dict={}
    for i in range(len(chr_dict[single_chr])):
        m=chr_dict[single_chr][i]
        if len(m)>=1000:
            a=m[0:1000].split('\t')
        else:
            a=m.split('\t')
        # .|. .| |. => 0|0 0| |0
        m=m.replace('.|','0|') 
        m=m.replace('|.','|0') 
        b=[m.index('\t0|') if '\t0|' in m else -1,m.index('\t1|') if '\t1|' in m else -1] 
        if -1 in b:
            b.remove(-1)
        b=min(b)
        b=m[b+1:-1].replace('|','').replace('\t','')
        bl=list(set(list(b)))
        chr_dict[single_chr][i]=None
        # Filter data
        if a[3] not in ['A','T','C','G'] or a[4] not in ['A','T','C','G'] or len(bl) > 2 : 
            continue
        # Splice data
        m='\t'.join(a[0:5])+'\t'+b+'\n'
        # Maker mutiple SNV
        if int(a[1]) in my_data_dict:
            if int(a[1]) not in multi_data_dict:
                multi_data_dict[int(a[1])]=[my_data_dict[int(a[1])]]
                multi_data_dict[int(a[1])].append(m)
            else:
                multi_data_dict[int(a[1])].append(m)
            kk=my_data_dict[int(a[1])].strip().split('\t')
            newb=getPoint([b,kk[5]])
            my_data_dict[int(a[1])]='\t'.join([kk[0],kk[1],kk[2],kk[3],'*',newb])+'\n'
        else:
            my_data_dict[int(a[1])]=m
    # Filter data again
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Filter distance >' + str(distance))
    point=list(my_data_dict.keys())
    point.sort()
    t=[]
    for i in range(len(point)):
        if i == 0:
            if point[i+1]-point[i]<=distance:
                t.append(point[i])
        elif i == len(point)-1:
            if point[i]-point[i-1]<=distance:
                t.append(point[i])
        else:
            if (point[i+1]-point[i]<=distance) or (point[i]-point[i-1]<=distance):
                t.append(point[i])
    # Write data
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Write chr' + single_chr)
    for i in t:
        result.write(my_data_dict[i])
    for i in multi_data_dict:
        for ii in multi_data_dict[i]:
            mresult.write(ii)
# Close files
result.close()
mresult.close()
####################################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'End VCFCheck.py')