import itertools,argparse,time

########Help############################################################
parser = argparse.ArgumentParser(description='this is a description')
parser.add_argument('-i',"--infile",required=True,dest="infile",help="The absolute path of the VCF file;")
parser.add_argument('-d',"--distance",dest="distance",help="Max identification distance of MNV; Default is 10;",default=10)
parser.add_argument('-a',"--ac",dest="ac",help="Min AC/adjust_AC of MNV; Default is 0 (means >=1);",default=0)
parser.add_argument('-m',"--maf",dest="maf",help="Min adjust MAF of MNV; Default is 0;",default=0)
parser.add_argument('-f',"--filter",dest="filterSwitch",help="If 'T', using adjust AC and MAF to filter data; if 'F', using AC to filter data; Default 'T';",default='T')
parser.add_argument('-s',"--sample",dest="sampleSwitch",help="If 'T', indexs of sample containing MNV are output; Default 'F';",default='F')
parser.add_argument('-j',"--joint",dest="joint",help="Max joint number for identifyMulti; Default is 6;",default=6)
parser.add_argument('-c',"--chr",required=True,dest="chr",help="Chromosomes list of vcf file;")
args=parser.parse_args()
####################################################################

print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Run MNVIdentityMulti.py')
########Load parameter############################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load parameter')
file_url = args.infile
mnv_distance = int(args.distance) 
mnv_ac = int(args.ac) 
mnv_maf = float(args.maf)
sample_switch = args.sampleSwitch
filter_switch = args.filterSwitch
multi_cutoff = int(args.joint) 
chr_list = args.chr
# joint number of MNV = max identification distance of MNV + 1
mnv_pnumber = mnv_distance+1
file_prefix='/'.join(file_url.split('/')[:-1])+'/'
output=file_prefix + 'mnv_multi.txt'
chr_list=chr_list.split(',')
####################################################################


########Functions############################################################
# Load MNV file
def readMNVData(chr_list,file_prefix):
    chr_dict={}
    for i in chr_list:
        chr_dict[i]=[]
    f=open(file_prefix + 'mnv.txt')
    for i in f:
        a=i.strip('\n').split('\t')[0:12]
        if a[11] == str(1):
            chr_dict[a[0]].append(a)
    f.close()
    return chr_dict
# Load VCF file
def readData(chr_list,file_url):
    # 创建以染色体为单元的字典
    chr_dict={}
    for i in chr_list:
        chr_dict[i]=[{},{},{},{}]
    f=open(file_url)
    for i in f:
        a=i.strip('\n').split('\t')
        if i[0]!='#':
            chr_dict[a[0]][0][int(a[1])]=a[5] # my_data_dict={}
            chr_dict[a[0]][1][int(a[1])]=a[3] # snv_ref={}
            chr_dict[a[0]][2][int(a[1])]=a[4] # snv_alt={}
            chr_dict[a[0]][3][int(a[1])]=a[2] # snv_rsid={}    
    f.close()
    return chr_dict
# Load multi.mnv file
def readMultiData(chr_list,file_prefix,chr_dict):
    multi_dict={}
    for i in chr_list:
        multi_dict[i]=[{},{},{},{}] # genetype ref alt rsid
    f=open(file_prefix + 'multiple_alleles.list')
    for i in f:
        a=i.strip('\n').split('\t')
        if int(a[1]) not in multi_dict[a[0]][0]:
            if int(a[1]) in chr_dict[a[0]][1]:
                multi_dict[a[0]][0][int(a[1])]=[a[5]] # snv_dict={}
                multi_dict[a[0]][1][int(a[1])]=[a[3]] # snv_ref={}
                multi_dict[a[0]][2][int(a[1])]=[a[4]] # snv_alt={}
                multi_dict[a[0]][3][int(a[1])]=[a[2]] # snv_rsid={}     
        else:
            multi_dict[a[0]][0][int(a[1])].append(a[5]) # snv_dict={}
            multi_dict[a[0]][1][int(a[1])].append(a[3]) # snv_ref={}
            multi_dict[a[0]][2][int(a[1])].append(a[4]) # snv_alt={}
            multi_dict[a[0]][3][int(a[1])].append(a[2]) # snv_rsid={}          
    f.close()
    return multi_dict
# get number of samples
def getSampleNum(file_url):
    f=open(file_url)
    k=0
    for i in f:
        if k>0:
            break
        else:
            sampleNum=len(i.strip('\n').split('\t')[5])/2 
            k=k+1
    f.close()
    return sampleNum
### Data Summation
def getPoint(all_point):
    sum_list = [0] * len(all_point[0])
    for i in range(len(all_point)):
        for ii,v in enumerate(all_point[i]):
            v_add = int(v) + sum_list[ii]
            sum_list[ii] = v_add
    return sum_list
# Expand ranges according to p and mnv_distance
def expandPoint(p,mnv_distance,pos_dict):
    pp=list(set(pos_dict[p[0]]+pos_dict[p[-1]]))
    pp.sort()
    start=p[-1]-mnv_distance
    end=p[0]+mnv_distance
    newpp=[]
    for n in pp:
        if n >=start and n <=end:
            newpp.append(n)
    return newpp
# Store N points (N=mnv_distance) of its left and right range for each point
def eachPointDict(pos_list,mnv_distance):
    pos_dict={}
    for i in range(len(pos_list)):
        if i-mnv_distance<0:
            pos_dict[pos_list[i]]=pos_list[0:(i+mnv_distance+1)]
        elif i+mnv_distance>len(pos_list)-1:
            pos_dict[pos_list[i]]=pos_list[(i-mnv_distance):len(pos_list)]
        else:
            pos_dict[pos_list[i]]=pos_list[(i-mnv_distance):(i+mnv_distance+1)]
    return pos_dict
# output one MNV
def writeData(p,chr_dict,single_chr,ac,adjust_ac,sampleNum,adjust_ac_index):
    ref=[]
    alt=[]
    rsid=[]
    for i in p:
        ref.append(chr_dict[single_chr][1][i])
        alt.append(chr_dict[single_chr][2][i])
        rsid.append(chr_dict[single_chr][3][i]) 
    return [single_chr,','.join(map(str, p)),'.',','.join(ref),','.join(alt),','.join(rsid),max(p)-min(p),len(ref),ac,adjust_ac,"%e"%(adjust_ac/(2*sampleNum)),1,','.join(str(e) for e in adjust_ac_index)]
####################################################################



########Data loading, processing and output################################################
# Open files
result=open(output,'w')
result.write('#chr\tpos\tMNVID\tref\talt\trsID\tdistance\tMNVType\tAC\tadjust_AC\tadjust_AF\tmultiple_alleles\tsample_index\n')
# Load data
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Load data')
mnv=readMNVData(chr_list,file_prefix)
chr_dict=readData(chr_list,file_url)
chr_multi_dict=readMultiData(chr_list,file_prefix,chr_dict)
sampleNum=getSampleNum(file_url)
# Cycle for chr
for single_chr in chr_list:
    print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'Process the chr' + single_chr)
    # Obtain positions of points
    pos_list=list(chr_dict[single_chr][0].keys())
    pos_list.sort()
    # Store N points (N=mnv_distance) of its left and right range for each point
    pos_dict=eachPointDict(pos_list,mnv_distance)
    # Extract and market mutiple SNV
    multi_list=[]
    for i in chr_multi_dict[single_chr][0]:
        for ii in range(len(chr_multi_dict[single_chr][0][i])):
            multi_list.append(str(i)+'.'+str(ii))
    multi_list = [ float(x) for x in multi_list ]
    # Cycle for MNV
    for i in mnv[single_chr]:
        if int(i[7])<=multi_cutoff:
            # Extract one MNV
            onemnv=i[1].split(',')
            A = [ float(x) for x in onemnv ]
            for ii in range(len(A)):
                a=A[ii]+0.1
                b=A[ii]+0.2
                c=A[ii]+0.3
                if a in multi_list:
                    A.append(a)
                if b in multi_list:
                    A.append(b)        
                if c in multi_list:
                    A.append(c)
            # Obtian all combinations for this MNV
            B=[]
            for ii in itertools.combinations(A,len(onemnv)):
                if len(set(list(map(int,ii)))) == len(onemnv):
                    a=list(ii)
                    a.sort(key=float) 
                    B.append(a)
            # Cycle for combinations
            for ii in B:
                bak=[[],[],[],[],[]]
                for iii in ii:
                    if int(iii) in chr_multi_dict[single_chr][0]:
                        bak[0].append(int(iii))
                        bak[1].append(chr_dict[single_chr][0][int(iii)])
                        bak[2].append(chr_dict[single_chr][1][int(iii)])
                        bak[3].append(chr_dict[single_chr][2][int(iii)])
                        bak[4].append(chr_dict[single_chr][3][int(iii)]) 
                        chr_dict[single_chr][0][int(iii)]=chr_multi_dict[single_chr][0][int(iii)][int(str(iii).split('.')[1])]
                        chr_dict[single_chr][1][int(iii)]=chr_multi_dict[single_chr][1][int(iii)][int(str(iii).split('.')[1])]
                        chr_dict[single_chr][2][int(iii)]=chr_multi_dict[single_chr][2][int(iii)][int(str(iii).split('.')[1])] 
                        chr_dict[single_chr][3][int(iii)]=chr_multi_dict[single_chr][3][int(iii)][int(str(iii).split('.')[1])]      
                         
                p=[ int(x) for x in ii ]
                k=len(p)
                # Extract genotype information for each point
                all_point=[]
                for n in p:
                    all_point.append(chr_dict[single_chr][0][n])
                # Sum data
                all_point=getPoint(all_point)
                # expand ranges according to p and mnv_distance
                p2=expandPoint(p,mnv_distance,pos_dict)
                # Extract genotype information for each point
                all_point2=[]
                for n in p2:
                    all_point2.append(chr_dict[single_chr][0][n])  
                # Sum data
                all_point2=getPoint(all_point2)
                # Calculate AC of MNV
                ac=all_point.count(k)
                # Calculate real AC of MNV
                tmp=[idx for idx,m in enumerate(all_point) if m == k] # 获取相等的索引
                adjust_ac=0
                adjust_ac_index=[]
                for n in tmp:
                    if all_point[n]==all_point2[n]:
                        adjust_ac=adjust_ac+1
                        adjust_ac_index.append(n)
                # Judge results
                if filter_switch=='T':
                    if adjust_ac > mnv_ac and min(adjust_ac/(2*sampleNum),(1-adjust_ac/(2*sampleNum))) >= mnv_maf:
                        if sample_switch=='T':
                            a=writeData(p,chr_dict,single_chr,ac,adjust_ac,sampleNum,adjust_ac_index)
                        else:
                            a=writeData(p,chr_dict,single_chr,ac,adjust_ac,sampleNum,'')
                        result.write('\t'.join(str(e) for e in a)+'\n')
                else:
                    if ac > mnv_ac:
                        if sample_switch=='T':
                            a=writeData(p,chr_dict,single_chr,ac,adjust_ac,sampleNum,adjust_ac_index)
                        else:
                            a=writeData(p,chr_dict,single_chr,ac,adjust_ac,sampleNum,'')
                        result.write('\t'.join(str(e) for e in a)+'\n')
                # Recovery chr_dict
                for n in range(len(bak[0])):
                    chr_dict[single_chr][0][int(bak[0][n])]=bak[1][n]
                    chr_dict[single_chr][1][int(bak[0][n])]=bak[2][n]
                    chr_dict[single_chr][2][int(bak[0][n])]=bak[3][n]
                    chr_dict[single_chr][3][int(bak[0][n])]=bak[4][n]
                    
# Close files
result.close()
####################################################################
print('[' + time.strftime("%H:%M:%S", time.localtime()) + '] ' + 'End MNVIdentifyMulti.py')



