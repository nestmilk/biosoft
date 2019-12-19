

# 将换行去除
# with open('/Users/nestmilk/Documents/aijin/ucsc.hg19.fasta','r') as f, open('/Users/nestmilk/Documents/aijin/hg19.fasta', 'w') as r:
#     seq = ""
#     for line in f:
#         line = line.strip('\n')
#         if line.find('>') == 0 :
#             if seq:
#                 r.write(seq + '\n')
#             seq = ""
#             r.write(line + '\n')
#         else:
#             seq += line
#     r.write(seq + '\n')


# 将chrY行后的所有其他contig去除
# with open('/Users/nestmilk/Documents/aijin/hg19.fasta','r') as f, open('/Users/nestmilk/Documents/aijin/hg19.fa', 'w') as r:
#     i = 0
#     for line in f:
#         i += 1
#         line = line.strip('\n')
#         if i <= 50:
#             r.write(line + '\n')
#         else:
#             break

import linecache
import math
import re
from tool.util import translate, reverse_complement

#hg19_refGene中start位置是从0开始的,end如果按0开始已经+1看起来像从1开始
#生成的vcf文件是从1开始的，end刚好与refGene中的end一致
chr_dic = ["chrM","chr1","chr2","chr3","chr4","chr5","chr6",
           "chr7","chr8","chr9","chr10","chr11","chr12","chr13",
           "chr14","chr15","chr16","chr17","chr18","chr19","chr20",
           "chr21","chr22","chrX","chrY"]


def read_chr(chr):
    line = linecache.getline('/Users/nestmilk/Documents/aijin/hg19.fa', (chr_dic.index(chr)+1)*2)
    return line


def get_transcript_list(query):
    list = []
    with open('hg19_ensGene.txt','r') as e, open('hg19_refGene.txt','r') as r:
        if query.find('ENST') != -1:
            for line in e:
                line = line.strip('\n')
                if line.find(query) > 0:
                    list = line.split('\t')
                    break
        else:
            for line in r:
                line = line.strip('\n')
                if line.find(query) != -1:
                    list = line.split('\t')
                    break
    if len(list) == 0 :
        raise RNAError("query nothing")
    if list[1] != query :
        raise RNAError("query is too short")
    return list


class RNAError(Exception):
    def __init__(self, ErrorInfo):
        self.errorInfo = ErrorInfo
    def __str__(self):
        return self.errorInfo


def get_cds_exon_start_end(list):
    begin = -1
    end = -1
    exon_num = len(list)
    for i in range(0, exon_num):
        if list[i] != -1:
            begin = i
            break

    for i in range(exon_num - 1, -1, -1):
        if list[i] != -1:
            end = i
            break
    return [begin, end]


def get_stop_index(seq):
    aa_seq = translate(seq)
    stop_index = aa_seq.find('*')
    if stop_index == -1:
        return -1
    else:
        return (stop_index + 1)*3

#获取原始多肽片段
def get_original(seq, pos1, pos2, left_aa, right_aa):
    if math.floor((pos1)/3) > left_aa:
        start = (math.floor(pos1/3) -left_aa)*3
    else:
        start = 0

    if math.ceil((pos2-1)/3) + right_aa > len(seq)/3:
        end = len(seq)
    else:
        end = (math.ceil((pos2-1)/3) + right_aa)*3
        pass
    return seq[start:end]

#获取突变多肽片段
def get_insert(seq, pos1, pos2, fragment, left_aa, right_aa):
    if math.floor(pos1/3) > left_aa:
        start = (math.floor(pos1/3) -left_aa)*3
    else:
        start = 0

    left_intact = seq[start:math.floor(pos1 / 3) * 3]
    left_hang = seq[math.floor(pos1 / 3) * 3:pos1]
    right_hang = seq[pos2 - 1:math.ceil((pos2-1)/3)*3]
    zone = left_hang + fragment + right_hang
    if (pos2-pos1-1)%3 != len(fragment)%3:
        zone_to_end = zone + seq[math.ceil((pos2-1)/3)*3:len(seq)]
        zone_to_end_intact = zone_to_end[0:math.floor(len(zone_to_end)/3)*3]
        stop_index = get_stop_index(zone_to_end_intact)
        if stop_index == -1:
            if right_aa < len(zone_to_end_intact) / 3:
                zone_right = zone_to_end_intact[0:right_aa * 3]
            else:
                zone_right = zone_to_end_intact
        else:
            zone_right = zone_to_end_intact[0:stop_index]
        contig = left_intact + zone_right
    else:
        stop_index = get_stop_index(zone)
        if stop_index == -1:
            if math.ceil((pos2-1)/3) + right_aa > len(seq)/3:
                end = len(seq)
            else:
                end = (math.ceil((pos2-1)/3)  + right_aa)*3
            right = seq[math.ceil((pos2-1)/3)*3:end]
            contig = left_intact + zone + right
        else:
            right = zone[0:stop_index]
            contig = left_intact + right

    return contig


class mRNA:
    def __init__(self, query):
        try:
            list = get_transcript_list(query)
        except RNAError as e:
            print("create mRNA: " + str(e) + ", query: " + query)
            exit(1)
        self.name = list[1]
        self.chr = list[2]
        if list[3] == '-':
            self.strand = False
        else:
            self.strand = True
        self.transcript_start = int(list[4])
        self.transcript_end = int(list[5])
        self.cds_start = int(list[6])
        self.cds_end = int(list[7])
        self.exon_num = int(list[8])
        self.exon_start_list = [int(i) for i in list[9].split(',')[:-1]]
        self.exon_end_list = [int(i) for i in list[10].split(',')[:-1]]
        self.exon_frame_list = [int(i) for i in list[15].split(',')[:-1]]
        #cds的染色体起始序号序列
        self.cds_start_list = []
        #cds的染色体终止序号序列
        self.cds_end_list = []
        #cds所在外显子的序列编号
        self.cds_exon_start_end = get_cds_exon_start_end(self.exon_frame_list)
        #获得cds_start_list,cds_end_list
        self.__get_cds_start_end_list()
        #utr5的序列
        self.utr5_seq = ''
        #cds的序列
        self.cds_seq = ''
        #utr3的序列
        self.utr3_seq = ''
        self.__set_utr5_cds_utr3_seq()


    def __get_cds_start_end_list(self):

        begin = self.cds_exon_start_end[0]
        end = self.cds_exon_start_end[1]
        if begin == -1:
            return ""

        if begin == end:
            self.cds_start_list.append(self.cds_start)
            self.cds_end_list.append(self.cds_end)
        else:
            for i in range(begin, end + 1):
                if i == begin:
                    self.cds_start_list.append(self.cds_start)
                    self.cds_end_list.append(self.exon_end_list[i])
                elif i == end:
                    self.cds_start_list.append(self.exon_start_list[i])
                    self.cds_end_list.append(self.cds_end)
                else:
                    self.cds_start_list.append(self.exon_start_list[i])
                    self.cds_end_list.append(self.exon_end_list[i])
    
    def __set_utr5_cds_utr3_seq(self):
        chr_seq = read_chr(self.chr)
        cds_seq = ''
        start_seq = ''
        end_seq = ''
        for i in range(len(self.cds_start_list)):
            cds_seq += chr_seq[self.cds_start_list[i]:self.cds_end_list[i]]

        for i in range(self.cds_exon_start_end[0] + 1):
            if i == self.cds_exon_start_end[0]:
                start_seq += chr_seq[self.exon_start_list[i]:self.cds_start]
            else:
                start_seq += chr_seq[self.exon_start_list[i]:self.exon_end_list[i]]

        for i in range(self.cds_exon_start_end[1], self.exon_num):
            if i == self.cds_exon_start_end[1]:
                end_seq += chr_seq[self.cds_end:self.exon_end_list[i]]
            else:
                end_seq += chr_seq[self.exon_start_list[i]:self.exon_end_list[i]]

        if bool(1-self.strand):
            cds_seq = reverse_complement(cds_seq.upper())
            utr3_seq = reverse_complement(start_seq.upper())
            utr5_seq = reverse_complement(end_seq.upper())
        else:
            utr5_seq = start_seq.upper()
            utr3_seq = end_seq.upper()

        self.cds_seq = cds_seq
        self.utr5_seq = utr5_seq
        self.utr3_seq = utr3_seq

    def get_cds_len(self):
        n = 0
        for i in range(len(self.cds_start_list)):
            n += self.cds_end_list[i]-self.cds_start_list[i]
        return n

    def get_pos_in_cds(self, start, end):
        pos_start_cds_index = -1
        pos_end_cds_index = -1
        for i in range(len(self.cds_start_list)):
            if self.cds_start_list[i] <= start < self.cds_end_list[i]:
                pos_start_cds_index = i
            if self.cds_start_list[i] <= end < self.cds_end_list[i]:
                pos_end_cds_index = i
        if pos_end_cds_index == -1 or pos_start_cds_index == -1:
            raise RNAError("get_pos_in_cds: outside cds")
        c_start = 0
        c_end = 0
        for i in range(pos_start_cds_index+1):
            if i==pos_start_cds_index:
                c_start += start - self.cds_start_list[i]
                break
            c_start += self.cds_end_list[i] - self.cds_start_list[i]
        for i in range(pos_end_cds_index+1):
            if i==pos_end_cds_index:
                c_end += end - self.cds_start_list[i]
                break
            c_end += self.cds_end_list[i] - self.cds_start_list[i]
        if self.strand:
            return [c_start, c_end]
        else:
            cds_len = self.get_cds_len()
            return [cds_len + 1 - c_end, cds_len + 1 - c_start]




#试用于annovar产生的数据，插入时start和end相同
def get_mutation(query, start, end, ori, sub, left_aa, right_aa):
    rna = mRNA(query)
    if bool(1-rna.strand):
        if ori != '-':
            ori = reverse_complement(ori)
        if sub != '-' :
            sub = reverse_complement(sub)
    try:
        start_end = rna.get_pos_in_cds(start, end)
    except RNAError as e:
        return [str(e),"","","",""]
    cds_start = start_end[0]
    cds_end = start_end[1]

    chgvs = ""
    ori_aa = ""
    mut_aa = ""
    if cds_start == cds_end:
        if ori == '-':
            if rna.strand:
                #正链的插入，即ins
                chgvs = "c." + str(cds_start) + "_" + str(cds_end+1) + "ins" + sub
                ori_seq = get_original(rna.cds_seq, cds_start, cds_end+1, left_aa, right_aa)
                mut_seq = get_insert(rna.cds_seq, cds_start, cds_end+1, sub, left_aa, right_aa)
            else:
                #负链的插入，即ins
                chgvs = "c." + str(cds_start-1) + "_" + str(cds_end) + "ins" + sub
                ori_seq = get_original(rna.cds_seq, cds_start-1, cds_end, left_aa, right_aa)
                mut_seq = get_insert(rna.cds_seq, cds_start-1, cds_end, sub, left_aa, right_aa)
        else:
            ori_finded = rna.cds_seq[cds_start-1]
            if ori_finded == ori:
                if sub == '-':
                    chgvs = "c." + str(cds_start) + "_" + str(cds_start) + "del" + ori_finded
                    sub = ''
                else:
                    if len(sub) > 1:
                        #单碱基的不等长替换，即del+ins
                        chgvs = "c." + str(cds_start) + "del" + ori_finded + "ins" + sub
                    else:
                        #单碱基的替换,即>
                        chgvs = "c." + ori_finded + str(cds_start) + ">" + sub
                ori_seq = get_original(rna.cds_seq, cds_start-1, cds_start+1, left_aa, right_aa)
                mut_seq = get_insert(rna.cds_seq, cds_start-1, cds_start+1, sub, left_aa, right_aa)
            else:
                raise RNAError("ori_finded is different from ori!")
    else:
        ori_finded = rna.cds_seq[cds_start-1:cds_end]
        if ori_finded == ori:
            if sub == '-':
                #单碱基的del
                chgvs = "c." + str(cds_start) + "_" + str(cds_end) + "del" + ori_finded
                ori_seq = get_original(rna.cds_seq, cds_start-1, cds_end+1, left_aa, right_aa)
                mut_seq = get_insert(rna.cds_seq, cds_start-1, cds_end+1, "", left_aa, right_aa)
            else:
                if len(ori_finded) != len(sub):
                    #多碱基的不等长替换，即del+ins
                    chgvs = "c." + str(cds_start) + "_" + str(cds_end) + "del" + ori_finded + "ins" + sub
                else:
                    #多碱基的等长替换，即>
                    chgvs = "c." + str(cds_start) + "_" + str(cds_end) + ori_finded +  ">" + sub
                ori_seq = get_original(rna.cds_seq, cds_start-1, cds_end+1, left_aa, right_aa)
                mut_seq = get_insert(rna.cds_seq, cds_start-1, cds_end+1, sub, left_aa, right_aa)
        else:
            raise RNAError("ori_finded is different from ori!")
    return [chgvs, ori_seq, mut_seq, translate(ori_seq), translate(mut_seq)]

def get_ensemble_gene():
    enst_name = {}
    with open("ensemblToGeneName.txt") as f:
        for line in f:
            line = line.strip("\n")
            list = line.split('\t')
            enst_name[list[0]]=list[1]
    return enst_name



enst_name = get_ensemble_gene()
with open("ensGene.txt", "r") as f, open("result.txt", "w") as r:
    i = 1
    for line in f:
        line = line.strip('\n')
        list = line.split(',')
        list[0] = list[0].split('\t')[-1]
        pre_list = list[:-1]
        last_list = list[-1].split('\t')
        chr = last_list[1]
        start = int(last_list[2])
        end = int(last_list[3])
        ori = last_list[4]
        sub = last_list[5]
        for item in pre_list:
            sub_list = item.split(':')
            enst = sub_list[1]
            if len(sub_list) == 3:
                enst_chgvs = ""
                enst_phgvs = ""
            elif len(sub_list) == 4:
                enst_chgvs = sub_list[3]
                enst_phgvs = ""
            else:
                enst_chgvs = sub_list[3]
                enst_phgvs = sub_list[4]
            result_list = get_mutation(enst, start, end, ori, sub, left_aa = 7, right_aa = 7)
            chgvs = result_list[0]
            ori_seq = result_list[1]
            mut_seq = result_list[2]
            ori_aa = result_list[3]
            mut_aa = result_list[4]
            r.write(str(i) + "\t" + enst_name[enst] + "\t" + enst + "\t" +
                    chr + "\t" + str(start) +"\t" + str(end) + "\t" +
                    enst_chgvs + "\t" + chgvs + "\t" + ori_seq + "\t" + mut_seq +
                    "\t" + enst_phgvs + "\t" + ori_aa + "\t" + mut_aa + "\n")
        i += 1



# get_mutation("NM_000314", 89711968, 89711968, "-", "G")
#
#
# rna = mRNA("NM_000314")
# b= get_original(rna.cds_seq, 586,587, 10, 10)
# print(b)
# print(translate(b))
#
# a = get_insert(rna.cds_seq, 586,587, "G",10,10)
# print(a)
# print(translate(a))




