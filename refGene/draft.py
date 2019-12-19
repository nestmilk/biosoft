from get_cds_for_ensembl_GRCh38_by_vep import get_mutation, RNAError
from tool.util import translate, reverse_complement

# def get_cds_by_mRNA1(query):
#     name = ''
#     chr = ''
#     strand = True
#     transcript_start = ''
#     transcript_end = ''
#     cds_start = ''
#     cds_end = ''
#     exon_num = 1
#     exon_start_list = []
#     exon_end_list = []
#     exon_frame_list = []
#     list = []
#
#     with open('hg19_ensGene.txt','r') as e, open('hg19_refGene.txt','r') as r:
#         if query.find('ENST') > 0:
#             for line in e:
#                 line = line.strip('\n')
#                 if line.find(query) > 0:
#                     list = line.split('\t')
#                     break
#         else:
#             for line in r:
#                 line = line.strip('\n')
#                 if line.find(query) > 0:
#                     list = line.split('\t')
#                     break
#     if len(list) == 0 :
#         return "query nothing"
#     if list[1] != query :
#         return "query is too short"
#
#     name = list[1]
#     chr = list[2]
#     if list[3] == '-':
#         strand = False
#     transcript_start = int(list[4])
#     transcript_end = int(list[5])
#     cds_start = int(list[6])
#     cds_end = int(list[7])
#     exon_num = int(list[8])
#     exon_start_list = [int(i) for i in list[9].split(',')[:-1]]
#     exon_end_list = [int(i) for i in list[10].split(',')[:-1]]
#     exon_frame_list = [int(i) for i in list[15].split(',')[:-1]]
#
#
#     begin = -1
#     end = -1
#     for i in range(0, exon_num):
#         if exon_frame_list[i] != -1:
#             begin = i
#             break
#
#     for i in range(exon_num-1, -1, -1):
#         if exon_frame_list[i] != -1:
#             end = i
#             break
#
#     line = read_chr(chr)
#
#     raw_cds = ''
#     if begin == end:
#         raw_cds = line[cds_start:cds_end]
#     else:
#         for i in range (begin, end+1):
#             if i == begin:
#                 seq = line[cds_start:exon_end_list[i]]
#             elif i == end:
#                 seq = line[exon_start_list[i]: cds_end]
#             else:
#                 seq = line[exon_start_list[i]:exon_end_list[i]]
#             raw_cds += seq
#     cds = ''
#     if strand:
#         cds = raw_cds
#     else:
#         cds = reverse_complement(raw_cds)
#
#     return cds

# import math
# print(str(math.ceil(0/3)))

# from refGene.get_cds import read_chr
# line = read_chr("chr7")
# print(line[55259513:55259515])

# from refGene.get_cds import mRNA
# from tool.util import translate
# rna1 = mRNA("ENST00000361452")
# cds1 = rna1.cds_seq
# print(cds1)
# cds1_before = cds1[:153]
# print(cds1_before)
# cds1_after = cds1[175:]
# print(cds1_after)
# print(translate(cds1_before + cds1_after))

# rna2 = mRNA("ENST00000397789")
# cds2 = rna2.cds_seq



with open('CaoGZ_HCC_output.vcf','r') as f, open('result.tsv', 'w') as r:
    line_num = 1
    for line in f:
        if line.find("#") == -1:
            line = line.strip("\n")
            list = line.split('\t')
            chr = list[0]
            start = int(list[1])
            ori = list[3]
            sub = list[4]
            CSQ_list = list[7].split(';')[-1].split(',')
            line_num_add = False
            for item in CSQ_list:
                transcript_list = item.split("|")
                symbol = transcript_list[3]
                ensg = transcript_list[4]
                enst = transcript_list[6]
                reported_exon = transcript_list[8]
                enst_version = transcript_list[10].split(":")[0] if transcript_list[10] else ""
                reported_CHGVS = transcript_list[10].split(":")[1] if transcript_list[10] and len(transcript_list[10].split(":")) > 1 else ""
                reported_PHGVS = transcript_list[11].split(":")[1] if transcript_list[11] and len(transcript_list[11].split(":")) > 1 else ""
                if reported_PHGVS != "":
                    try:
                        list = get_mutation(enst, start, ori, sub, 15, 15)
                        strand = list[0]
                        chgvs = list[1]
                        ori_seq = list[2]
                        mut_seq = list[3]
                        ori_aa = list[4]
                        mut_aa = list[5]
                        # finded_ori = list[2] if strand == "+" else reverse_complement(list[2])
                    except RNAError as e:
                        strand = "na"
                        chgvs = str(e)
                        ori_seq = "na"
                        mut_seq = "na"
                        ori_aa = "na"
                        mut_aa = "na"
                        # finded_ori = "na"
                    # if strand != "na":
                    r.write(str(line_num) + "\t" + symbol + "\t" + enst + "\t" + chr + "\t" +
                            str(start) + "\t" + ori + "\t" + sub + "\t" +
                            strand + "\t" + reported_CHGVS + "\t" + chgvs + "\t" +
                            reported_PHGVS + "\t" + ori_seq + "\t" + mut_seq  + "\t" +
                            ori_aa + "\t" + mut_aa + "\n")
                    line_num_add = True
            if line_num_add:
                line_num += 1



