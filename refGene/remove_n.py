#remove \n
# with open('/Users/nestmilk/Documents/aijin/Homo_sapiens.GRCh38.dna.primary_assembly.fa','r') as f, open('/Users/nestmilk/Documents/aijin/GRCh38.fasta', 'w') as r:
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

with open('/Users/nestmilk/Documents/aijin/GRCh38.fasta','r') as f, open('/Users/nestmilk/Documents/aijin/GRCh38.fa', 'w') as r:
    i = 0
    for line in f:
        i += 1
        line = line.strip('\n')
        if i <= 50:
            r.write(line + '\n')
        else:
            break
