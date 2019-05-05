import math, sys, getopt, time

condon_table = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R", "AAT": "N", "AAC": "N",
    "GAT": "D", "GAC": "D", "TGT": "C", "TGC": "C", "CAA": "Q", "CAG": "Q",
    "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "AAA": "K", "AAG": "K", "TTT": "F", "TTC": "F", "CCT": "P", "CCC": "P",
    "CCA": "P", "CCG": "P", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "AGT": "S", "AGC": "S", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGG": "W", "TAT": "Y", "TAC": "Y", "GTT": "V", "GTC": "V", "GTA": "V",
    "GTG": "V", "TAA": "*", "TGA": "*", "TAG": "*"
}

def get_time():
    return time.strftime('%H:%M:%S------', time.localtime(time.time()))


def reverse_complement (sequence):
    reverse_seq = sequence[::-1]
    reverse_complement_seq = ""
    for i in reverse_seq:
        if i == "A":
            reverse_complement_seq += "T"
        elif i == "T":
            reverse_complement_seq += "A"
        elif i == "C":
            reverse_complement_seq += "G"
        elif i == "G":
            reverse_complement_seq += "C"
    return reverse_complement_seq


def read_fasta (file, rev_com):
    sequences_list = []
    with open(file, 'r') as f:
        i = 0
        for line in f:
            i += 1
            if i % 2 == 0:
                line = line.strip('\n')
                sequences_list.append(line)
                if rev_com:
                    rev_com_seq = reverse_complement(line)
                    sequences_list.append(rev_com_seq)
    print (get_time() + "DNA sequences: "+ str(int(i/2)) + ", ("+str(len(sequences_list))+ ", if double as reverse_complemented).")
    return sequences_list


def raw_peptide (sequence, start_num):
    n = int(math.floor((len(sequence) - start_num) / 3))
    new_sequence = sequence[start_num:]
    peptide = ""
    for i in range(0, n):
        codon = new_sequence[i*3:i*3+3]
        peptide += condon_table[codon]
    return peptide


def all_pure_peptides (file, peptide_len, rec_com):
    if rec_com:
        print (get_time() + "6 frames translation is starting.")
    else:
        print (get_time() + "3 frames translation is starting.")
    sequences_list = read_fasta(file, rec_com)
    pure_peptides = set()
    num = 0
    ignored_num = 0
    for sequence in sequences_list:
        peptides = []
        peptides.append(raw_peptide(sequence, 0))
        peptides.append(raw_peptide(sequence, 1))
        peptides.append(raw_peptide(sequence, 2))
        for peptide in peptides:
            inner_peptide_list = peptide.split('*')
            for inner_peptide in inner_peptide_list:
                if len(inner_peptide) >= peptide_len:
                    pure_peptides.add(inner_peptide)
                    num += 1
                else:
                    if len(inner_peptide)>0:
                        ignored_num += 1
    print (get_time() + "saved peptides: "+ str(num) + ", deleted peptides: " + str(ignored_num) + ".")
    return list(pure_peptides)


def main(argv):
    inputfile = ""
    outputfile = "peptides.txt"
    len = 8
    rev_com_required = False
    try:
        opts, args = getopt.getopt(argv, "hi:o:l:r:", ["ifile","ofile"])
    except getopt.GetoptError:
        print (get_time() + 'python translate.py -i <inputfile> -o <outputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print (get_time() + '''
             python main.py -i <inputfile> -o <outputfile>
            -h: help;
            -i or --ifile: inputfile of assembled sequences as fasta format;
            -o or --ofile: outputfile of finded peptides, default as peptides.txt;
            -l: integer, the minimal length of required peptides, default as 8;
            -r: boolean, if reverse_complemented sequeces are required, default as False.
            ''')
            sys.exit(1)
        elif opt in ('-i', "--ifile"):
            inputfile = arg
        elif opt in ('-o', "--ofile"):
            outputfile = arg
        elif opt in ('-l'):
            len = int(arg)
        elif opt in ('-r'):
            if arg.upper() == "TRUE":
                rev_com_required = True
    pure_peptides = all_pure_peptides(inputfile, len, rev_com_required)
    with open(outputfile,'w') as f:
        for peptide in pure_peptides:
            f.write(peptide + "\n")
    print (get_time() + "translation is completed!!!")


if __name__ == "__main__":
    main(sys.argv[1:])
