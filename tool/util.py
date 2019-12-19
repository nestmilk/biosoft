import math

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

def translate(sequence, step=0):
    seq = sequence.upper()[step:]
    n = int(math.floor((len(seq)) / 3))
    peptide = ""
    for i in range(0, n):
        codon = seq[i * 3:i * 3 + 3]
        peptide += condon_table[codon]
    return peptide

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

