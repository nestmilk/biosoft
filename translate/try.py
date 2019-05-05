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


