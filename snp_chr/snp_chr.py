with open("snp151Common.txt", "r") as r, open("snp151Chr.txt", "w") as w:
    for line in r:
        list=line.split("\t")
        list[1]=list[1].replace("chr","")
        line="\t".join(list)
        w.write(line)

