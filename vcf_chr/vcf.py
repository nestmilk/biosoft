with open("/Users/nestmilk/Downloads/small_exac_common_3_b37.vcf", 'r') as f, \
    open("/Users/nestmilk/Downloads/small_exac_common_3_b37.chr.vcf", 'w') as r:
    for line in f:
        if line.find('#') == 0:
            r.write(line)
        else:
            r.write("chr" + line )
