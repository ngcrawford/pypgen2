from ngs_parsers import VCF

vcf = VCF.VCF('/Users/testudines/DATA/vcf_5M_50k_SNPs/RNAseq_sample_lab_ids_v2.vcf.gz')


for c, i in enumerate(vcf.vcf_file_iterator()):
    if c > 5: break
    print i

    # ac = {}
    # for s in samples:

    #     if i[s] is None:
    #         ac[s] = None

    #     elif i[s]['GT'] == '0/0':
    #         ac[s] = 0

    #     elif i[s]['GT'] == '1/0' or i[s]['GT'] == '0/1':
    #         ac[s] = 1

    #     elif i[s]['GT'] == '1/1':
    #         ac[s] = 2

    # df.append(ac)