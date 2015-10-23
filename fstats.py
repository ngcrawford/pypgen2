import gzip
import pandas
from ngs_parsers import VCF

vcf = VCF.VCF('ngs_parsers/test/example.vcf.gz')
vcf.populations = {'pop1':['c511','c512','c513','c514','c515',
                           'c563','c614','c630','c639','c640'],
                   'pop2':['m523','m524','m525','m589','m675',
                           'm676','m682','m683','m687','m689']}


for c, i in enumerate(vcf.vcf_file_iterator()):
    print vcf.calc_mean_heterozygosity(i)
    #if c > 15: break



