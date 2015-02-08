import sys
import operator
import rebuildanal_helper
from collections import Counter
import Queue
import logging
import simplejson as json
import os
import bisect
import solvebio

solvebio.login()

###
# do it from json files
# chrom, start_position, rsid, reference_allele, alt (list), hugo_symbol,
# variant_classification
# amino acid residue (gene:pos)
###
#mutations = json.load(open('mutations.json'))
residues = json.load(open('residues.json')) 

### dictionary format
file = open('tcga_residues_dictionary')
tcga_residues = json.load(file)
file.close()

class Rows(object):
    def __init__(self,object):
        for x in object.keys():
            setattr(self, x, object[x])

tcga_rows = map(lambda x: Rows(x), tcga_residues)
counts_aa = json.load(open('counts_aa.json'))


residues_7 = json.load(open('tcga_residues_7.json'))
residues_only = map(lambda x: x['residue'], residues_7)





















'''
tcga_silent = filter(lambda x: x[6] == 'Silent', mutations)
tcga_nonsilent = filter(lambda x: x[6] != 'Silent', mutations)        


for r in tcga_nonsilent:
        if r[0] + ':' + str(r[1]) in snps:
                        r.append( snp_to_aa[r[0] + ':' + str(r[1])])
                        tcga_residues.append(r)

'''



    




