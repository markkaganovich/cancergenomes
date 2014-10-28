import sys
import operator
import rebuildanal_helper
from collections import Counter
import Queue
import logging
import simplejson as json
import os


###
# do it from json files
# chrom, start_position, rsid, reference_allele, alt (list), hugo_symbol,
# variant_classification
# amino acid residue (gene:pos)
###
#mutations = json.load(open('mutations.json'))
residues = json.load(open('residues.json')) #nonsilent mutations

### dictionary format
residues_dicts = []
for v in residues:
    residues_dicts.append({'chrom': v[0], 'start_position' : v[1], 'rsid' :
v[2], 'ref': v[3], 'alt' : v[4], 'hugo_symbol': v[5], 'variant_classification' :
v[6], 'residue': v[7]})

class Rows(object):
    def __init__(self,object):
        for x in object.keys():
            setattr(self, x, object[x])

tcga_rows = map(lambda x: Rows(x), residues_dicts)
counts_aa = json.load('counts_aa.json')




















'''
tcga_silent = filter(lambda x: x[6] == 'Silent', mutations)
tcga_nonsilent = filter(lambda x: x[6] != 'Silent', mutations)        


for r in tcga_nonsilent:
        if r[0] + ':' + str(r[1]) in snps:
                        r.append( snp_to_aa[r[0] + ':' + str(r[1])])
                        tcga_residues.append(r)

'''



    




