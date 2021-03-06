#anal.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import commands
import os
import csv 
import headers
import operator
import matplotlib.pyplot as plt
import numpy as np


db = create_engine('sqlite:///tcga_somatic.db', echo = False)

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select().execute()

class Rows(object):
    def __init__(self,object):
        for x in object.keys():
            setattr(self, x, object[x])

f = filter(lambda x: x.variant_classification in first_pass_mutations, m)
rows = map(lambda x: Rows(x), f)

def count_snps(results_list, outputfile = 'snpcount'):
    snpcount = {}
    keys = set([])
    for x in results_list:
        snp = x.chrom+':'+x.start_position 
        if snp in keys:
            snpcount[snp] = snpcount[snp] + 1
        else:
            snpcount[snp]= 1
            keys.add(snp)
    out = open(outputfile, 'w')
    json.dump(snpcount, out)
    return snpcount

#get snp count
snpcountfile = 'snpcount'
if snpcountfile in os.listdir('./'):
    snpcount = json.load(open(snpcountfile, 'r'))
else:
    snpcount = count_snps(rows)

def count_genes(rows, outputfile = 'gene_snp_count'):
    gene_snp_count = {}
    for r in rows:
        gene = r.hugo_symbol
        try:
            gene_snp_count[gene] = gene_snp_count[gene] + 1
        except KeyError:
            gene_snp_count[gene] = 1
    out = open(outputfile, 'w')
    json.dump(gene_snp_count, out)
    return gene_snp_count

#get gene_snp count
genecountfile = 'gene_snp_count'
if genecountfile in os.listdir('./'):
    gene_snp_count = json.load(open(genecountfile, 'r'))
else:
    gene_snp_count = count_genes(rows)

def explore_snp_counts(snpcount, outputfile = 'counts.png'):
    a = sorted(snpcount.iteritems(), key = operator.itemgetter(1), reverse=True)
    counts = map(lambda x: x[1], a)
    snps = map(lambda x: x[0], a)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ind = list(range(0, len(counts)))
    width = .2
    print "making bar"
    bar = ax.bar(ind, counts, width, color="r")
    print "saving ..."
    plt.savefig(outputfile)

def count_types(rows):
    types = list(set(map(lambda x: x.cancer_type, rows)))
    for t in types:
        n = filter(lambda x: x.cancer_type == t, rows)
        samples = map(lambda x: x.tumor_sample_barcode, n)
        print t + ':' + str(len(set(samples)))

#filter snps by their frequency. only use those that appear > .....

print "about to make matrix...."


matrixfile = 'matrix'

def get_matrix(results_list, outputfile = 'matrix'):
    sample_matrix = {}
    #mutations = list(set(map(lambda x: x.chrom + ':' + x.start_position, f)))
    matrix_samples = set([])
    for x in results_list:
        sample = x.tumor_sample_barcode
        print sample
        snp = x.chrom + ':' + x.start_position
        if sample in matrix_samples:
            if snp in sample_matrix[sample]:
                sample_matrix[sample][snp] = sample_matrix[sample][snp] + 1
            else:
                sample_matrix[sample][snp] = 1
        else:
            sample_matrix[sample] = {}
            sample_matrix[sample][snp] = 1
            matrix_samples.add(sample)
    out = open(outputfile, 'w')
    json.dump(sample_matrix, out)

if matrixfile not in os.listdir('./'):
     get_matrix(rows, outputfile = matrixfile)
sample_matrix = json.load(open(matrixfile, 'r'))

def get_snp_sample_matrix(results_list, outputfile = 'snp_sample'):
    snp_sample = {}
    matrix_samples = set([])
    for x in results_list:
        sample = x.tumor_sample_barcode
        print sample
        snp = x.chrom + ':' + x.start_position
        if snp not in matrix_samples:
            snp_sample[snp] = []
            matrix_samples.add(snp)
        if sample not in snp_sample[snp]:
            snp_sample[snp].append(sample)
    out = open(outputfile, 'w')
    json.dump(snp_sample, out)
        
snp_sample_file = 'snp_sample'        
if snp_sample_file not in os.listdir('./'):
    get_snp_sample_matrix(rows, outputfile = snp_sample_file)
snp_sample = json.load(open(snp_sample_file, 'r'))

snp_gene = dict(map(lambda x: (x.chrom +':' + x.start_position, x.hugo_symbol), rows))

def get_gene_sample(snp_sample, snp_gene, outputfile = 'gene_sample'):
    gene_sample = {}
    snps = snp_sample.keys()
    for s in snps:
        gene = snp_gene[s]
        try:
            gene_sample[gene].extend(snp_sample[s])
        except KeyError:
            gene_sample[gene] = snp_sample[s]
    for k in gene_sample.keys():
        gene_sample[k] = list(set(gene_sample[k]))
    out = open(outputfile, 'w')
    json.dump(gene_sample, out)
    return gene_sample

gene_sample_file = 'gene_sample'
if gene_sample_file not in os.listdir('./'):
    get_gene_sample(snp_sample, snp_gene, outputfile = gene_sample_file)
gene_sample = json.load(open(gene_sample_file, 'r'))

genes = gene_sample.keys()

#plot gene_samples
'''
gene_sample_count = {}
for g in genes:
    gene_sample_count[g] = gene_sample[g].__len__()
explore_snp_counts(gene_sample_count, 'gene_sample_counts.png')
'''

genecount = {}
for g in genes:
    genecount[g] = len(gene_sample[g])

def convert_to_set(dic):
    for k in dic.keys():
        dic[k] = set(dic[k])
    return dic

gene_sample_set = convert_to_set(gene_sample)

def run_co_occur(gene_sample_set, outputfile = 'co_occur_np'):
    genes = gene_sample_set.keys()
    co = np.identity(len(genes))
    #co_occur = {}
    for i, gi in enumerate(genes):
        #co_occur[i] = {}
        for j, gj in enumerate(genes):
            #co_occur[i][j] = gene_sample_set[i].intersection(gene_sample[j]).__len__()
            co[i,j] = gene_sample_set[gi].intersection(gene_sample[gj]).__len__()
    #json.dump(co_occur, open(outputfile, 'w'))
    np.save(open(outputfile, 'w'), co)
    return co

'''
if 'co_occur_np' in os.listdir('./'):
    co = np.load(open('co_occur_np'))
else:
    co = run_co_occur(gene_sample_set)
'''
'''
prob = {}
for g in genes:
    prob[g] = gene_sample_set[g].__len__()
json.dump(prob, open('Prob_genes', 'w'))
'''

prob = json.load(open('Prob_genes'))

# translate co_occur matrix into conditional probabilityes

cond_co_occur = np.load(open('cond_co_occur'))
'''
cond_co_occur = np.identity(len(genes))
for i, gi in enumerate(genes):
    for j, gj in enumerate(genes): 
        cond_co_occur[i,j] = co[i,j] / float(prob[gi])

np.save(open('cond_co_occur', 'w'), cond_co_occur)
        
        #if i not in cond_co_occur.keys():
        #    cond_co_occur[i] = {}
        #cond_co_occur[i][j] = float(co_occur[i][j]) /  float(prob[i])
    
'''

# find most mutated but least likely for memo
p53 = 17721
np.where(cond_co_occur[p53] == 0.)
b = map(lambda x: len(gene_sample[genes[x]]), a[0])


bp_to_aa = json.load(open('bp_to_aa'))
# convert rows to amino acids, then run the co_occurrence stuff
for r in rows:
    try:
        setattr(r, 'residue', bp_to_aa[r.chrom+':'+r.start_position])
    except KeyError:
        print r.chrom +':' + r.start_position
        rows.remove(r)


#given a group of samples and genes, calculate the co occurrence stuff      
