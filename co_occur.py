#co_occur.py
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

bp_to_aa = json.load(open('bp_to_aa'))
# convert rows to amino acids, then run the co_occurrence stuff
badrows = []
for r in rows:
    try:
        setattr(r, 'residue', bp_to_aa[r.chrom+':'+r.start_position])
    except KeyError:
        print r.chrom +':' + r.start_position
        badrows.append(r)
badrows = set(badrows)
print "about to make matrix...."


matrixfile = 'matrix'

def get_matrix(results_list, outputfile = 'matrix'):
    '''
    sample x snp matrix (for each sample, what are the residues)
    '''
    sample_matrix = {}
    #mutations = list(set(map(lambda x: x.chrom + ':' + x.start_position, f)))
    matrix_samples = set([])
    for x in results_list:
        sample = x.tumor_sample_barcode
        print sample
        try:
            residue = x.residue['pos']
        except AttributeError:
            continue
        if sample in matrix_samples:
            if residue in sample_matrix[sample]:
                sample_matrix[sample][residue] = sample_matrix[sample][residue] + 1
            else:
                sample_matrix[sample][residue] = 1
        else:
            sample_matrix[sample] = {}
            sample_matrix[sample][residue] = 1
            matrix_samples.add(sample)
    out = open(outputfile, 'w')
    json.dump(sample_matrix, out)

if matrixfile not in os.listdir('./'):
     get_matrix(rows, outputfile = matrixfile)
sample_matrix = json.load(open(matrixfile, 'r'))

'''
the next two functions are what you need to make the matrix that the co_occur function uses

'''
def get_snp_sample_matrix(results_list, outputfile = 'snp_sample'):
    '''
    make residue dictionary
    residue1: sample1, sample2, sample3
    '''
    snp_sample = {}
    matrix_samples = set([])
    for x in results_list:
        sample = x.tumor_sample_barcode
        print sample
        try:
            residue = str(x.residue['pos']) + ':'+ x.hugo_symbol
        except AttributeError:
            continue
        if residue not in matrix_samples:
            snp_sample[residue] = []
            matrix_samples.add(residue)
        if sample not in snp_sample[residue]:
            snp_sample[residue].append(sample)
    out = open(outputfile, 'w')
    json.dump(snp_sample, out)
    return snp_sample
        
snp_sample_file = 'residue_sample'        
if snp_sample_file not in os.listdir('./'):
    get_snp_sample_matrix(rows, outputfile = snp_sample_file)
snp_sample = json.load(open(snp_sample_file, 'r'))


def get_gene_sample(snp_sample, outputfile = 'gene_sample'):
    '''
    group previous dictionary by genes
    '''
    gene_sample = {}
    snps = snp_sample.keys()
    for s in snps:
        gene = s.split(':')[1]
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
    get_gene_sample(snp_sample, outputfile = gene_sample_file)
gene_sample = json.load(open(gene_sample_file, 'r'))


genes = gene_sample.keys()

def convert_to_set(dic):
    for k in dic.keys():
        dic[k] = set(dic[k])
    return dic

gene_sample_set = convert_to_set(gene_sample)
snp_sample_set = convert_to_set(snp_sample)



def run_co_occur(gene_sample_set, genes, samples, outputfile = 'co_occur_np'):
    #genes = gene_sample_set.keys()
    #samples = set(samples)
    co = np.identity(len(genes))
    if len(samples) < len(all_samples):
        print len(samples)
        flag = 1
    else:
        flag = 0
    #co_occur = {}
    for i, gi in enumerate(genes):
        print gi
        #co_occur[i] = {}
        for j, gj in enumerate(genes):
            if flag:
            #co_occur[i][j] = gene_sample_set[i].intersection(gene_sample[j]).__len__()
                set1 = gene_sample_set[gi].intersection(samples)
                set2 = gene_sample_set[gj].intersection(samples)
                co[i,j] = len(set1.intersection(set2))
            else:
                co[i,j] = len(gene_sample_set[gi].intersection(gene_sample_set[gj]))
            #co[i,j] = gene_sample_set[gi].intersection(gene_sample[gj]).__len__()
    #json.dump(co_occur, open(outputfile, 'w'))
    print len(samples)
    np.save(open(outputfile, 'w'), co)
    return co

global all_samples
all_samples = set(map(lambda x: x.tumor_sample_barcode, rows))


if 'co_occur_np' in os.listdir('./'):
    co = np.load(open('co_occur_np'))
else:
    co = run_co_occur(gene_sample_set, genes, all_samples)


prob = {}
for g in genes:
    prob[g] = gene_sample_set[g].__len__()
json.dump(prob, open('Prob_genes', 'w'))
'''
prob = json.load(open('Prob_genes'))

# translate co_occur matrix into conditional probabilityes
cond_co_occur = np.identity(len(genes))
for i, gi in enumerate(genes):
    for j, gj in enumerate(genes): 
        cond_co_occur[i,j] = co[i,j] / float(prob[gi])

np.save(open('cond_co_occur', 'w'), cond_co_occur)
        
        #if i not in cond_co_occur.keys():
        #    cond_co_occur[i] = {}
        #cond_co_occur[i][j] = float(co_occur[i][j]) /  float(prob[i])
    
'''


'''
make the same co_occur function for residues, need to input residue frequency to determine cutoffs
and only use residues that are relatively frequent
'''
if 'counts_aa' in os.listdir('./'):
    counts_aa = json.load(open('counts_aa'))
else:
    make_counts_aa(genes, counts)

f=[]
for g in genes:
    f.extend([k+':'+ g for k in counts_aa[g].keys() if counts_aa[g][k] >2])
f_set = set(f)

f_rows = filter(lambda x: x not in badrows and str(x.residue['pos']) +':' +x.hugo_symbol in f_set, rows)    

res_sample = get_snp_sample_matrix(f_rows, 'f_rows.snp_sample')
#rs = convert_to_set(res_sample)
res_genes = get_gene_sample(res_sample, 'f_res_sample')
rg = convert_to_set(res_genes)
co_f = run_co_occur(rg, rg.keys(), all_samples, 'filter_gene_co')

gbm = filter(lambda x: x.cancer_type == 'GBM', rows)  
res_sample = get_snp_sample_matrix(gbm, 'f_rows.snp_sample')
res_genes = get_gene_sample(res_sample, 'f_res_sample')
rg = convert_to_set(res_genes)
rs = convert_to_set(res_sample)
co_f = run_co_occur(rg, rg.keys(), all_samples, 'filter_gene_co')

s = []
map(lambda x: s.extend(res_sample[x]), res_sample.keys())
len(set(s))


def get_co(gene1, gene2, gene_sample_set = gene_sample_set, samples = set(s)):
    set1 = gene_sample_set[gene1].intersection(samples)
    set2 = gene_sample_set[gene2].intersection(samples)
    sam = set1.intersection(set2)
    n = len(sam)
    return [n, sam]

def get_co2(gene1, gene_sample_set= snp_sample_set, samples=samples, genes = f, co=co_r):
    a =np.where(co[genes.index(gene1)] >0)
    g = map(lambda x: genes[int(x)], a[0])
    return g


#m = map(lambda x: len(filter(lambda y: y.tumor_sample_barcode == x, rows)), all_samples)
p53samples = map(lambda x: x.tumor_sample_barcode, filter(lambda x: x.hugo_symbol =='TP53', rows))
p53_m = map(lambda x: len(filter(lambda y: y.tumor_sample_barcode == x, rows)), list(set(p53samples)))

a = get_co('TP53', 'MDM2', gene_sample_set, all_samples)
a_m = map(lambda x: len(filter(lambda y: y.tumor_sample_barcode == x, rows)), list(set(a[1])))
scipy.stats.ttest_1samp(a_m, np.mean(m))


