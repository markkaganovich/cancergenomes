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
        
snp_sample_file = 'residue_sample'        
if snp_sample_file not in os.listdir('./'):
    get_snp_sample_matrix(rows, outputfile = snp_sample_file)
snp_sample = json.load(open(snp_sample_file, 'r'))

#snp_gene = dict(map(lambda x: (x.residue['pos'], x.hugo_symbol), rows))

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



def run_co_occur(gene_sample_set, genes, samples, outputfile = 'co_occur_np'):
    #genes = gene_sample_set.keys()
    samples = set(samples)
    co = np.identity(len(genes))
    #co_occur = {}
    for i, gi in enumerate(genes):
        #co_occur[i] = {}
        for j, gj in enumerate(genes):
            #co_occur[i][j] = gene_sample_set[i].intersection(gene_sample[j]).__len__()
            set1 = gene_sample_set[gi].intersection(samples)
            set2 = gene_sample_set[gj].intersection(samples)
            co[i,j] = len(set1.intersection(set2))
            #co[i,j] = gene_sample_set[gi].intersection(gene_sample[gj]).__len__()
    #json.dump(co_occur, open(outputfile, 'w'))
    np.save(open(outputfile, 'w'), co)
    return co


if 'co_occur_np' in os.listdir('./'):
    co = np.load(open('co_occur_np'))
else:
    co = run_co_occur(gene_sample_set)


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