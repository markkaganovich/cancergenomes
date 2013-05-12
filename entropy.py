#calculate entropy

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
from collections import Counter
import numpy as np
import math


db = create_engine('sqlite:///tcga_somatic.db', echo = False)

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select().execute()
rows = filter(lambda x: x.variant_classification in first_pass_mutations, m)

def make_gene_snp(rows):
    gene_snp = {}
    for r in rows:
        gene = r.hugo_symbol
        snp = r.chrom + ':' + r.start_position
        try:
            gene_snp[gene].append(snp)
        except KeyError:
            gene_snp[gene] = [snp]
    return gene_snp

gene_snp = make_gene_snp(rows)
counts = {}
for g in gene_snp.keys():
    counts[g] = Counter(gene_snp[g])

def poisson(k, l):
    return math.pow(l, k)/math.factorial(k) * np.exp(-1*l)


# caculate entropty
gene_sequences = json.load(open('gene_sequences'))
genes = gene_snp.keys()
entropy_possion = {}
for g in genes:
    v = counts[g].values()
    #v = map(lambda x: long(x), v)
    try:
        #entropy[g] = -1 * sum(map(lambda x: x * float(len(v))/len(gene_sequences[g]) * np.log(x* (float(len(v))/len(gene_sequences[g]))), v))
        #entropy[g] = -1 * float(sum(map(lambda x: np.log(x), v))) / (len(gene_sequences[g])*2) 
        
        total = []
        for k in v:
            l = float(sum(v)-k)/ len(gene_sequences[g])
            try:
                p = poisson(k, l)
            except OverflowError:
                p = 0
            total.append(p)
        entropy_possion[g] =  product(total)
    except KeyError:
        continue

a1 = sorted(entropy_possion.iteritems(), key = operator.itemgetter(1))
ents1 = map(lambda x: x[1], a1)
d1 = map(lambda x: x[0], a1)
'''
entropy_normalized = {}
for g in genes:
    if sum(counts[g].values()) > 30:
        entropy_normalized[g] = entropy[g] / sum(counts[g].values())
'''
hack = {}
hack_max = {}
keys = set(gene_sequences.keys())
for g in genes:
    l = counts[g].values()    
    v = np.array(counts[g].values())
    hack_max[g] = max(v/np.mean(v))   
    

a = sorted(hack_max.iteritems(), key = operator.itemgetter(1), reverse=True)
ents = map(lambda x: x[1], a)
d = map(lambda x: x[0], a)


'''
fig = plt.figure()
ax = fig.add_subplot(111)
ind = list(range(0, len(ents)))
width = .2
print "making bar"
bar = ax.bar(ind, ents, width, color="r")
print "saving ..."
plt.savefig('hack.png')
'''
#if g in keys:
#hack[g] = np.dot(v,v)/float(sum(v))/len(gene_sequences[g])


m_mut = sorted(counts.iteritems(), key = lambda x: sum(x[1].values()), reverse=True)
r_mut = map(lambda x: x[0], m_mut)
distances = {}
for g in genes:
    distances[g] = r_mut.index(g) - d.index(g) 
e = sorted(distances.iteritems(), key = operator.itemgetter(1), reverse=True)


def entr(v):
    return sum(map(lambda x: (float(x)/1000) * np.log(float(x)/1000), v))
