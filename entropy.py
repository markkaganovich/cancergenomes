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

# caculate entropty
genes = gene_snp.keys()
entropy = {}
for g in genes:
    v = counts[g].values()
    entropy[g] = -1 * sum(map(lambda x: (float(x)/sum(v)) * np.log(float(x)/sum(v)), v))

entropy_normalized = {}
for g in genes:
    if sum(counts[g].values()) > 30:
        entropy_normalized[g] = entropy[g] / len(counts[g].values())

hack = {}
for g in genes:
    v = np.array(counts[g].values())
    hack[g] = np.dot(v,v)-float(sum(v))

a1 = sorted(entropy_normalized.iteritems(), key = operator.itemgetter(1), reverse=True)
ents1 = map(lambda x: x[1], a)
d1 = map(lambda x: x[0], a)
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

a = sorted(hack.iteritems(), key = operator.itemgetter(1), reverse=True)
ents = map(lambda x: x[1], a)
d = map(lambda x: x[0], a)

distances = {}

