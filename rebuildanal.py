#rebuildanal.py
#
# August 12, 2013
#
from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
from sqlalchemy.sql import and_, or_
import sys
import json
import commands
import os
import operator
import matplotlib.pyplot as plt
import numpy as np
from pysqlite2 import dbapi2 as sqlite
import rebuildanal_helper
import math
import random
import numpy
from collections import Counter

db = create_engine('sqlite:///tcga_somatic.db', echo = False)
metadata = MetaData(db)

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select(Mutations.c.hugo_symbol != 'Unknown').execute().fetchall()   # filter out IGR variants; those that don't map to known genes

class Rows(object):
    def __init__(self,object):
        for x in object.keys():
            setattr(self, x, object[x])

tcga_rows = map(lambda x: Rows(x), m)
genes = list(set(map(lambda x: x.hugo_symbol, tcga_rows)))    

if 'snp_to_aa' in os.listdir('./'):
	snp_to_aa = json.load(open('snp_to_aa'))
else:	
	snp_to_aa = rebuildanal_helper.get_snp_to_aa(tcga_rows)

snps = set(snp_to_aa.keys())
tcga_residues = []
for r in tcga_rows:
	if r.chrom + ':' + r.start_position in snps:
		r.residue = snp_to_aa[r.chrom + ':' + r.start_position]
		tcga_residues.append(r)

if 'counts_aa' in os.listdir('./'):
	counts_aa = json.load(open('counts_aa'))
else:
	counts_aa = rebuildanal_helper.make_counts_aa(tcga_residues)


#########################################################################################
'''
poission analysis of peaks per genes

'''

prtn_len = json.load(open('prtn_len'))

def poisson(k, l):
    try:
        p = math.pow(l, k)/math.factorial(k) * np.exp(-1*l)
    except OverflowError:
        p = 0
    return p


poisson_residues = {}
for g in counts_aa:
    m = sum(counts_aa[g].values())
    try:
        p = float(m) / prtn_len[g]
    except KeyError:
        continue
    for aa in counts_aa[g]:
        if counts_aa[g][aa] > 1:
            poisson_residues[g+':'+aa] = poisson(counts_aa[g][aa], p)


avg = numpy.array([0.0] * prtn_len['BRAF'])

peak_poisson_values = []
while sim < 1000:
	null = []
	sum_mut = 0
	for p in range(1, prtn_len['BRAF']):
   		if sum_mut == sum(counts_aa['BRAF'].values()):
   			break
    	else:
        	mut = random.randint(0, sum(counts_aa['BRAF'].values()) - sum_mut) 
        	sum_mut += mut
        	null.append(mut)

    peak_poisson_values.append(map(lambda x: poisson(x, sum(null)), null))
    sim += 1


sim_peaks = []
peaks_max = []
peaks_mean = []
peaks_sd = []
sim = 0
while sim < 10000:
	gene_muts = sum(counts_aa['BRAF'].values())
	pos = []
	for p in range(1, gene_muts):
		# pick a position randomly
		pos.append(random.randint(0, prtn_len['BRAF']))
	sim +=1
	counts = Counter(pos)
	peaks_max.append(max(counts.values()))	
	#peaks_sd.append(numpy.std(counts.values()))
	#peaks_mean.append(numpy.mean(counts.values()))
s = numpy.std(peaks_max)
m = numpy.mean(peaks_max)

(numpy.max(counts_aa['BRAF'].values()) - m) / s











