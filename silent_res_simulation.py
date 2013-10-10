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
import Queue
import threading
import logging

db = create_engine('sqlite:///tcga_somatic.db', echo = False)
metadata = MetaData(db)

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select(Mutations.c.hugo_symbol != 'Unknown').execute().fetchall()   # filter out IGR variants; those that don't map to known genes

class Rows(object):
    def __init__(self,object):
        for x in object.keys():
            setattr(self, x, object[x])

tcga_silent = filter(lambda x: x.variant_classification == 'Silent')
genes = list(set(map(lambda x: x.hugo_symbol, tcga_silent)))


if 'snp_to_aa' in os.listdir('./'):
	snp_to_aa = json.load(open('snp_to_aa'))
else:	
	snp_to_aa = rebuildanal_helper.get_snp_to_aa(tcga_silent)

snps = set(snp_to_aa.keys())
tcga_residues = []
for r in tcga_silent:
	if r.chrom + ':' + r.start_position in snps:
		r.residue = snp_to_aa[r.chrom + ':' + r.start_position]
		tcga_residues.append(r)

if 'counts_aa_silent' in os.listdir('./'):
	counts_aa= json.load(open('counts_aa_silent'))
else:
	counts_aa  = rebuildanal_helper.make_counts_aa(tcga_residues)
	json.dump(counts_aa, open('counts_aa_silent'))



def sim(gene):
	peaks_max = []
	sim = 0
	while sim < 1E4:
		gene_muts = sum(counts_aa[gene].values())
		pos = []
		for p in range(0, gene_muts):
			pos.append(random.randint(0, prtn_len[gene]))
		sim +=1
		counts = Counter(pos)
		peaks_max.append(max(counts.values()))	
	std = numpy.std(peaks_max)
	mean_peak = numpy.mean(peaks_max)
	metric = (numpy.max(counts_aa[gene].values()) - mean_peak) / std
	return mean_peak, std, metric


logging.basicConfig(filename='simulations_silent.log',level=logging.DEBUG)

prtn_len = json.load(open('prtn_len'))

def do_work(item):
	print "Worker running: %s" % item
	result = sim(item)
	#peak_stds[item] = result
	logging.info('\t' + str(item) + ' \t ' + str(result[0]) + '\t' + str(result[1]) + '\t' + str(result[2]))
	return result
	

def worker():
    while True:
        item = q.get()
        do_work(item)	
        q.task_done()

q = Queue.Queue()
for i in range(10):
     t = threading.Thread(target=worker)
     t.daemon = True
     t.start()

#sim7genes = json.load(open('sim7genes'))
#sim8genes = json.load(open('sim8genes'))
for gene in genes:
	if gene in counts_aa.keys() and sum(counts_aa[gene].values()) > 4 and gene in prtn_len.keys():
		print "Queuing %s" % gene
		q.put(gene)

q.join()




