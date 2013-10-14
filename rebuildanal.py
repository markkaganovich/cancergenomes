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
from collections import Counter
import Queue
import threading
import logging
import scipy
from scipy import stats

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

tcga_nonsilent = filter(lambda x: x.variant_classification != 'Silent', tcga_rows)

tcga_silent = filter(lambda x: x.variant_classification == 'Silent', tcga_rows)

if 'snp_to_aa' in os.listdir('./'):
	snp_to_aa = json.load(open('snp_to_aa'))
else:	
	snp_to_aa = rebuildanal_helper.get_snp_to_aa(tcga_rows)

snps = set(snp_to_aa.keys())
tcga_residues = []
for r in tcga_nonsilent:
	if r.chrom + ':' + r.start_position in snps:
		r.residue = snp_to_aa[r.chrom + ':' + r.start_position]
		tcga_residues.append(r)

if 'counts_aa' in os.listdir('./'):
	counts_aa = json.load(open('counts_aa'))
else:
	counts_aa = rebuildanal_helper.make_counts_aa(tcga_residues)
	json.dump(counts_aa, open('counts_aa', 'w'))


residues = list(set(map(lambda x: x.residue, tcga_residues)))

tcga_residues_silent = []
for r in tcga_silent:
	if r.chrom + ':' + r.start_position in snps:
		r.residue = snp_to_aa[r.chrom + ':' + r.start_position]
		tcga_residues_silent.append(r)

silent_residues = list(set(map(lambda x: x.residue, tcga_residues_silent)))

'''
suppr_types = map(lambda x: x.variant_classification, filter(lambda y: y.hugo_symbol in suppressors, tcga_rows))
suppr_types[0]
suppr_types[1]
from collections import Counter
suppr_counter = Counter(suppr_types)
suppr_counter.keuys()
suppr_counter.keys()
suppr_counter
onco_types = map(lambda x: x.variant_classification, filter(lambda y: y.hugo_symbol in oncogenes, tcga_rows))
onco_counter = Counter(onco_types)
'''

## read in simulations.log result

def sim_output(simulation_file = 'simulations.log'):

	sim_gene_results = {}
	lines = open(simulation_file).readlines()
	simgenes = map(lambda x: x.split('\t')[1].strip('  '), lines)
	for gene in simgenes:
		sim_gene_results[gene] = {}

	for l in lines:
		tabs = l.split('\t')
		g = tabs[1].strip('  ')
		sim_gene_results[g]= {}
		sim_gene_results[g]['mean'] = float(tabs[2])
		sim_gene_results[g]['std'] = float(tabs[3])
		sim_gene_results[g]['metric'] = float(tabs[4].strip('\n'))

	return sim_gene_results

## per residue

def pile_up(sim_gene_results, residues, counts_aa):
	pile_ups = {}
	what = []
	for res in residues:
		gene = res.split(':')[0]
		try:
			mean = sim_gene_results[gene]['mean']
			std = sim_gene_results[gene]['std']
		except KeyError:
			continue
		if res not in counts_aa[gene].keys():
			what.append(res)
		else:
			val = counts_aa[gene][res]
			pile_ups[res] = (val-mean)/std

	return pile_ups


sim_gene_results = sim_output('simulations.log')
pile_ups = pile_up(sim_gene_results, residues, counts_aa)
pileupsorted = sorted(pile_ups.iteritems(), key=operator.itemgetter(1), reverse=True)


counts_aa_silent = json.load(open('counts_aa_silent'))
sim_gene_results_silent = sim_output('simulations_silent.log')
pile_up_silent = pile_up(sim_gene_results_silent, silent_residues, counts_aa_silent) 

##########################################################################
# cancer type distribution

cancer_types = map(lambda x: x.cancer_type, tcga_residues)
cancer_types_silent = map(lambda x: x.cancer_type, tcga_residues_silent)
cancer_types_counts = Counter(cancer_types)
cancer_types_silent_count = Counter(cancer_types_silent)
map(lambda x: (x, cancer_types_counts[x]/float(sum(cancer_types_counts.values()))), cancer_types_counts.keys())
map(lambda x: (x, cancer_types_silent_count[x]/float(sum(cancer_types_silent_count.values()))), cancer_types_silent_count.keys())

cancers = cancer_types_counts.keys()

'''
total = float(sum(cpc.values()))
distr = []
for c in cancers:	
	if c in cpc.keys():
		distr.append(cpc[c]/total)
	else:
		distr.append(0.0)
'''


def get_np_array(cancers, cancer_counts):
	distr = []
	for c in cancers:	
		if c in cancer_counts.keys():
			distr.append(float(cancer_counts[c]))
		else:
			distr.append(0.0)
	distr = np.array(distr)
	#distr = distr/sum(distr)
	return distr

overall = get_np_array(cancers, cancer_types_counts)
silent_distr = get_np_array(cancers, cancer_types_silent_count)

scipy.stats.ks_2samp(overall, silent_distr)


expected_freq = overall / sum(overall) 
def peak_chisq(peak, tcga_residues, expected_freq):
	peak_rows = filter(lambda x: x.residue == peak[0], tcga_residues)
	peak_cancers = Counter(map(lambda x: x.cancer_type, peak_rows))
	distr = get_np_array(cancers, peak_cancers) 
	chisq = scipy.stats.chisquare(distr, expected_freq*sum(distr))
	return chisq

#distr_ls = []
'''
chisq = []
for peak in pileupsorted:
	peak_rows = filter(lambda x: x.residue == peak[0], tcga_residues)
	peak_cancers = Counter(map(lambda x: x.cancer_type, peak_rows))
	distr = get_np_array(cancers, peak_cancers) 
	chisq.append(scipy.stats.chisquare(distr, expected_freq*sum(distr)))
'''


logging.basicConfig(filename='pile_up_chisq.log',level=logging.DEBUG)

def do_work(peak):
	print "Worker running: %s" % item
	result = peak_chisq(peak,tcga_residues,expected_freq)
	logging.info('\t' + str(item) + ' \t ' + str(result[0]) + '\t' + str(result[1]))
	return result
	

def worker():
    while True:
        peak = q.get()
        do_work(peak)	
        q.task_done()

q = Queue.Queue()
for i in range(14):
     t = threading.Thread(target=worker)
     t.daemon = True
     t.start()

for peak in pileupsorted:
	print "Queuing %s" % peak[0]
	q.put(peak)

q.join()





