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


db = create_engine('sqlite:///tcga_somatic.db', echo = False)

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select().execute()
rows = filter(lambda x: x.variant_classification in first_pass_mutations, m)

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

def explore_snp_counts(snpcount):
    a = sorted(snpcount.iteritems(), key = operator.itemgetter(1), reverse=True)
    counts = map(lambda x: x[1], a)[0:50000]
    snps = map(lambda x: x[0], a)[0:50000]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ind = list(range(0, len(counts)))
    width = .2
    print "making bar"
    bar = ax.bar(ind, counts, width, color="blue")
    print "saving ..."
    plt.savefig('counts.png')

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
        


snp_gene = dict(map(lambda x: (x.chrom +':' + x.start_position, x.hugo_symbol), rows))


def make_matrix(outputfile = 'genotype_matrix.temp', snpcountfile = 'snpcount.temp'):
    out = open(outputfile, 'w')
    snpcountout = open(snpcountfile, 'w')

    sam = session.query(Mutations).all()
    a = list(set(filter(lambda x: x.Variant_Classification != 'Silent', sam)))
    allsnps = list(set(map(lambda x: str(x.Chromosome) + ':' + str(x.Start_Position), a)))
    samples = list(set(map(lambda x: x.Tumor_Sample_Barcode, a)))

    line = 'SAMPLE,'
    snpcount = {}
    for snp in allsnps:
        line = line + snp + ','
        snpcount[snp] = 0
    out.write(line.strip(',')+'\n')

    for s in samples:
        try:
            out.write(s + ',')
            rows = filter(lambda x: x.Tumor_Sample_Barcode == s, a)
            rowsnps = set(map(lambda x: str(x.Chromosome) + ':' + str(x.Start_Position), rows))
            snps = ''
            for snp in allsnps:
                if snp in rowsnps:
                    snps = snps + '1' + ','
                    snpcount[snp] = snpcount[snp] + 1
                else:
                    snps = snps + '0' + ','
                    snpcount[snp] = snpcount[snp] + 0
            out.write(snps.strip(',')+'\n')
        except:
            continue
    json.dump(snpcount, snpcountout)
