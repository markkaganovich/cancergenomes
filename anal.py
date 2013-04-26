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


db = create_engine('sqlite:///tcga_somatic.db', echo = False)

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select().execute()
f = filter(lambda x: x.variant_classification in first_pass_mutations, e)

def count_snps(results_list):
    snpcount = {}
    for x in results_list:
        snp = x.chrom+':'+x.start_position 
        if snp in snpcount.keys():
            snpcount[snp] = snpcount[snp] + 1
        else:
            snpcount[snp]= 0
    return snpcount

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
