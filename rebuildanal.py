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











