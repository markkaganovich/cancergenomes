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

db = create_engine('sqlite:///tcga_somatic.db', echo = False)
metadata = MetaData(db)
#Session = sessionmaker(db)
#session = Session()

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select(Mutations.c.hugo_symbol != 'Unknown').execute().fetchall()   # filter out IGR variants; those that don't map to known genes

class Rows(object):
    def __init__(self,object):
        for x in object.keys():
            setattr(self, x, object[x])

rows = map(lambda x: Rows(x), m)
gene_names = list(set(map(lambda x: x.hugo_symbol, rows)))    

polyphen = create_engine('sqlite:///polyphen-2.2.2-whess-2011_12.sqlite', echo = False)
polyphen_metadata = MetaData(polyphen)
polyphen_features = Table('features', polyphen_metadata, autoload = True)

for gene in gene_names:
	gene_rows = filter(lambda x: x.hugo_symbol == gene, rows)
	gene_features = polyphen_features.select(and_(polyphen_features.c.gene == gene, polyphen_features.c.based_on == 'alignment')).execute().fetchall()   # filter only entries in polyphen that come from alignment not alignment_mz
	for gr in gene_rows:
		filter(lambda x: x.chrpos == int(gr.start_position), gene_features)

conn = sqlite.connect("SIFT/Human_CHR1.sqlite")
tables = cursor.execute("select name from sqlite_master where type = 'table'").fetchall()

sorted_rows = sorted(rows, key=lambda x: x.chrom)

def dispatch(chrom, pos):
	pos = int(pos)
	chrom = chrom.strip(' ')
	db_name = 'SIFT/Human_CHR' + chrom + '.sqlite'
	conn = sqlite.connect(db_name)
	cursor = conn.cursor()
	tables = map(lambda x: x[0], cursor.execute("select name from sqlite_master where type = 'table'").fetchall())
	table_ranges = map(lambda x: (int(x.split('_')[1]), int(x.split('_')[2])), tables)
	select_range = filter(lambda x: pos >= x[0] and pos <= x[1], table_ranges)[0]
	select_table = 'chr' + chrom + '_' + str(select_range[0])+ '_' + str(select_range[1])
	sift_row = cursor.execute("select AAPOS2 from " + select_table + "  where CHR == ? and COORD2 == ?", ['chr'+chrom, str(pos)]) 
	return sift_row









