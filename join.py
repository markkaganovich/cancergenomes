#join.py

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
import numpy as np


db = create_engine('sqlite:///tcga_somatic.db', echo = False)

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

metadata = MetaData(db)
Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select().execute()
snps = filter(lambda x: x.variant_classification in first_pass_mutations, m)

db = create_engine('sqlite:///polyphen-2.2.2-whess-2011_12.sqlite', echo = False)

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

features = Table('features', metadata, autoload = True)
bp_to_aa = {}

nones = []
for s in snps:
    print s.hugo_symbol
    try:
        q = session.query(features).filter(and_(features.c.chrom == 'chr'+s.chrom, features.c.chrpos == int(s.start_position))).first()
    except ValueError:
        continue
    if q is None:
        nones.append(q)
        print "None  " +s.chrom+':'+s.start_position
    else:
        bp_to_aa[s.chrom+':'+s.start_position] = {'pos' : q.pos, 'aa': q.aa1}

json.dump(bp_to_aa, open('bp_to_aa', 'w'))

