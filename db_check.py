from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import and_, or_
import sys
import json
import commands
import os
import operator

db = create_engine('sqlite:///tcga.db', echo = False)
metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

db2 = create_engine('sqlite:///tcga_somatic.db', echo = False)
metadata2 = MetaData(db2)
Session2 = sessionmaker(db2)
session2 = Session2()

Mutations1 = Table('Mutations', metadata, autoload = True)
Mutations2 = Table('mutations_v1', metadata2, autoload = True)
m1 = Mutations1.select().execute()
rows1 = map(lambda x: x, m1)

m2 = Mutations2.select().execute()
rows2 = map(lambda x: x, m2)

rows2_unique = set(map(lambda x: x.tumor_sample_barcode + ':' + x.chrom + ':' + x.start_position, rows2 ))
where_at = set([])
for r in rows1:
    obj = r.Tumor_Sample_Barcode +':'+r.Chrom+':'+str(r.start_position)
    if obj not in rows2_unique:
        where_at.add(r)

json.dump(list(where_at), open('db_differences', 'w'))

