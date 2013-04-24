# import 1KG low cov data

import db_importer
import os
import re
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
from sqlalchemy.orm import sessionmaker

files = ['../1000GenomesData/CEU.lowcov.header', '../1000GenomesData/YRI.lowcov.header', '../1000GenomesData/CHBJPT.lowcov.header', '../1000GenomesData/YRI.trio.header', '../1000GenomesData/CEU.trio.header']

def extract_header(files):
    headers = []
    for f in files:
        l = open(f, 'r').readline().split('\t')
        map(lambda x: headers.append(x), l)
    headers = sorted(headers, key=lambda x: x.startswith('N'))
    return list(set(headers))

db = create_engine('sqlite:///GENOTYPES.db', echo = True)
Session = sessionmaker(db)
session = Session()
metadata = MetaData(db)

headers = extract_header(files)

key_columns = ['rsid']
columns = map(lambda x: db_importer.headers.synonyms(x), headers)

kg_table = Table('kg_lowcov', metadata, 
    *(Column(rowname, String(), primary_key = db_importer.get_key_columns(rowname, key_columns)) for rowname in headers))


for filename in files:
    print filename
    csvfile = open(filename, 'r')
    reader = csv.DictReader(csvfile, fieldnames = columns, delimiter = '\t')
    for row in reader:
        if extra_columns is not None:
            row.update(extra_columns)
        table.insert(prefixes=['OR IGNORE']).values(**row).execute()
    session.commit()

