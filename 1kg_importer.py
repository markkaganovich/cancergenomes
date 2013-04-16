# import 1KG low cov data

import db_importer
import os
import re
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table

files = ['../1000GenomesData/CEU.lowcov.header', '../1000GenomesData/YRI.lowcov.header', '../1000GenomesData/CHBJPT.lowcov.header', '../1000GenomesData/YRI.trio.header', '../1000GenomesData/CEU.trio.header']

db = create_engine('sqlite:///GENOTYPES.db', echo = True)

for f in files:
    print f
    db_importer.import_data(filename = f, tablename = '1kg_lowcov_1.0', db = db)