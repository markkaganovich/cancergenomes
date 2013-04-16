# import 1KG low cov data

import db_importer
import os
import re
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table

files = ['../1000GenomesData/CEU.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/YRI.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/CHBJPT.low_coverage.2010_09.genotypes.vcf', '../1000GenomesData/YRI.trio.2010_09.genotypes.vcf', '../1000GenomesData/CEU.trio.2010_09.genotypes.vcf']

b = create_engine('sqlite:///GENOTYPES.db', echo = True)

for f in files:
    print f
    db_importer.import_data(filename = f, tablename = '1kg_lowcov_1.0', db = db)