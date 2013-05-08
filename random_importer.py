#random_importer.py

from sqlalchemy import *
import os
import re

import headers
import db_importer

db = create_engine('sqlite:///tcga_somatic.db', echo = True)

directory = 'mark/cancerdata/'
filename = 'knownGene_join'
#cancer = 'BRCA'
#f = directory + filename

#primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
primary_keys = ['hg19.knownGene.name']

db_importer.import_data(filename = f, tablename = 'mutations_v1', db = db, extra_columns = add_columns, key_columns = primary_keys)
