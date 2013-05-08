#random_importer.py

from sqlalchemy import *
import os
import re

import headers
import db_importer

db = create_engine('sqlite:///tcga_somatic.db', echo = True)

directory = 'mark/cancerdata/'
filename = directory +'knownGene_join'
#cancer = 'BRCA'
#f = directory + filename

#primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
primary_keys = ['hg19.knownGene.name']

db_importer.import_data(filename = filename, tablename = 'knownGene', db = db, key_columns = primary_keys)
