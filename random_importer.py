#random_importer.py

from sqlalchemy import *
import os
import re

import headers
import db_importer

db = create_engine('sqlite:///tcga_somatic.db', echo = True)

directory = 'mark/cancerdata/'
filename = 'COLON.illumina.maf'
cancer = 'COLON'
f = directory + filename

primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
add_columns = {'cancer_type': cancer, 'filename' : filename}
db_importer.import_data(filename = f, tablename = 'mutations_colon', db = db, extra_columns = add_columns, key_columns = primary_keys)
