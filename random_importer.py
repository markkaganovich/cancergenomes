#random_importer.py

from sqlalchemy import *
import os
import re

import headers
import db_importer

db = create_engine('sqlite:///tcga_somatic.db', echo = True)

directory = 'mark/cancerdata/'
filename = 'genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.3.2.0.somatic.maf'
cancer = 'BRCA'
f = directory + filename

primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
add_columns = {'cancer_type': cancer, 'filename' : filename}
db_importer.import_data(filename = f, tablename = 'mutations_v1', db = db, extra_columns = add_columns, key_columns = primary_keys)
