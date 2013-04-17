# tcga importer

from sqlalchemy import *
import os
import re

import headers
import db_importer

db = create_engine('sqlite:///tcga_somatic.db', echo = False)

directory = 'mark/cancerdata/'
matched_files = map(lambda x: re.match(r'[a-zA-Z.]*maf(\Z)', x), os.listdir(directory))
files = []
for m in matched_files:
    if m:
        f = directory+m.group()
        files.append(f)
        cancer= re.search('[A-Z]+', f).group() 
        primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
        add_columns = {'cancer_type': cancer, 'filename' : m.group()}
        db_importer.import_data(filename = f, tablename = 'mutations_v1', db = db, extra_columns = add_columns, key_columns = primary_keys)


#headers.explore_headers(files)

