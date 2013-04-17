# tcga importer

from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import and_, or_
import sys
import json
import commands
import os
import operator
import re

import headers
import db_importer

db = create_engine('sqlite:///tcga_somatic.db', echo = False)

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

directory = 'mark/cancerdata/'
matched_files = map(lambda x: re.match(r'[a-zA-Z.]*maf(\Z)', x), os.listdir(directory))
files = []
for m in matched_files:
    if m:
        f = directory+m.group()
        files.append(f)
        cancer = re.search('[A-Z]+', f)    
        db_importer.import_data(filename = f, tablename = 'mutations_v1', db = db, extra_columns = {'cancer_type': cancer, 'filename' : m.group()})

headers.explore_headers(files)


