#  make hapmap tables

import db_importer
import os
import re
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table

directory = 'mark/genotypes/'
matched_files = map(lambda x: re.match(r'hapmapchr+(\d*\Z)', x), os.listdir(directory))
files = []
for m in matched_files:
	if m:
		f = directory+m.group()
		files.append(f)


db = create_engine('sqlite:///GENOTYPES.db', echo = True)

for f in files:
	print f
	db_importer.import_data(filename = f, tablename = 'hapmap_raw_1.0', db = db)