#  make hapmap tables

import db_importer
import os
import re

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
	db_importer.make_table(filename = f, tablename = 'hapmap_raw', db = db)
