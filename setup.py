'''
TCGA Somatic Variation setup

mark, Aug 2013
'''

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import commands
import os
import csv 

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

syn = {'chrom' : ['chromosome', 'chr'],
	'cell_line_id' : ['iid'],
	'rsid' : ['id'],
	'hugo_symbol': ['name2']
}

def synonyms(fieldname):
	'''
	return key in syn table that is synonyms with the given fieldname
	'''
	col = None
	for k in syn.keys():
		if fieldname.lower() in syn[k]:
			col = k
			break
	if col is None:
		col = fieldname.lower()

	return col

file_to_table = {'draft2_samples_QC+.txt': 'hapmap_draft2_samples'}

def get_tablename(filename):
	try:
		tn = file_to_table[filename]
	except KeyError:
		tn = filename
	return tn	

def find_fieldnames(infile, not_header_flag = None):
	reader = open(infile, 'r')
	h = reader.next().strip('\n')
	if not_header_flag is not None:
		while(h.startswith(not_header_flag)): # for example '#'
			h = h[1:]
	dialect = csv.Sniffer().sniff(h)
	delim = dialect.delimiter
	fieldnames = h.split(delim)

	for f in range(0, len(fieldnames)):
		if fieldnames[f].startswith('#'):
			fieldnames[f] = fieldnames[f][1:]

	return fieldnames, delim


def get_key_columns(colname, key_columns):
	if key_columns is None:
		return False
	if colname in key_columns:
		return True
	else:
		return False

def make_table(tablename = None, db = None, columns = [], key_columns = None):

	metadata = MetaData(db)
	Session = sessionmaker(db)
	session = Session()

	if db.dialect.has_table(db.connect(), tablename):
		table = Table(tablename, metadata, autoload = True)
	else:
		table = Table(tablename, metadata, 
	        	#Column('id', Integer, primary_key=True),
	            *(Column(rowname, String(), primary_key = get_key_columns(rowname, key_columns)) for rowname in columns))
		table.create()
	session.commit()

	return table


def import_data(filename = 'testheader2', tablename = None, db = None, extra_columns = None, key_columns = None):

	Session = sessionmaker(db)
	session = Session()

	fieldnames, delim = headers.find_fieldnames(filename)
	
	#check that fieldnames are desired column names
	columns = map(lambda x: headers.synonyms(x), fieldnames)
	if extra_columns is not None:
		columns.extend(extra_columns.keys())

	if tablename is None:
		tablename = headers.get_tablename(filename)

	csvfile = open(filename, 'r')
	reader = csv.DictReader(csvfile, fieldnames = columns, delimiter = delim)

	table = make_table(tablename, db, columns, key_columns)

	print columns

	for row in reader:
		print row
		print row.values()
		if extra_columns is not None:
			row.update(extra_columns)
		table.insert(prefixes=['OR IGNORE']).values(**row).execute()

	session.commit()


if __name__ == "__main__":

	db = create_engine('sqlite:///tcga_somatic.db', echo = False)
	directory = 'mark/cancerdata/'
	matched_files = map(lambda x: re.match(r'[a-zA-Z.]*maf(\Z)', x), os.listdir(directory))
	files = []
	for m in matched_files:
    	if m:
        	f = directory+m.group()
        	files.append(f)
        	cancer = re.search('[A-Z]+', f).group() 
        	primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
        	add_columns = {'cancer_type': cancer, 'filename' : m.group()}
        	import_data(filename = f, tablename = 'mutations_v1', db = db, extra_columns = add_columns, key_col
umns = primary_keys)

# convert to hg19
# db_cleanup.py


