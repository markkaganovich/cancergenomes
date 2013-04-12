'''
import data into databases
file, fieldnames
make a table

'''
from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import commands
import os
import csv 
import headers

#db = create_engine('sqlite:///GENOTYPES.db', echo = True)

def make_table(filename = 'testheader2', tablename = None, db = None):

	metadata = MetaData(db)
	Session = sessionmaker(db)
	session = Session()

	print metadata.tables.keys()

	fieldnames, delim = headers.find_fieldnames(filename)
	#check that fieldnames are desired column names
	columns = map(lambda x: headers.synonyms(x), fieldnames)

	if tablename is None:
		tablename = headers.get_tablename(filename)

	csvfile = open(filename, 'r')
	reader = csv.DictReader(csvfile, fieldnames = columns, delimiter = delim)

	if tablename in metadata.tables.keys():
		#table = Table(tablename, metadata, autoload = True)
		table = metadata.tables[tablename]
	else:
		table = None
	
	print "here"

	for row in reader:
		if table is None:
			#create the table
			print metadata.tables.keys()
			table = Table(tablename, metadata, 
	        	Column('id', Integer, primary_key=True),
	            *(Column(rowname, String()) for rowname in columns))
			table.create()
		else:
			table.insert().values(**row).execute()

	session.commit()

	class Genotypes(object):
		pass
	mapper(Genotypes, table)




