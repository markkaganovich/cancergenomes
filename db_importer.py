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

def make_table(filename = 'testheader2'):

	db = create_engine('sqlite:///GENOTYPES.db', echo = True)
	metadata = MetaData(db)
	Session = sessionmaker(db)
	session = Session()	

	fieldnames, delim = headers.find_fieldnames(filename)
	#check that fieldnames are desired column names
	columns = map(lambda x: headers.synonyms(x), fieldnames)

	try:
		tablename = headers.file_to_table[filename]
	except KeyError:
		tablename = 'test'
		print "Need table name, change this to assert when this becomes a function"

	csvfile = open(filename, 'r')
	reader = csv.DictReader(csvfile, fieldnames = fieldnames, delimiter = delim)

	table = None
	for row in reader:
		if table is None:
			#create the table
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




