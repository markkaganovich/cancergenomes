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

def make_table(tablename = None, db = None, columns = []):

	metadata = MetaData(db)
	Session = sessionmaker(db)
	session = Session()

	if db.dialect.has_table(db.connect(), tablename):
		table = Table(tablename, metadata, autoload = True)
	else:
		table = Table(tablename, metadata, 
	        	Column('id', Integer, primary_key=True),
	            *(Column(rowname, String()) for rowname in columns))
		table.create()
	session.commit()

	return table



def import_data(filename = 'testheader2', tablename = None, db = None):

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

	table = make_table(tablename, db, columns)

	for row in reader:
		table.insert().values(**row).execute()

	session.commit()

	class Genotypes(object):
		pass
	mapper(Genotypes, table)




