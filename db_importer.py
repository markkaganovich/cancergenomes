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

def make_table(tablename = None, db = None, columns = [], key_columns = None):

	metadata = MetaData(db)
	Session = sessionmaker(db)
	session = Session()

	print columns

	if db.dialect.has_table(db.connect(), tablename):
		table = Table(tablename, metadata, autoload = True)
	else:
		table = Table(tablename, metadata, 
	        	Column('id', Integer, primary_key=True),
	            *(Column(rowname, String(), primary_key = get_key_columns(rowname, key_columns)) for rowname in columns))
		table.create()
	session.commit()

	return table

def get_key_columns(colname, key_columns):
	if key_columns is None:
		return False
	if colname in key_columns:
		return True
	else:
		return False


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

	for row in reader:
		print row
		if extra_columns is not None:
			row.update(extra_columns)
		print('\n')
		print values(**row)
		table.insert().values(**row).execute()

	session.commit()

	class Genotypes(object):
		pass
	mapper(Genotypes, table)




