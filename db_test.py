#db_test.py

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table

db = create_engine('sqlite:///GENOTYPES.db', echo = True)

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

tablename = "hapmap_raw_1.0"
table = Table(tablename, metadata, autoload = True)
