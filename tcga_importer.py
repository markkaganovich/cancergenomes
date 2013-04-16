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

db = create_engine('sqlite:///tcga_2.db', echo = False)

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

headers.explore_headers(files)


