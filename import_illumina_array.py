#import Illumina array data

from sqlalchemy import *
import os
import re

import headers
import db_importer

directory = 'mark/arraydata/'
arrayfile = 'Array1S_finalreport'

db = create_engine('sqlite:///mark/arraydata/arrays.db', echo = False)
db_importer.import_data(filename = directory+arrayfile, tablename = arrayfile, db = db)