#random_importer.py

from sqlalchemy import *
import os
import re

import headers
import db_importer

#db = create_engine('sqlite:///tcga_somatic.db', echo = True)
db = create_engine('sqlite:///tcga_somatic.db', echo = True)

directory = 'mark/cancerdata/'
filename = directory +'refseq_frames.txt'
#cancer = 'BRCA'
#f = directory + filename

#primary_keys = ['chrom', 'start_position', 'tumor_sample_barcode']  
primary_keys = ['hugo_symbol', 'exoncount']

db_importer.import_data(filename = filename, tablename = 'refseq_frames', db = db, key_columns = primary_keys)

def make_protein_sequence(fastafilename):
    f = open(fastafilename)
    protein_sequence = {}
    for l in f:
        if l.startswith('>'):
            gene = re.search('GN=[a-zA-Z0-9]*', l)
            if gene is not None:
                gene = gene.group().split('=')[1]
                protein_sequence[gene] = ''
        else:
            if gene is not None:
                protein_sequence[gene] = protein_sequence[gene] + l.strip('\n')
    return protein_sequence


