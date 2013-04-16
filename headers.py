#handle headers

'''
1) need function to figure out which line is the header, if there is one
2) make header the fieldnames, else need input for fieldnames

'''

import csv

def find_fieldnames(infile):
	reader = open(infile, 'r')
	h = reader.next().strip('\n')
	while(h.startswith('#')):
		h = h[1:]
	dialect = csv.Sniffer().sniff(h)
	delim = dialect.delimiter
	fieldnames = h.split(delim)

	return fieldnames, delim


def synonyms(fieldname):
	try:
		col = syn[fieldname.lower()]
	except KeyError:
		col = fieldname.lower()
	return col

def get_tablename(filename):
	try:
		tn = file_to_table[filename]
	except KeyError:
		tn = filename
	return tn

###############################################################################
syn = {'chromosome': 'chrom',
		'chr':'chrom',
		'IID' : 'cell_line_id',
		'id' : 'rs#'   # this is going to be an issue: ID fieldname messes with 'id' column name 
		}

file_to_table = {'draft2_samples_QC+.txt': 'hapmap_draft2_samples'}