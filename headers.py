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
		col = syn[fieldname]
	except KeyError:
		col = fieldname
	return col


###############################################################################
syn = {'Chromosome': 'Chrom',
		'IID' : 'cell_line_id'
		}

file_to_table = {'draft2_samples_QC+.txt': 'hapmap_draft2_samples'}