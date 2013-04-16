#handle headers

'''
1) need function to figure out which line is the header, if there is one
2) make header the fieldnames, else need input for fieldnames

'''

import csv

def find_fieldnames(infile, not_header_flag = None):
	reader = open(infile, 'r')
	h = reader.next().strip('\n')
	if not_header_flag is not None:
		while(h.startswith(not_header_flag)): # for example '#'
			h = h[1:]
	dialect = csv.Sniffer().sniff(h)
	delim = dialect.delimiter
	fieldnames = h.split(delim)

	return fieldnames, delim


def synonyms(fieldname):
	'''
	try:
		col = syn[fieldname.lower()]
	except KeyError:
		col = fieldname.lower()
	return col
	'''
	'''
	return key in syn table that is synonyms with the given fieldname
	'''
	for k in syn.keys():
		if fieldname.lower() in syn[k]:
			col = k
			break
	if col is None:
		col = fieldname.lower()

	return col

def get_tablename(filename):
	try:
		tn = file_to_table[filename]
	except KeyError:
		tn = filename
	return tn

###############################################################################
syn = {'chrom' : ['chromosome', 'chr'],
	'cell_line_id' : ['iid'],
	'rs#' : ['id']
}

file_to_table = {'draft2_samples_QC+.txt': 'hapmap_draft2_samples'}

def explore_headers(files):
	out = open('explore_headers_out', 'w')
	for f in files:
		fieldname, delim = find_fieldnames(f)
		for fn in fieldname:
			out.write(fn + '\t')
		out.write('\n')


