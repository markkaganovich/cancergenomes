from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import operator
from pysqlite2 import dbapi2 as sqlite
import rebuildanal_helper
from collections import Counter

import Queue
import logging
import simplejson


db = create_engine('sqlite:///tcga_somatic.db', echo = False)
metadata = MetaData(db)

mutations_table = Table('mutations_v1', metadata, autoload = True)
mutation_entries= mutations_table.select(mutations_table.c.hugo_symbol != 'Unknown').execute().fetchall() 

# make it into an explicit dictionary
variants = []
for m in mutation_entries:
	variants.append(tuple([m['chrom'], m['start_position'], m['dbsnp_rs'], m['reference_allele'], tuple(set([m['tumor_seq_allele1'], m['tumor_seq_allele2']]))]))


###
# do it from json files
###
def import_from_json(file):
	f = open(file)
	variants = []
	lines = f.readlines()
	for l in lines:
		m = simplejson.loads(l)
		variants.append(tuple([m['chromosome'], m['start'], m['dbsnp_rs'], m['reference_allele'], tuple(set([m['tumor_seq_allele1'], m['tumor_seq_allele2']]))]))
	return variants

v1 = import_from_json('somatic_mutations1.json')
v2 = import_from_json('somatic_mutations2.json')
v3 = import_from_json('somatic_mutations3.json')
v4 = import_from_json('somatic_mutations4.json')
variants = v1+v2+v3+v4
#variants = sorted(variants, key = lambda x: (x[0], x[1]))

vcf_format = sorted(list(set(variants)), key = lambda x: (x[0], x[1]))

variants_json = []
for variant in vcf_format:
	d = {'chrom': variant[0], 'pos': variant[1], 'id': variant[2], 'ref': variant[3], 'alt': list(set(variant[4]))}
	variants_json.append(d)
	

def print_vcf(variants_json, file):
	output = open(file, 'w')
	for variant in variants_json:
			alt = ''
			for i in variant['alt']:
				alt = alt + i + ','
			alt = alt.strip(',')
			output.write(variant['chrom'] + '\t' + str(variant['pos']) + '\t' + variant['id'] + '\t' + variant['ref'] + '\t' + alt + '\n')


print_vcf(variants_json, 'tcga.vcf')



