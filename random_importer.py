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
    prtn_len = {}
    for l in f:
        if l.startswith('>'):
            gene = re.search('geneSymbol=[a-zA-Z0-9]*', l)
            if gene is not None:
                gene = gene.group().split('=')[1]
                length = re.search('Len=[0-9]*', l)
                if length is not None:
                    length = length.group().split('=')[1]
                    prtn_len[gene] = length
            '''
            t = l.split(' ')
            name = t[0].split('|')[2].split('_')[0]
            pl = re.search('Length:\s[0-9]*', l)
            if pl is not None:
                length = int(pl.group().split(':')[1]) 
                prtn_len[name] = length
            '''
    return prtn_len

'''

            if gene is not None:
                gene = gene.group().split('=')[1]
                protein_sequence[gene] = ''
        else:
            if gene is not None:
                protein_sequence[gene] = protein_sequence[gene] + l.strip('\n')
    return protein_sequence
'''



def get_start(starts, frames, snp, strand, stops, rem):
    '''
    starts must be reversed for - strand inputs
    '''
    if len(starts) != 1:
        next_start = starts[1]
        next_frame = frames[1]
        if (strand == '+' and snp >= next_start) or (strand =='-' and snp <= next_start):
            if strand == '+':
                rem = rem + stops[0] - starts[0] + 1
            if strand == '-':
                rem = rem+starts[0] - stops[0] 
            print rem
            return [next_start, next_frame, starts[1:], frames[1:], stops[1:], rem]
    return [starts[0], frames[0], starts, frames, stops, rem]

''' 
do this for every gene: get an rseq thing from the refseq db, then 
consolidate that counts into residues
'''
def codon(rseq, count):
    codon_counts = {}
    starts = map(lambda x: int(x), rseq.exonstarts.strip(',').split(','))
    frames = map(lambda x: int(x), rseq.exonframes.strip(',').split(','))
    stops = map(lambda x: int(x), rseq.exonends.strip(',').split(','))
    snps = map(lambda x: int(x.split(':')[1]), count.keys())
    snps.sort()
    rem = 0
    chrom = count.keys()[0].split(':')[0]
    if rseq.strand == '-':
        temp = starts
        starts = stops
        stops = temp
        starts.reverse()
        frames.reverse()
        stops.reverse()
        snps.sort(reverse=True)
    for snp in snps:
        [start, frame, starts, frames, stops, rem] = get_start(starts, frames, snp, rseq.strand, stops, rem)
        if rseq.strand == '+':
            base = (snp - start + 1 - frame)
        if rseq.strand == '-':
            print start
            print snp
            base = (start - snp )
            print base
        residue = (rem + base) / 3 
        print("residue:%s", residue)
        if base % 3 == 3:
            print "wtf"
        k = str(residue)
        print("k",k)
        if k in codon_counts.keys():
            codon_counts[k] = codon_counts[k] + count[chrom+':'+str(snp)]
            print "double"
            print k
        else:
            codon_counts[k] = count[chrom+':'+str(snp)]
    return codon_counts









