from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import and_, or_
import sys
import json
import commands
import os
import operator

db = create_engine('sqlite:///tcga.db')
db.echo = False

metadata = MetaData(db)

Session = sessionmaker(db)
session = Session()

Mutations = Table('Mutations', metadata, 
	Column('file', String),
	Column('cancer', String),
	Column('Hugo_Symbol', String),
	Column('Entrez_Gene_Id', String),
	Column('Center', String),
	Column('NCBI_Build', String),
	Column('Chromosome', String),
	Column('Start_Position', Integer),
	Column('End_Position', Integer),
	Column('Variant_Classification', String), # Frame_Shift_Del, Frame_Shift_Ins, In_Frame_Del, In_Frame_Ins, Missense_Mutation, Nonsense_Mutation, Silent, Splice_Site, Translation_Start_Site, Nonstop_Mutation, 3'UTR, 3'Flank, 5'UTR, 5'Flank, IGR1 , Intron, RNA, Targeted_Region, De_novo_Start_InFrame, or De_novo_Start_OutOfFrame
	Column('Variant_Type', String), # SNP, DNP, TNP, ONP, INS, DEL, or Consolidated
	Column('Reference_Allele', String), # + strand ref allele
	Column('Tumor_Seq_Allele1', String),
	Column('Tumor_Seq_Allele2', String),
	Column('Tumor_Sample_Barcode', String),
	Column('Matched_Norm_Sample_Barcode', String),
	Column('Matched_Norm_Seq_Allele1', String),
	Column('Matched_Norm_Seq_Allele2', String),
	Column('Tumor_Validation_Allele1', String),
	Column('Tumor_Validation_Allele2', String),
	Column('Matched_Norm_Validation_Allele1', String),
	Column('Matched_Norm_Validation_Allele2', String),
	Column('Validation_Status', String),
	Column('Mutation_Status', String),
	Column('Sequencer', String),
	)


synonyms = {'Chromosome': 'Chrom'}

Snps = Table('Snps', metadata,
	Column('chromosome', String),
	Column('position', String),
	Column('gene', String),
	)

metadata.create_all(db)

def make_matrix(outputfile = 'genotype_matrix.temp', snpcountfile = 'snpcount.temp'):
	out = open(outputfile, 'w')
	snpcountout = open(snpcountfile, 'w')

	sam = session.query(Mutations).all()
	a = list(set(filter(lambda x: x.Variant_Classification != 'Silent', sam)))
	allsnps = list(set(map(lambda x: str(x.Chromosome) + ':' + str(x.Start_Position), a)))
	samples = list(set(map(lambda x: x.Tumor_Sample_Barcode, a)))

	line = 'SAMPLE,'
	snpcount = {}
	for snp in allsnps:
		line = line + snp + ','
		snpcount[snp] = 0
	out.write(line.strip(',')+'\n')

	for s in samples:
		try:
			out.write(s + ',')
			rows = filter(lambda x: x.Tumor_Sample_Barcode == s, a)
			rowsnps = set(map(lambda x: str(x.Chromosome) + ':' + str(x.Start_Position), rows))
			snps = ''
			for snp in allsnps:
				if snp in rowsnps:
					snps = snps + '1' + ','
					snpcount[snp] = snpcount[snp] + 1
				else:
					snps = snps + '0' + ','
					snpcount[snp] = snpcount[snp] + 0
			out.write(snps.strip(',')+'\n')
		except:
			continue
	json.dump(snpcount, snpcountout)

def select_top_genes(genotype_matrix_file = 'genotype_matrix.temp'):
	genotype_matrix = open(genotype_matrix_file)
	g = genotype_matrix.next()
	snps = g.strip('\n').split(',')[1:]
	alls = session.query(Snps).all()
	snptogene = {}
	for a in alls:
		snp =  a.chromosome + ':' + a.position 
		snptogene[snp] = a.gene 

	genecount = {}
	allgenes = list(set(snptogene.values()))
	for g in allgenes:
		genecount[g] = 0

	for g in genotype_matrix:
		l = g.strip('\n').split(',')[1:]
		for i in range(0, len(l)):
			genei = snptogene[snps[i]]
			genecount[genei] = genecount[genei] + int(l[i])

	sorted_genes = sorted(genecount.iteritems(), key=operator.itemgetter(1))
	json.dump(sorted_genes, open('sorted_genes', 'w'))

	return sorted_genes






def co_occur_gene(genotype_matrix_file = 'genotype_matrix.temp', genecofile = 'geneco', sorted_genes_file = 'sorted_genes'):
	'''
	take file of genotype matrix and combine each line to make a co-occurrence matrix for each snp by snp pair
	group snps by genes and calculate the gene-based co-occurrence
	'''
	genotype_matrix = open(genotype_matrix_file)
	g = genotype_matrix.next()
	snps = g.strip('\n').split(',')[1:]
	alls = session.query(Snps).all()
	snptogene = {}
	for a in alls:
		snp =  a.chromosome + ':' + a.position 
		snptogene[snp] = a.gene 

	sorted_genes = json.load(open(sorted_genes_file))

	#initialize gene by gene dictionary
	geneco = {}
	gene_genotype = {}
	#genecount = {}
	select_genes = map(lambda x: x[0], sorted_genes)[0:100]
	allgenes = list(set(snptogene.values()))[0:1000]
	allgenes.extend(select_genes)
	allgenes = list(set(allgenes))
	for g in allgenes:
		gene_genotype[g] = 0
	for g in select_genes:
		geneco[g] = {}

	#initialize co_occur count: the phastniel algorithm
	g = genotype_matrix.next()
	l = g.strip('\n').split(',')[1:]
	for i in range(0, len(l)):
		genei = snptogene[snps[i]]
		#if genei in select_genes:
		if genei in allgenes:
			if int(l[i]) == 1:
				gene_genotype[genei] = 1
		else:
			continue
		#genecount[genei] = gene_genotype[genei]
	for i in select_genes:
		for j in allgenes:
			try:
				geneco[i][j] = gene_genotype[i] * gene_genotype[j]
			except KeyError:
				continue


		#for j in range(0, len(l)):
	#		genej = snptogene[snps[j]]
	#		geneco[genei][genej] = int(l[i]) * int(l[j])



	#continue for rest of snp lines
	for g in genotype_matrix:
		print g.split(',')[0]
		l = g.strip('\n').split(',')[1:]
		for ga in allgenes:
			gene_genotype[ga] = 0
		for i in range(0, len(l)):
			genei = snptogene[snps[i]]
			#try:
			if genei in allgenes:
				if int(l[i]) == 1:
					gene_genotype[genei] = 1
				#gene_genotype[genei] = gene_genotype[genei] + int(l[i])
				#genecount[genei] = genecount[genei] + gene_genotype[genei]
		for gi in select_genes:
			for gj in allgenes:
				try:
					geneco[gi][gj] = geneco[gi][gj] + gene_genotype[gi] * gene_genotype[gj]
				except KeyError:	
					continue	
			#except KeyError:
			#	continue	


		#for i in range(0, len(l)):
		#	genei = snptogene[snps[i]]
		#	for j in range(0, len(l)):
		#		genej = snptogene[snps[j]]
		#		geneco[genei][genej] = geneco[genei][genej] + (int(l[i]) * int(l[j]))

	out = open(genecofile, 'w')
	line = 'GENES,'
	genes = geneco.keys()
	for gene in allgenes:
		line = line + gene + ','
	out.write(line.strip(',') + '\n')

	for i in select_genes:
		line =''
		for j in all_genes:
			line = i + str(geneco[i][j]) + ','
 		out.write(line.strip(',') + '\n')

 	return geneco



def count(Query):
	rows = Query.all()
	print "here"
	duplicateinds = []
	for i in range(0, len(rows)):	
		l = range(0, len(rows))
		l.remove(i)
		newrows = map(lambda x: rows[x], l)
		if rows[i] in newrows:
			duplicateinds.append(i)
	num = len(rows) - len(duplicateinds)/2
	return num


def make_snptogenetable(genotable = Mutations):
	'''
	take table from tcga maf files, and make a new table of just snp to gene relationships

	'''
	a = session.query(Mutations).all()
	snp_gene= list(set(map(lambda x: (str(x.Chromosome) + ':' + str(x.Start_Position), x.Hugo_Symbol), a)))
	snp_dict = dict(snp_gene)
	for k in snp_dict.keys():
		chrom = k.split(':')[0]
		pos = k.split(':')[1]
		inputdic = {'chromosome' : chrom, 'position' : pos, 'gene' : snp_dict[k]}
		i =Snps.insert()
		i.execute(inputdic)


#build db from maf files
def insert_from_file(filename = 'ov_liftover.aggregated.capture.tcga.uuid.somatic.maf', cancer = 'None'):
	f = open(filename)
	for h in f:
		if not h.startswith('#'):
			header = h.split('\t')
			header = map(lambda x: x.lower(), header)
			break
	print header

	for l in f:
		l = l.split('\t')
		inputdic = {}
		for c in Mutations.c:
			if c.name.lower() in header:
				inputdic.update({c.name: l[header.index(c.name.lower())]})
			elif c.name in synonyms.keys():
				inputdic.update({c.name: l[header.index(synonyms[c.name].lower())]})
		inputdic.update({'cancer' : cancer})
		i = Mutations.insert()
		i.execute(inputdic)

def convert_hg18tohg19(liftoverdir = '/home/mkagan/liftover/', chainfilename = 'hg18tohg19.over.chain'):
	Session = sessionmaker(db)
	session = Session()
	
	all36 = session.query(Mutations).filter_by(NCBI_Build = '36').all()

	chainfile = liftoverdir+chainfilename
	bed18 = open('hg18.bed', 'w')
	for a in all36:
		bed18.write('chr'+str(a.Chromosome) + '\t' + str(a.Start_Position) + '\t' + str(int(a.End_Position)+1) + '\n')
	commands.getstatusoutput("%s hg18.bed %s hg19.bed unmapped" % (liftoverdir+'liftOver', chainfile))

	unmapped = []
	lines = open('unmapped').readlines()
	for l in lines:
		if not l.startswith('#'):
			snppos = l.split('\t')[0] + ':' + l.split('\t')[1]
			unmapped.append(snppos)
	
	hg19 = open('hg19.bed')
	maf19temp = open('maf19temp', 'w')
	keys = all36[0].keys()
	for k in keys:
		maf19temp.write(k + '\t')
	maf19temp.write('\n')	

	for a in all36:
		snppos = 'chr' + str(a.Chromosome) + ':' + str(a.Start_Position)
		if snppos not in unmapped:
			try:
				l = hg19.next()
			except StopIteration:
				break
			newchrom = l.split('\t')[0].split('chr')[1]
			newstart = l.split('\t')[1]
			newend = str(int(newstart) + 1)
			a.Chromosome = newchrom
			a.Start_Position = newstart
			a.End_Position = newend
			newline = ''
			for k in keys:
				newline = newline + str(getattr(a, k)) + '\t'
			maf19temp.write(newline.strip('\t') + '\n')

			#newa = session.query(mutations).filter_by(Tumor_Sample_Barcode = a.Tumor_Sample_Barcode, Start_Position = a.Start_Position)
			#newa.update({"Chromosome": newchrom, "Start_Position": newstart, "End_Position": newend, "NCBI_Build" : 37}, synchronize_session=False)
		else:
			print snppos
			continue
	maf19temp.close()
	
	session.query(Mutations).filter_by(NCBI_Build = '36').delete(synchronize_session=False)
	session.commit()
	maf19 = open('maf19temp').readlines()
	for line in maf19[1:]:
		inputdic = {}
		for i,k in enumerate(Mutations.c):
			l = line.split('\t')
			try:
				inputdic.update({k.name: l[i]})
				inputdic.update({'NCBI_Build': '37'})
				i = Mutations.insert()
				i.execute(inputdic)
			except IndexError:
				continue
 
def co_occurr(genotype_matrix_file = 'genotype_matrix.temp'):
	'''
	take file of genotype matrix and combine each line to make a co-occurrence matrix for each snp by snp pair
	'''
	genotype_matrix = open(genotype_matrix_file)
	g = genotype_matrix.next()
	snps = g.strip('\n').split(',')[1:]
	snpco = {}
	for s in snps:
		snpco[s] = {}

	#initialize k=1 iteration	
	g = genotype_matrix.next()
	l = g.strip('\n').split(',')[1:]
	for i in range(0, len(l)):
		for j in range(0, len(l)):
			snpco[snps[i]][snps[j]] = int(l[i]) * int(l[j])

	#continue for rest of the snp lines
	for g in genotype_matrix:
		l = g.strip('\n').split(',')[1:]
		for i in range(0, len(l)):
			for j in range(0, len(l)):
				snpco[snps[i]][snps[j]] = snpco[snps[i]][snps[j]] + (int(l[i]) * int(l[j]))

 	out = open('snpco', 'w')
 	line = 'SNPs,'
 	for i in snps:
 		line = line + i + ','
 	out.write(line.strip(',') + '\n')

 	for i in snps:
 		out.write(i + ',')
 		line = ''
 		for j in snps:
 			line = line + str(snpco[i][j]) + ','
 		out.write(line.strip(',') + '\n')

	return snpco




if __name__ == "__main__":
	
	args = sys.argv[1:]
	print args 
	if '-insert' in args:
		f = args[1]
		cancertype = args[2]
		insert_from_file(filename = f, cancer = cancertype)
	if '-make_matrix' in args:
		make_matrix(args[1], args[2])
	if '-co_occur_gene' in args:
		co_occur_gene(args[1], args[2])

