from sqlalchemy import *
import sys
from sqlalchemy.orm import sessionmaker
import commands
import os

db = create_engine('sqlite:///tcga.db')
db.echo = True

metadata = MetaData(db)

Session = sessionmaker(db)
session = Session()

mutations = Table('mutations', metadata, 
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

metadata.create_all(db)

synonyms = {'Chromosome': 'Chrom'}

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
		for c in mutations.c:
			if c.name.lower() in header:
				inputdic.update({c.name: l[header.index(c.name.lower())]})
			elif c.name in synonyms.keys():
				inputdic.update({c.name: l[header.index(synonyms[c.name].lower())]})
		inputdic.update({'cancer' : cancer})
		i = mutations.insert()
		i.execute(inputdic)

def convert_hg18tohg19(liftoverdir = '/home/mkagan/liftover/', chainfilename = 'hg18tohg19.over.chain'):
	Session = sessionmaker(db)
	session = Session()
	'''
	all36 = session.query(mutations).filter_by(NCBI_Build = '36').all()

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
	'''
	#session.query(mutations).filter_by(NCBI_Build = '36').delete(synchronize_session=False)
	maf19 = open('maf19temp').readlines()
	for line in maf19[1:]:
		inputdic = {}
		for i,k in enumerate(mutations.c):
			l = line.split('\t')
			inputdic.update({k.name: l[i]})
		inputdic.update({'NCBI_Build': '37'})
		i = mutations.insert()
		i.execute(inputdic)







def bed_to_mafstyle():
	'''
	convert bed file to something like a maf file, basically write the db to a file, then read it back
	'''




if __name__ == "__main__":
	
	args = sys.argv[1:]
	print args 
	if '-insert' in args:
		f = args[1]
		cancertype = args[2]
		insert_from_file(filename = f, cancer = cancertype)

