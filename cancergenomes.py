from sqlalchemy import *
import sys

db = create_engine('sqlite:///tcga.db')
db.echo = True

metadata = MetaData(db)

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


if __name__ == "__main__":
	
	args = sys.argv[1:]
	print args 
	if '-insert' in args:
		f = args[1]
		cancertype = args[2]
		insert_from_file(filename = f, cancer = cancertype)

