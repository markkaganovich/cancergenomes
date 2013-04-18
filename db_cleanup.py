# clean up DB functions
from sqlalchemy import *
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import and_, or_
import sys
import json
import commands
import os
import operator


def convert_hg18tohg19(db, tablename, build_col = 'ncbi_build', liftoverdir = '/home/mkagan/liftover/', chainfilename = 'hg18tohg19.over.chain'):
    db.dialect.get_table_names(db.connect())

    Session = sessionmaker(db)
    session = Session()
    metadata = MetaData(db)
    table = Table(tablename, metadata, autoload = True)
   
    all36 = session.query(table).filter(getattr(table.c, build_col) == '36').all()

    chainfile = liftoverdir+chainfilename
    bed18 = open('hg18.bed', 'w')
    for a in all36:
        bed18.write('chr'+str(a.chrom) + '\t' + str(a.start_position) + '\t' + str(int(a.end_position)+1) + '\n')
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
        snppos = 'chr' + str(a.chrom) + ':' + str(a.start_position)
        if snppos not in unmapped:
            try:
                l = hg19.next()
            except StopIteration:
                break
            newchrom = l.split('\t')[0].split('chr')[1]
            newstart = l.split('\t')[1]
            newend = str(int(newstart) + 1)
            a.chrom = newchrom
            a.start_position = newstart
            a.end_position = newend
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

    session.query(table).filter(getattr(table.c, build_col) == '36').delete(synchronize_session=False)
    session.commit()
    
    maf19 = open('maf19temp').readlines()
    for line in maf19[1:]:
        inputdic = {}
        for i,k in enumerate(table.c):
            l = line.split('\t')
            try:
                inputdic.update({k.name: l[i]})
                inputdic.update({build_col: '37'})
                i = table.insert()
                i.execute(inputdic)
            except IndexError:
                continue


build_col = 'ncbi_build'
tablename = 'mutations_v1'
db = create_engine('sqlite:///tcga_somatic.db', echo = False)

convert_hg18tohg19(db, tablename, build_col = build_col)
