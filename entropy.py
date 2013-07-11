#calculate entropy

from sqlalchemy.orm import backref, mapper, relation, sessionmaker
from sqlalchemy.sql import and_, or_
from sqlalchemy import create_engine, Column, String, Integer, MetaData, Table
import sys
import json
import commands
import os
import csv 
import headers
import operator
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
import math
import scipy

db = create_engine('sqlite:///tcga_somatic.db', echo = False)

first_pass_mutations = ['Frame_Shift_Del', 'Frame_Shift_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Translation_Start_Site']

metadata = MetaData(db)
Session = sessionmaker(db)
session = Session()

Mutations = Table('mutations_v1', metadata, autoload = True)
m = Mutations.select().execute()
rows = filter(lambda x: x.variant_classification in first_pass_mutations, m)

def make_gene_snp(rows):
    gene_snp = {}
    for r in rows:
        gene = r.hugo_symbol
        snp = r.chrom + ':' + r.start_position
        try:
            gene_snp[gene].append(snp)
        except KeyError:
            gene_snp[gene] = [snp]
    return gene_snp

gene_snp = make_gene_snp(rows)
counts = {}
for g in gene_snp.keys():
    counts[g] = Counter(gene_snp[g])

def poisson(k, l):
    try:
        p = math.pow(l, k)/math.factorial(k) * np.exp(-1*l)
    except OverflowError:
        p = 0
    return p


# caculate entropty
gene_sequences = json.load(open('gene_sequences'))
genes = gene_snp.keys()
entropy_possion = {}
for g in genes:
    v = counts[g].values()
    #v = map(lambda x: long(x), v)
    try:
        if gene_sequences[g].__len__() > 3:
        #entropy[g] = -1 * sum(map(lambda x: x * float(len(v))/len(gene_sequences[g]) * np.log(x* (float(len(v))/len(gene_sequences[g]))), v))
        #entropy[g] = -1 * float(sum(map(lambda x: np.log(x), v))) / (len(gene_sequences[g])*2) 
            total = []
            for k in v:
                l = float(sum(v))/ len(v)
                p = poisson(k, l)
                total.append(p)
            entropy_possion[g] =  np.product(total)
    except KeyError:
        continue

a1 = sorted(entropy_possion.iteritems(), key = operator.itemgetter(1))
ents1 = map(lambda x: x[1], a1)
d1 = map(lambda x: x[0], a1)
'''
entropy_normalized = {}
for g in genes:
    if sum(counts[g].values()) > 30:
        entropy_normalized[g] = entropy[g] / sum(counts[g].values())

'''

def combination(n,c):
    return np.product(range(n-c+1, n+1))/float(math.factorial(c))

def binomial(n, k, p):
    return combination(n,k) * math.pow(p,k) * math.pow(1-p, n-k)

'''
comb = {}
for g in genes:
    v = counts[g].values()
    try: 
        l = gene_sequences[g].__len__()*2
        if l > 6:
            p = len(v) 
            print l
            print p
            den = combination(l,p)
            print den
            total = map(combination, [p]*len(v), v)
            print total
            comb[g] = np.product(total)/math.pow(den, len(v))
    except KeyError:
        continue
'''


if 'counts_aa' in os.listdir('./'):
    counts_aa = json.load(open('counts_aa'))
else:
    make_counts_aa(genes, counts)

prtn_len = json.load(open('prtn_len'))

hack = {}
hack_max = {}
keys = set(gene_sequences.keys())
prtn_len_keys = set(prtn_len.keys())
for g in genes:
    l = counts[g].values()    
    v = np.array(counts_aa[g].values())
    if len(v) > 0 and g in prtn_len_keys:
        hack_max[g] = max(v/np.mean(v))   
    

poisson_residues = {}
for g in counts_aa:
    m = sum(counts_aa[g].values())
    try:
        p = float(m) / prtn_len[g]
    except KeyError:
        continue
    for aa in counts_aa[g]:
        if counts_aa[g][aa] > 1:
            poisson_residues[g+':'+aa] = poisson(counts_aa[g][aa], p)

'''
poisson_by_gene = {}
for k in poisson_residues:
    g = k.split(':')[0]
    if g in poisson_by_gene.keys():
        poisson_by_gene[g].append(poisson_residues[k])
    else:
        poisson_by_gene[g] = [poisson_residues[k]]
'''


for g in counts_aa:     
    try:        
        gi[g] = genes_in_poisson.index(g)
    except ValueError:
        continue


a = sorted(hack_max.iteritems(), key = operator.itemgetter(1), reverse=True)
ents = map(lambda x: x[1], a)
d = map(lambda x: x[0], a)
mut_rate = map(lambda x: sum(counts_aa[x].values())/float(prtn_len[x]), d)
prod = map(lambda x: ents[x] * mut_rate[x], range(0, len(d)))
div = map(lambda x: ents[x] / mut_rate[x], range(0, len(d)))
    

#plt.scatter(ents[0:30], mut_rate[0:30]) 
f = [d[i] for i in range(0,len(d)) if mut_rate[i] >= .2 and ents[i]<2]




'''
fig = plt.figure()
ax = fig.add_subplot(111)
ind = list(range(0, len(ents)))
width = .2
print "making bar"
bar = ax.bar(ind, ents, width, color="r")
print "saving ..."
plt.savefig('hack.png')
'''
#if g in keys:
#hack[g] = np.dot(v,v)/float(sum(v))/len(gene_sequences[g])


m_mut = sorted(counts.iteritems(), key = lambda x: sum(x[1].values()), reverse=True)
'''r_mut = map(lambda x: x[0], m_mut)

distances = {}
for g in genes:
    distances[g] = r_mut.index(g) - d.index(g) 
e = sorted(distances.iteritems(), key = operator.itemgetter(1), reverse=True)
'''

def entr(v):
    return sum(map(lambda x: (float(x)/1000) * np.log(float(x)/1000), v))



#convert coutns to counts_aa by hitting bp_to_aa hash
bp_to_aa = json.load(open('bp_to_aa'))
def make_counts_aa(genes, counts):
    silent = []
    counts_aa = {}
    for g in genes:
        counts_aa[g] = {}
        for c in counts[g].keys():
            try:
                aa = bp_to_aa[c]['pos']
                if aa in counts_aa[g].keys():
                    counts_aa[g][aa] = counts_aa[g][aa] + counts[g][c]
                else:
                    counts_aa[g][aa] = counts[g][c]
            except KeyError:
                silent.append(c)
    return counts_aa

'''

binomial_genes = {}
for g in genes:
    v = counts_aa[g].values()
    n = sum(v)
    try:
        p = float(n)/prtn_len[g]
    except KeyError:
        continue
    binomial_genes[g] = map(binomial, [n]*len(v), v, [p]*len(v))
'''

s = set(poisson_residues.keys())
peak_rows = filter(lambda x: x.hugo_symbol +':'+str(x.residue['pos']) in s, rows)
types = list(set(map(lambda x: x.cancer_type, rows)))
a = sorted(poisson_residues.iteritems(), key = operator.itemgetter(1))
residues = map(lambda x: x[0], a)
residue_by_type = {}

singletons = filter(lambda x: x.hugo_symbol +':'+str(x.residue['pos']) not in s, rows)


#initialize
for r in residues:
    residue_by_type[r] = {}
    for t in types:
        residue_by_type[r][t] = 0

for p in peak_rows:
    r = p.hugo_symbol +':' + str(p.residue['pos'])
    print r
    if r in residues:
        t = p.cancer_type
        residue_by_type[r][t] = residue_by_type[r][t] + 1


res2 = filter(lambda x: sum(residue_by_type[x].values()) == 2, residue_by_type.keys())
res3 = filter(lambda x: sum(residue_by_type[x].values()) == 3, residue_by_type.keys())
res4 = filter(lambda x: sum(residue_by_type[x].values()) == 4, residue_by_type.keys())
res5 = filter(lambda x: sum(residue_by_type[x].values()) == 5, residue_by_type.keys())
res6 = filter(lambda x: sum(residue_by_type[x].values()) == 6, residue_by_type.keys())
res7 = filter(lambda x: sum(residue_by_type[x].values()) == 7, residue_by_type.keys())
res8plus = filter(lambda x: sum(residue_by_type[x].values()) > 7, residue_by_type.keys()) 

res = res8plus
res_type = {}
for r in res:
    total = sum(residue_by_type[r].values())
    res_type[r] = {}
    for t in types:
        res_type[r][t] = float(residue_by_type[r][t])/ total


single_types = {}
for t in types:
    single_types[t] = 0

for s in singletons:
    single_types[s.cancer_type] = single_types[s.cancer_type] + 1

#print to file
output = open('residue_by_type_frac','w')
output = open('res_type_8plus', 'w')
c = '\t'
for t in types:
    c = c+t+','

output.write(c.strip(',')+'\n')
for r in res:
    l = str(r) + ','
    for t in types:
        l = l+str(res_type[r][t]) + ','
    l = l.strip(',')
    l = l+'\n'
    output.write(l)

total = sum(single_types.values())
single_types_frac = {}
for t in single_types:
    single_types_frac[t] = single_types[t]/float(total)

cancer_distr = {}
for t in types:
    cancer_distr[t] = 0

for r in rows:
    cancer_distr[r.cancer_type] = cancer_distr[r.cancer_type] + 1

total = sum(cancer_distr.values())
for t in types:
    cancer_distr[t] = float(cancer_distr[t])/total


res_type_norm = {}
for r in res:
    res_type_norm[r] = {}
    for t in types:
        res_type_norm[r][t] = res_type[r][t]/single_types_frac[t]
    

#print to file
output = open('res_type_8plus_norm', 'w')
c = '\t'
for t in types:
    c = c+t+','

output.write(c.strip(',')+'\n')
for r in res:
    l = str(r) + ','
    for t in types:
        l = l+str(res_type_norm[r][t]) + ','
    l = l.strip(',')
    l = l+'\n'
    output.write(l)

def table(data):
    output = {}
    print type(data)
    if type(data) is list:
        for t in data:
            output[t] = 0
        for t in data:
            output[t] = output[t] + 1 
    
    if type(data) == dict:     
        for t in data.values():
            output[t] = 0
        for k in data.keys():
            output[data[k]] = output[data[k]] + 1
    
    return output





tumor_suppr_rankings = []

for t in tumor_suppr:
    try:
        tumor_suppr_rankings.append(g.index(t))
    except:
        continue


tumor_suppr_mean_rankings = []

g_mean = map(lambda x: x[0], min_res_by_sample_sorted)

for t in tumor_suppr:
    try:
        tumor_suppr_mean_rankings.append(g_mean.index(t))
    except:
        continue

oncogene_mean_rankings = []
for o in oncogenes:
    try:
        oncogene_mean_rankings.append(g_mean.index(o))
    except:
        continue


onco_min = map(lambda x: x[1], filter(lambda x: x[0] in oncogenes, min_res_sorted))
tumor_min = map(lambda x: x[1], filter(lambda x: x[0] in tumor_suppr, min_res_sorted))

onco_mean = map(lambda x: x[1], filter(lambda x: x[0] in oncogenes, mean_res_sorted))
tumor_mean = map(lambda x: x[1], filter(lambda x: x[0] in tumor_suppr, mean_res_sorted))

cancer_genes = []
cancer_genes.extend(oncogenes)
cancer_genes.extend(tumor_suppr)

cancer_gene_min = filter(lambda x: x[0] in cancer_genes, min_res_sorted)
cancer_gene_mean = filter(lambda x: x[0] in cancer_genes, mean_res_sorted)

genes = poisson_by_gene.keys()

sim_min = []
sim_mean = []
sim_inds = []
for j in range(0, 1000):
    random_genes = map(lambda x: genes[x], [np.random.randint(len(genes)) for i in range(0, 105)])
    #random_min = filter(lambda x: x[0] in random_genes, min_res_sorted)
    #r_min = np.mean(map(lambda x: x[1], random_min))
    #random_mean = filter(lambda x: x[0] in random_genes, mean_res_sorted)
    #r_mean = np.mean(map(lambda x: x[1], random_mean))
    #sim_min.append(r_min)
    #sim_mean.append(r_mean)
    inds_rand = map(lambda x: g.index(x), filter(lambda y: y in g, random_genes))
    sim_inds.append(inds_rand)

np.mean(sim_min)
np.mean(sim_mean)


s = sum(map(lambda x: np.array(sorted(x)), sim_inds))


p_s = sorted(poisson_residues.iteritems(), key = operator.itemgetter(1))

ps = map(lambda x: x[0], p_s)
setps = set(ps[0:1000]) 

peak_rows = filter(lambda x: x.hugo_symbol +':'+str(x.residue['pos']) in setps, rows)




types = list(set(map(lambda x: x.cancer_type, rows)))
#a = sorted(poisson_residues.iteritems(), key = operator.itemgetter(1))
#residues = map(lambda x: x[0], a)
residue_by_type = {}

singletons = filter(lambda x: x.hugo_symbol +':'+str(x.residue['pos']) not in setps, rows)


#initialize
res = list(setps)
for r in res:
    residue_by_type[r] = {}
    for t in types:
        residue_by_type[r][t] = 0

for p in peak_rows:
    r = p.hugo_symbol +':' + str(p.residue['pos'])
    print r
    if r in residues:
        t = p.cancer_type
        residue_by_type[r][t] = residue_by_type[r][t] + 1

#general distribution of cancers
cancer_distr = {}
for t in types:
    cancer_distr[t] = 0

for r in rows:
    cancer_distr[r.cancer_type] = cancer_distr[r.cancer_type] + 1

total = sum(cancer_distr.values())
for t in types:
    cancer_distr[t] = float(cancer_distr[t])/total


res_type_norm = {}
for r in res:
    res_type_norm[r] = {}
    for t in types:
        res_type_norm[r][t] = residue_by_type[r][t]/cancer_distr[t]


#print to file
output = open('residue_by_type_frac_ps_norm','w')

c = '\t'
for t in types:
    c = c+t+','

output.write(c.strip(',')+'\n')
for r in res:
    l = str(r) + ','
    for t in types:
        l = l+str(res_type_norm[r][t]) + ','
    l = l.strip(',')
    l = l+'\n'
    output.write(l)

single_genes_in_res = map(lambda x: x.split(':')[0], res) 
singletongenes = map(lambda x: x.hugo_symbol, singletons)
single_genes_in_res = filter(lambda x: x in singletongenes, single_genes_in_res)
single_res_rows = filter(lambda x: x.hugo_symbol in single_genes_in_res, singletons)
#single_res = map(lambda x: x.hugo_symbol + ':' + str(x.residue['pos']), single_res_rows)

res.remove('WASH1:134')   #pseudogene, and it has a peak but no singletons! which is crazy

single_by_type = {}
for r in single_genes_in_res:
    single_by_type[r] = {}
    for t in types:
        single_by_type[r][t] = 0.00001

for p in single_res_rows:
    r = p.hugo_symbol
    print r
    t = p.cancer_type
    single_by_type[r][t] = single_by_type[r][t] + 1

for r in single_by_type:
    total = sum(single_by_type[r].values())
    for t in types:
        single_by_type[r][t] = single_by_type[r][t]/float(total)

res_type_norm_gene = {}
for r in res:
    res_type_norm_gene[r] = {}
    gene = r.split(':')[0]
    for t in types:
        res_type_norm_gene[r][t] = residue_by_type[r][t]/single_by_type[gene][t]


output = open('residue_by_type_norm_gene','w')

c = '\t'
for t in types:
    c = c+t+','

output.write(c.strip(',')+'\n')
for r in res:
    l = str(r) + ','
    for t in types:
        l = l+str(res_type_norm_gene[r][t]) + ','
    l = l.strip(',')
    l = l+'\n'
    output.write(l)

output.close()



'''
from anal.py
'''
gm = sorted(genecount.iteritems(), key=operator.itemgetter(1), reverse=True)
mutation_genes = map(lambda x: x[0], gm) 


>>> np.mean(map(lambda x: mutation_genes.index(x), oncogenes))
7646.3090909090906
>>> np.mean(map(lambda x: oncos.index(x), oncogenes))         
7764.8000000000002
>>> np.mean(map(lambda x: mutation_genes.index(x), tumor_suppr))
4558.3442622950815
>>> np.mean(map(lambda x: oncos.index(x), tumor_suppr))         
9250.0983606557384






total = sum(single_types.values())
single_types_frac = {}
for t in single_types:
    single_types_frac[t] = single_types[t]/float(total)





    

