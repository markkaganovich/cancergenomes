import operator

execfile('init.py')
from rebuildanal_helper import *

if 'snp2aa.json' in os.listdir('./'):
	snp2aa = json.load(open('snp2aa.json'))
else:	
	snp2aa = rebuildanal_helper.get_snp_to_aa(tcga_rows)



bps = map(lambda x: x.chrom + ':' + str(x.start_position), tcga_rows)

#############
# do it just for individual genes
#
# 
gene = 'BRAF'
braf = gene_bps('BRAF', tcga_rows)

'''
oncogenes = json.load(open('oncogenes'))
for g in oncogenes:
	bps = gene_bps(g, tcga_rows)
	try:
		file = open('./gene_bps/' + g + '_bps.json', 'w')
		json.dump(dict(Counter(bps)), file)
		file.close()
	except IOError:
		continue


for g in oncogenes:
	file = g + '_bps.json'
	try:
		scp('/genes_bps/'+file)
		plot_file(file, './graphs/')
	except:
		continue
	
'''

def select(key, value, dataset):
	return filter(lambda x : getattr(x, key) == value, dataset)


def gene_bps(gene, dataset):
	gene_rows = select('hugo_symbol', gene, dataset)
	bps = sorted(map(lambda x: x.start_position, gene_rows))
	return bps

def scp(remotefile, remotehost = "mkagan@scg3", localfile = None):
	print remotefile, localfile, remotehost
	remotefile_path = "~/mark/cancergenomes/"
	if not localfile:
		localfile = remotefile
	os.system('scp  "%s:%s" "%s"' % (remotehost, remotefile_path + remotefile, localfile))

def plot_file_old(json_file):
	file = open(json_file)
	data = json.load(file)
	sorted_data = sorted(data.iteritems(), key = operator.itemgetter(0))
	genome_x = map(lambda x: int(x[0]), sorted_data)
	gene_x = map(lambda x: x - min(genome_x), genome_x)
	y = map(lambda x: x[1], sorted_data)
	plt.figure(1)
	plt.axes(frameon = False)
	plt.plot(gene_x, y, 'bo', alpha=.5)
	plt.plot(gene_x, y, 'k', alpha=.2)
	plt.grid(True)
	plt.show()


def plot_file(json_file):
	fig1 = plt.figure(facecolor=".75")
	ax1 = plt.axes(frameon=False)
	ax1.get_xaxis().tick_bottom()   # Turn off ticks at top of plot

	ax1.axes.get_xaxis().set_visible(False)
	#ax1.axes.get_yaxis().set_visible(False)     # Hide y axis

	# Add a plot
	file = open(json_file)
	data = json.load(file)
	sorted_data = sorted(data.iteritems(), key = operator.itemgetter(0))
	genome_x = map(lambda x: int(x[0]), sorted_data)
	gene_x = map(lambda x: x - min(genome_x), genome_x)
	y = map(lambda x: x[1], sorted_data)
	ax1.plot(gene_x, y, 'bo', alpha=.4)

	#plt.legend(loc='lower right')

	# Draw the x axis line
	# Note that this must be done after plotting, to get the correct
	# view interval
	xmin, xmax = ax1.get_xaxis().get_view_interval()
	ymin, ymax = ax1.get_yaxis().get_view_interval()
	ax1.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))

	plt.show()



###########
# compare bps
#
#
file = open('aa2snp.json')
aa2snp = json.load(file)
file.close()

bps = map(lambda x: x.hugo_symbol + ':' +x.chrom + ':' + str(x.start_position), tcga_rows)
counts_bp = json.load(open('counts_bp'))
sim_gene_results = sim_output('simulations_bp_total.log')
pile_ups = pile_up_bp(sim_gene_results, bps, counts_bp)
pileupsorted = sorted(pile_ups.iteritems(), key=operator.itemgetter(1), reverse=True)

aa2snp = json.load(open('aa2snp.json'))
snp2aa = json.load(open('snp2aa.json'))

# get top pileups for bps. scatterplot the top 1000
bp_pairs = {}
ba = []
peaks = []
for p in pileupsorted[0:1000]:
	aa = snp2aa[p[0]]
	both_snps = aa2snp[aa]
	if len(both_snps) > 2:
		other_bp = filter(lambda x: x != p[0], both_snps)[0]
		other_peak = filter(lambda x: x[0] == other_bp, pileupsorted)[0][1]
		bp_pairs[aa] = (p[1], other_peak)
		ba.append(aa)
		peaks.append((p[1], other_peak))
	else:
		continue

tcga_rows_nonsyn = filter(lambda x: x.variant_classification != 'Silent', tcga_rows)
residues = map(lambda s: s.residue, tcga_rows_nonsyn)
#counts_aa = json.load(open('counts_aa'))
sim_gene_results = sim_output('simulations_aa.log')
#pile_ups = pile_up(sim_gene_results, residues, counts_aa)
#pileupsorted = sorted(pile_ups.iteritems(), key=operator.itemgetter(1), reverse=True)

counts_nonsyn  = json.load(open('counts_non_syn.json'))
pileups = pile_up(sim_gene_results, residues, counts_nonsyn)


## nonsyn (7:37)
tcga_rows_nonsyn = filter(lambda x: x.variant_classification != 'Silent', tcga_rows)
sim_gene_results = sim_output('simulations_aa_nonsyn.log')
counts_nonsyn  = json.load(open('counts_non_syn.json'))
residues = map(lambda s: s.residue, tcga_rows_nonsyn)
pileups = pile_up(sim_gene_results, residues, counts_nonsyn)

#silent
tcga_rows_silent = filter(lambda x: x.variant_classification == 'Silent', tcga_rows)
sim_gene_results = sim_output('simulations_aa_silent.log')
counts_silent  = json.load(open('counts_silent.json'))
residues = map(lambda s: s.residue, tcga_rows_silent)
pileups_silent = pile_up(sim_gene_results, residues, counts_silent)





def plot_sim(gene, rows, pileups, counts):
	braf_residues = map(lambda y: y.residue, filter(lambda x: x.hugo_symbol == gene, rows))  #tcga_rows_nonsyn
	#t = map(lambda x:  pileups[x], braf_residues)  #pile_ups_test
	r = map(lambda x: int(x.split(':')[1]), counts[gene].keys())  #counts_nonsyn
	x_range = list(range(min(r), max(r)))

	v = []
	for i in x_range:
		if gene+':'+str(i) not in counts[gene].keys():
			v.append(0)
		else:
			v.append(counts[gene][gene +':'+str(i)])

	return x_range, v, t

trace = Scatter(x = x_range, y = v)
data = Data([trace])




fig1 = plt.figure(facecolor=".65")

ax1 = plt.axes(frameon=False)
ax1.axes.get_xaxis().set_visible(False)
ax1.tick_params(axis='y', direction='out', color='#ffffff')
ax1.bar(np.arange(19), bin_values2, color='#ffffff', edgecolor="None")
plt.savefig('test.png', transparent="True")




x_nonsyn = list(range(c_sorted[0][0], c_sorted[-1][0]))
y_nonsyn = []
for i in range(0,len(x_nonsyn)):
    if x_nonsyn[i] not in c.keys():
        y_nonsyn.append(0)
    else:
        y_nonsyn.append(c[x_nonsyn[i]])

tcga_rows_nonsyn = filter(lambda x: x.variant_classification != 'Silent', tcga_rows)






