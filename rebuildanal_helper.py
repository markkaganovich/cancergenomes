import os
import matplotlib.pyplot as plt
import simplejson as json
import operator
from collections import Counter
import plotly.plotly as py
from plotly.graph_objs import *
py.sign_in("mark.kaganovich", "62d1k8xa77")
import numpy as np


def sim_output(simulation_file = 'simulations.log'):

	sim_gene_results = {}
	lines = open(simulation_file).readlines()
	simgenes = map(lambda x: x.split('\t')[1].strip('  '), lines)
	for gene in simgenes:
		sim_gene_results[gene] = {}

	for l in lines:
		tabs = l.split('\t')
		g = tabs[1].strip('  ')
		sim_gene_results[g]= {}
		sim_gene_results[g]['mean'] = float(tabs[2])
		sim_gene_results[g]['std'] = float(tabs[3])
		sim_gene_results[g]['metric'] = float(tabs[4].strip('\n'))

	return sim_gene_results

## per residue

def pile_up(sim_gene_results, residues, counts_aa):
	pile_ups = {}
	what = []
	for res in residues:
		gene = res.split(':')[0]
		try:
			mean = sim_gene_results[gene]['mean']
			std = sim_gene_results[gene]['std']
		except KeyError:
			continue
		try:
			if res not in counts_aa[gene].keys():
				what.append(res)
			else:
				val = counts_aa[gene][res]
				pile_ups[res] = (val-mean)/std
		except:
			continue

	return pile_ups



## HACK. pile_up exception for BP case (the keys of counts_bp are not gene:position, but rather just position)
def pile_up_bp(sim_gene_results, residues, counts_aa):
	pile_ups = {}
	what = []
	for res in residues:
		gene = res.split(':')[0]
		chrom = res.split(':')[1]
		start_position = str(res.split(':')[2])
		chrom_pos = chrom + ':' + start_position
		try:
			mean = sim_gene_results[gene]['mean']
			std = sim_gene_results[gene]['std']
		except KeyError:
			continue
		try:
			if start_position not in counts_aa[gene].keys():
				what.append(chrom_pos)
			else:
				val = counts_aa[gene][start_position]
				pile_ups[chrom_pos] = (val-mean)/std
		except:
			continue

	return pile_ups




def get_table(chrom, pos):
	pos = int(pos)
	chrom = chrom.strip(' ')
	db_name = 'SIFT/Human_CHR' + chrom + '.sqlite'
	conn = sqlite.connect(db_name)
	cursor = conn.cursor()
	tables = map(lambda x: x[0], cursor.execute("select name from sqlite_master where type = 'table'").fetchall())
	if tables is None:
		return None
	table_ranges = map(lambda x: (int(x.split('_')[1]), int(x.split('_')[2])), tables)
	range_list = filter(lambda x: pos >= x[0] and pos <= x[1], table_ranges)
	if range_list == []:
		return None
	select_range = range_list[0]
	select_table = 'chr' + chrom + '_' + str(select_range[0])+ '_' + str(select_range[1])
	table = cursor.execute("select COORD2, AAPOS2 from " + select_table + " where SNP = 'Reference'").fetchall()
	table_dic = {}
	for t in table:
		table_dic[t[0]] = t[1]
	positions = set(table_dic.keys())
	return select_range, table_dic, positions


def get_snp_to_aa(rows):
	snp_to_aa = {}
	r = rows[1]
	curr_chrom = str(r.chrom).strip(' ')
	pos = int(r.start_position) 
	current_range, current_table, positions = get_table(curr_chrom, pos)
	for c in map(lambda x: str(x), range(1,22)) + ['X', 'Y']:
		select_rows = filter(lambda x: x.chrom == c, rows)
		sorted_rows = sorted(select_rows, key = lambda x: int(x.start_position))
		print c
		for r in sorted_rows:
			chrom = str(r.chrom).strip(' ')
			pos = int(r.start_position)
			if not (pos > current_range[0] and pos < current_range[1]):
				try:
					[current_range, current_table, positions] = get_table(chrom, pos)
					print "changing tables"
				except TypeError:
					continue
			if pos in positions:
				aa = current_table[pos]
				snp_to_aa[chrom +':' + str(pos)] = r.hugo_symbol + ':' + str(aa)
			print chrom + ':' + str(pos)
	return snp_to_aa


def make_counts_aa(rows):
	counts_aa = {}
	for r in rows:
		if r.hugo_symbol not in counts_aa.keys():
			counts_aa[r.hugo_symbol] = {}
		try:
			if r.residue in counts_aa[r.hugo_symbol].keys():
				counts_aa[r.hugo_symbol][r.residue] += 1
			else:
				counts_aa[r.hugo_symbol][r.residue] = 1
		except AttributeError:
			continue
	return counts_aa


def query_position(chrom_pos, rows):
	input_chrom = chrom_pos.split(':')[0]
	input_pos = chrom_pos.split(':')[1]
	return filter(lambda x: x.chrom == input_chrom and str(x.start_position) == input_pos, rows)


def select(key, value, dataset):
	return filter(lambda x : getattr(x, key) == value, dataset)


def gene_bps(gene, dataset):
	gene_rows = select('hugo_symbol', gene, dataset)
	bps = map(lambda x: x.start_position, gene_rows)
	return bps

def scp(remotefile, remotehost = "mkagan@scg3", localfile = None):
	print remotefile, localfile, remotehost
	remotefile_path = "~/mark/cancergenomes/"
	if not localfile:
		localfile = remotefile
	os.system('scp  "%s:%s" "%s"' % (remotehost, remotefile_path + remotefile, localfile))


def plot_file(json_file, show_or_save = 'save', input_path = './', output_path = './'):
	fig1 = plt.figure(facecolor=".75")
	ax1 = plt.axes(frameon=False)
	ax1.get_xaxis().tick_bottom()   # Turn off ticks at top of plot

	#ax1.axes.get_xaxis().set_visible(False)
	#ax1.axes.get_yaxis().set_visible(False)     # Hide y axis

	# Add a plot
	file = open(input_path + json_file)
	data = json.load(file)
	sorted_data = sorted(data.iteritems(), key = operator.itemgetter(0))
	genome_x = map(lambda x: int(x[0]), sorted_data)
	gene_x = map(lambda x: x - min(genome_x), genome_x)
	y = map(lambda x: x[1], sorted_data)
	ax1.plot(gene_x, y, 'bo', alpha=.7, clip_on = False)

	#plt.legend(loc='lower right')

	# Draw the x axis line
	# Note that this must be done after plotting, to get the correct
	# view interval
	xmin, xmax = ax1.get_xaxis().get_view_interval()
	ymin, ymax = ax1.get_yaxis().get_view_interval()
	ax1.add_artist(plt.Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=2))

	if show_or_save == 'save':
		plt.savefig(output_path + json_file + '.png', transparent = True)
	if show_or_save == 'show':
		plt.show()

def sim_output(simulation_file = 'simulations.log'):

	sim_gene_results = {}
	lines = open(simulation_file).readlines()
	simgenes = map(lambda x: x.split('\t')[1].strip('  '), lines)
	for gene in simgenes:
		sim_gene_results[gene] = {}

	for l in lines:
		tabs = l.split('\t')
		g = tabs[1].strip('  ')
		sim_gene_results[g]= {}
		sim_gene_results[g]['mean'] = float(tabs[2])
		sim_gene_results[g]['std'] = float(tabs[3])
		sim_gene_results[g]['metric'] = float(tabs[4].strip('\n'))

	return sim_gene_results

def pile_up(sim_gene_results, residues, counts_aa):
	pile_ups = {}
	what = []
	for res in residues:
		gene = res.split(':')[0]
		try:
			mean = sim_gene_results[gene]['mean']
			std = sim_gene_results[gene]['std']
			n = len(counts_aa[gene])
		except KeyError:
			continue
		try:
			if res not in counts_aa[gene].keys():
				what.append(res)
			else:
				val = counts_aa[gene][res]
				pile_ups[res] = (val-mean)/(std * np.sqrt(n))
		except:
			continue

	return pile_ups



