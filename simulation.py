import json
import threading
import logging
import Queue
import random
import numpy
from collections import Counter


def simulate(gene, counts, length):
	peaks_max = []
	sim = 0
	while sim < 1E4:
		gene_muts = sum(counts[gene].values())
		pos = []
		for p in range(0, gene_muts):
			pos.append(random.randint(0, length[gene]))
		sim +=1
		counts_pos = Counter(pos)
		peaks_max.append(max(counts_pos.values()))	
	std = numpy.std(peaks_max)
	mean_peak = numpy.mean(peaks_max)
	metric = (numpy.max(counts[gene].values()) - mean_peak) / std
	return mean_peak, std, metric, peaks_max




def do_work(item, counts, length):
	print "Worker running: %s" % item
	result = simulate(item, counts, length)
	#peak_stds[item] = result
	logging.info('\t' + str(item) + ' \t ' + str(result[0]) + '\t' + str(result[1]) + '\t' + str(result[2]))
	return result
	

def run_simulation(dna_element = 'prtn', counts_file = 'counts_aa.json', filename = 'simulation14.log', genes_file = 'genes.json'):
	
	def worker():
		while True:
			print "working..."
			item = q.get()
			do_work(item, counts, length)	
			q.task_done()

	logging.basicConfig(filename=filename,level=logging.DEBUG)

	counts = json.load(open(counts_file))
	length = json.load(open('prtn_len'))
	genes = json.load(open(genes_file))

	if dna_element == 'dna':
		for g in genes:
			if g in length.keys():
				length[g] = length[g] * 3

	q = Queue.Queue()
	for i in range(10):
		t = threading.Thread(target=worker)
    	t.daemon = True
    	t.start()

	for gene in genes:
		if gene in counts.keys() and sum(counts[gene].values()) > 4 and gene in length.keys():
			print "Queuing %s" % gene
			q.put(gene, counts, length)

	q.join()



