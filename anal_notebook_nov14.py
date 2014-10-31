import operator

execfile('init.py')
from rebuildanal_helper import *

if 'snp_2aa' in os.listdir('./'):
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

	#ax1.axes.get_xaxis().set_visible(False)
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



