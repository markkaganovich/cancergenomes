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

braf_bps = map(lambda y: y.start_position, filter(lambda x: x.hugo_symbol == 'BRAF', tcga_rows))

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



