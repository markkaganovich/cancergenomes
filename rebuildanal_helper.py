

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
	select_range = filter(lambda x: pos >= x[0] and pos <= x[1], table_ranges)[0]
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
				[current_range, current_table, positions] = get_table(chrom, pos)
				print "changing tables"
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

