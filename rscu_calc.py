from codon_table_class import Codon_Table




filename = "findings.txt"
def make_usage_table(filename):
	# print 'opening files.....'
	with open(filename, "r") as infile, open("rscu_score.txt", "w") as outfile:
		# print '    files open!'
		all_tables = []
		headers = []
		# print 'reading from files...'
		for line in infile:
			if ">" in line:
				# print '    found header line'
				outfile.write(line)
			else:
				# print '    making codon table....'
				x = Codon_Table(line)#, codon_db_order)
				all_tables.append(x)
				outfile.write(str(x))



make_usage_table(filename)





