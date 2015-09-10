aa_to_codon = {'I':['AUA', 'AUC', 'AUU'], 'M':['AUG'], 'T':['ACA', 'ACC', 'ACG', 'ACU'],
	'N':['AAC', 'AAU'], 'K':['AAA', 'AAG'], 'S':['AGC', 'AGU'], 'R':['AGA', 'AGG'],
	'L':['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'], 'P':['CAA', 'CCC', 'CCG', 'CCU'], 'H':['CAC', 'CAU'],
	'Q':['CAA', 'CAG'], 'R':['CGG', 'CGU', 'CGA', 'CGC'], 'V':['GUA', 'GUC', 'GUG', 'GUU'],
	'A':['GCA', 'GCC', 'GCG', 'GCU'], 'D':['GAC', 'GAU'], 'E':['GAA', 'GAG'],
	'G':['GGA', 'GGC', 'GGG', 'GGU'], 'S':['UAC', 'UCC', 'UCG', 'UCU'],
	'F':['UUC', 'UUU'], 'Y':['UAC', 'UAU'], 'C':['UGC', 'UGU'],
	'W':['UGG'], '*':['UAA', 'UAG', 'UGA']}

codon_to_aa = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'*', 'UAG':'*',
    'UGC':'C', 'UGU':'C', 'UGA':'*', 'UGG':'W',
    }

codon_db_order = ['CGA', 'CGC', 'CGG', 'CGU', 'AGA', 'AGG', 'CUA', 'CUC', 
'CUG', 'CUU', 'UUA', 'UUG', 'UCA', 'UCC', 'UCG', 'UCU', 'AGC', 'AGU', 'ACA', 
'ACC', 'ACG', 'ACU', 'CCA', 'CCC', 'CCG', 'CCU', 'GCA', 'GCC', 'GCG', 'GCU', 
'GGA', 'GGC', 'GGG', 'GGU', 'GUA', 'GUC', 'GUG', 'GUU', 'AAA', 'AAG', 'AAC', 
'AAU', 'CAA', 'CAG', 'CAC', 'CAU', 'GAA', 'GAG', 'GAC', 'GAU', 'UAC', 'UAU', 
'UGC', 'UGU', 'UUC', 'UUU', 'AUA', 'AUC', 'AUU', 'AUG', 'UGG', 'UAA', 'UAG', 'UGA']

class Codon_Table():
	"""docstring for ClassName"""
	def __init__(self, codon_count):
		#super(ClassName, self).__init__()
		#self.arg = arg
		self.codon_count = codon_count.replace('    ', '').split(' ')
		self.codon_order = codon_db_order
		self.codon_table = self._make_table()
		self.aa_count = self._count_aa()
		self.rscu_table = self._rscu()

	def _make_table(self):
		table = {}
		for i in xrange(len(self.codon_order)):
			table[self.codon_order[i]] = int(self.codon_count[i])
		return table

	def _count_aa(self):
		table = {}
		for codon in self.codon_table:
			aa = codon_to_aa[codon]
			if aa not in table:
				table[aa] = int(self.codon_table[codon])
			else:
				table[aa] += int(self.codon_table[codon])
		return table

	def _rscu(self):
		table = {}
		for codon in self.codon_table:
			aa = codon_to_aa[codon]
			total = 0.0
			for sym in aa_to_codon[aa]:
				total += self.codon_table[sym]
			if total != 0:
				table[codon] = (self.codon_table[codon]/total)*len(aa_to_codon[aa])	
			else:
				table[codon] = 0
		
		# for a in aa_to_codon:
		# 	print a
		# 	for c in aa_to_codon[a]:
		# 		print c, ":", str(table[c]) + ",",
		# 	print 
		# print table
		return table


	def compare(self, other):
		total = 0
		for aa in self.rscu_table:
			total += self.rscu_table[aa] - other.rscu_table[aa]
			print aa, self.rscu_table[aa] - other.rscu_table[aa]
		print "Difference:", total 

	def __str__(self):
		string = ""

		for aa in aa_to_codon:
			string += "    Amino Acid: " + aa + "  Count: " + str(self.aa_count[aa]) +"\n"
			string += "    Codon Usage: "
			for codon in aa_to_codon[aa]:
				string += codon + " : " + str(self.codon_table[codon]) + "  "
			string += "\n    RSCU Values: "
			for codon in aa_to_codon[aa]:
				string += codon + " : " + str(self.rscu_table[codon]) + "  "
			string += "\n\n"

		return string 












