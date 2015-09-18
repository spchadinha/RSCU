
from scipy.stats import gmean
import math
from Bio import SeqIO


# dictionary that maps single-letter amino acid codes to codons that code for the amino acid
aa_to_codon = {'I':['AUA', 'AUC', 'AUU'], 'M':['AUG'], 'T':['ACA', 'ACC', 'ACG', 'ACU'],
	'N':['AAC', 'AAU'], 'K':['AAA', 'AAG'], 'S':['AGC', 'AGU'], 'R':['AGA', 'AGG'],
	'L':['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'], 'P':['CAA', 'CCC', 'CCG', 'CCU'], 'H':['CAC', 'CAU'],
	'Q':['CAA', 'CAG'], 'R':['CGG', 'CGU', 'CGA', 'CGC'], 'V':['GUA', 'GUC', 'GUG', 'GUU'],
	'A':['GCA', 'GCC', 'GCG', 'GCU'], 'D':['GAC', 'GAU'], 'E':['GAA', 'GAG'],
	'G':['GGA', 'GGC', 'GGG', 'GGU'], 'S':['UAC', 'UCC', 'UCG', 'UCU'],
	'F':['UUC', 'UUU'], 'Y':['UAC', 'UAU'], 'C':['UGC', 'UGU'],
	'W':['UGG'], '*':['UAA', 'UAG', 'UGA']}

# dictionary that maps codons to the amino acid they code for
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

# ordered list of codons that represents the order in which codons are referenced
# in the codon usage database FTP files
codon_db_order = ['CGA', 'CGC', 'CGG', 'CGU', 'AGA', 'AGG', 'CUA', 'CUC', 
'CUG', 'CUU', 'UUA', 'UUG', 'UCA', 'UCC', 'UCG', 'UCU', 'AGC', 'AGU', 'ACA', 
'ACC', 'ACG', 'ACU', 'CCA', 'CCC', 'CCG', 'CCU', 'GCA', 'GCC', 'GCG', 'GCU', 
'GGA', 'GGC', 'GGG', 'GGU', 'GUA', 'GUC', 'GUG', 'GUU', 'AAA', 'AAG', 'AAC', 
'AAU', 'CAA', 'CAG', 'CAC', 'CAU', 'GAA', 'GAG', 'GAC', 'GAU', 'UAC', 'UAU', 
'UGC', 'UGU', 'UUC', 'UUU', 'AUA', 'AUC', 'AUU', 'AUG', 'UGG', 'UAA', 'UAG', 'UGA']

class Codon_Table():
	"""
	Represents a codon table from given codon usage information (taken from the 
	codon usage database)

	Parameters: 
		- name : String - a name for the table being created
		- codon_count : String - a string of all the codon count values for a given sequence
	"""
	def __init__(self, name, codon_count):
		# format line of input to make list of codon usage
		self.codon_count = codon_count.replace('    ', '').split(' ')
		# used ordered list of codons as reference
		self.codon_order = codon_db_order
		# assign the object a name
		self.name = name
		# make codon usage table 
		self.codon_table = self._make_table()
		# make amino acid usage table
		self.aa_count = self._count_aa()
		# calcualte RSCU values for each codon
		self.rscu_table = self._rscu()
		# calculate the ENC for the given table
		self.enc = self._enc_calc()
		# perform nucleotide composition analysis 
		self._comp_analysis() #self.per_comp, self.pos_comp = 

	# def _find_orf(self, fastafile):
	# 	record = SeqIO.read("NC_005816.fna", "fasta")

	# 	for strand, nuc in [(+1, record.seq), (-1, record.seq.reverse_complement())]:
	# 		all_coding_seq = ""
	# 		for i in xrange(len(nuc)):
	# 	    	if nuc[i:i+3] == "AUG":
	# 	    		j = i+3
	# 	    		gene = nuc[i:i+3]
	# 	    		while nuc[j:j+3] != "UAG" and nuc[j:j+3] != "UGA" and nuc[j:j+3] != "UAA" and j < len(nuc):
	# 	    			gene += nuc[j:j+3]
	# 	    			j += 3
	# 	    			count += 1
	# 	    		gene += nuc[j:j+3]
	# 	    		if len(gene) > 300:
	# 	    			all_coding_seq += gene

	# 	    #print all_coding_seq


	def _make_table(self):
		"""
		Given a codon ordering and list of integers, returns a dictionary mapping 
		codons to their frequency

		Returns:
			- table : {String:Int} - a dictionary mapping codons to their count
		"""
		table = {}
		for i in xrange(len(self.codon_order)):
			table[self.codon_order[i]] = int(self.codon_count[i])
		return table

	def _count_aa(self):
		"""
		Given a codon usage table and codon to amino acid mapping, returns a dictionary
		mapping amino acids to their count

		Returns:
			- table : {String:Int} - a dictionary mapping amino acids to their count
		"""
		table = {}
		for codon in self.codon_table:
			aa = codon_to_aa[codon]
			if aa not in table: # check to see if the amino acid has already been added to the table
				table[aa] = int(self.codon_table[codon]) # if so, update the existing count
			else: # if the amino acid hasn't been added to the table
				table[aa] += int(self.codon_table[codon]) # add the amino acid to the table and assign a count
		return table

	def _rscu(self):
		"""
		Given a codon usage table and codon to amino acid mapping, returns a dictionary
		mapping codons to their Relative Synonymous Codon Usage value

		Returns:
			- table : {String:Float} - a dictionary mapping codons to their RSCU values 
		"""
		table = {}
		for codon in self.codon_table:
			aa = codon_to_aa[codon]
			if self.aa_count[aa] != 0:
				# RSCU calculation: divide codon count by number of occurences of the amino acid it 
				# represents, then multiply by number of codons that code for the amino acid
				table[codon] = (float(self.codon_table[codon])/self.aa_count[aa])*len(aa_to_codon[aa])	
			else:
				# if the amino acid is never used, assign 0 as RSCU value
				table[codon] = 0
		return table


	def compare(self, other):
		"""
		Given another Codon_Table object, makes a simple, non-statistical comparison
		between this object and the given object

		Parameters:
			- other : Codon_Table - a Codon_Table object to which the current object will
			is compared to
		"""
		total = 0
		for aa in self.rscu_table:
			total += self.rscu_table[aa] - other.rscu_table[aa]
			# print the difference between individual rscu values
			print aa, self.rscu_table[aa] - other.rscu_table[aa]
		# print the total difference between all rscu values
		print "Difference:", total 

	def CAI(self, ref_set):
		"""
		Given a reference set of Codon_Table objects, calculates the Codon Adaptation Index
		for the current Codon_Table object

		Parameters:
			- ref_set : [Codon_Table] - a list of Codon_Table objects from which the ACI is 
			calculated for the current object

		Returns:
			- cai : Float - the CAI score for the current object and given reference set
		"""
		w_lst = []
		count = 0
		max_table = self._max_rscu_table(ref_set)
		for codon in self.rscu_table:
			if self.rscu_table[codon] > 0:
				count += 1
				if max_table[codon] > 0:
					w_lst.append(self.rscu_table[codon]/max_table[codon])
				else:
					w_lst.append(0)
		cai = gmean(w_lst)#/count
		return cai

	def _max_rscu_table(self, ref_set):
		"""
		Given a set of Codon_Table objects, returns the maximum RSCU value for each codon
		as a dictionary mapping codon to max value

		Parameters:
			- ref_set : [Codon_Table] - list of Codon_Table objects 

		Returns:
			- max_table : {String:Float} - dictionary mapping codons to the max RSCU value 
			from the input set
		"""
		max_table = {}
		for codon in ref_set[0].rscu_table:
			rscu_list = []
			for ref in ref_set:
				rscu_list.append(ref.rscu_table[codon])
			max_table[codon] = max(rscu_list)
		return max_table

	def _enc_calc(self):
		"""
		Calculates the effective number of codons from a codon usage table

		Returns:
			- enc : Float - the effective number of codons 
		"""

		indiv_fk = {}
		# generate the sum of rscu values squared for each amino acid
		for aa in aa_to_codon:
			if aa != '*' and aa != 'W' and aa != 'M': # excluding stop codons and aa with one codon
				accum = 0
				num_aa = self.aa_count[aa] 
				for codon in aa_to_codon[aa]:
					# use pseudocount to control for the absence of an amino acid family in small sequences
					accum += ((float(self.codon_table[codon])+1)/(num_aa + len(aa_to_codon[aa])))**2

				if num_aa > 1:
					indiv_fk[aa] = (num_aa*accum - 1)/(num_aa - 1)
				elif num_aa == 0:
					indiv_fk[aa] = 0
				else:
					indiv_fk[aa] = (num_aa*accum - 1)/(num_aa)
		
		# find the mean of the sums above for each amino acid family
		mean_fk = {2:[], 3:[], 4:[], 6:[]}
		for fk_aa in indiv_fk:
			mean_fk[len(aa_to_codon[fk_aa])].append(indiv_fk[fk_aa])

		# compute the ENC from the mean of rscu sums
		fk2 = float(sum(mean_fk[2]))/len(mean_fk[2])
		fam2 = 9
		if fk2 == 0: fam2 = 0
		else: fam2 /= fk2 
		fk3 = float(sum(mean_fk[3]))/len(mean_fk[3])
		fam3 = 1
		if fk3 == 0: fam3 = 0
		else: fam3 /= fk3
		fk4 = float(sum(mean_fk[4]))/len(mean_fk[4])
		fam4 = 5
		if fk4 == 0: fam4 = 0
		else: fam4 /= fk4
		fk6 = float(sum(mean_fk[6]))/len(mean_fk[6])
		fam6 = 3
		if fk6 == 0: fam6 = 0
		else: fam6 /= fk6

		enc = 2 + fam2 + fam3 + fam4 + fam6
		return enc

	def host_usage_effect(self, host):
		"""
		Estimates the effects of host codon usage on an individual's codon usgae

		Parameters:
			- host : Codon_Table - list of Codon_Table objects that represent the individual's host(s)

		Returns:
			- dab : Float - float value between 0.0 and 1.0 representing the host's influence on individual's codon usage
		"""

		num_sum = 0
		a_sum = 0
		b_sum = 0

		for codon in self.rscu_table:
			if codon_to_aa[codon] != 'M' and codon_to_aa[codon] != 'W' and codon_to_aa[codon] != '*':
				num_sum += self.rscu_table[codon] * host.rscu_table[codon]
				a_sum += self.rscu_table[codon]**2
				b_sum += host.rscu_table[codon]**2
		rab = num_sum/math.sqrt(a_sum * b_sum)
		dab = (1-rab)/2
		return dab


	def _comp_analysis(self):
		"""
		Using the codon usage data generates statistics about the composition of the sequence
		eg. percent A, U, C, and G and other like statistics

		Returns:

		"""
		table = {1:{'A':0, 'U':0, 'C':0, 'G':0}, 2:{'A':0, 'U':0, 'C':0, 'G':0}, 3:{'A':0, 'U':0, 'C':0, 'G':0}}
		for codon in self.codon_table:
			count = 1
			for n in codon:
				table[count][n] += self.codon_table[codon]
				count += 1
		
		per_comp = {}
		pos_comp = {1:{'A':0, 'U':0, 'C':0, 'G':0}, 2:{'A':0, 'U':0, 'C':0, 'G':0}, 3:{'A':0, 'U':0, 'C':0, 'G':0}}
		total = 0.0
		for i in table:
			pos_count = 0.0
			for n in table[i]:
				pos_count += table[i][n]
				if n in per_comp:
					per_comp[n] += table[i][n]
					total += table[i][n]
				else:
					per_comp[n] = table[i][n]
					total += table[i][n]
			for n in pos_comp[i]:
				pos_comp[i][n] = table[i][n]/pos_count

		for p in per_comp:
			per_comp[p] = per_comp[p]/total

		# print "Percentage composition"
		# print per_comp
		# print "Percent composition of each nucleotied at each position"
		# print pos_comp
		# print 
		
	def add(self, other):
		for codon in self.codon_table:
			self.codon_table[codon] += other.codon_table[codon]

		self.aa_count = self._count_aa()
		self.rscu_table = self._rscu()
		self.enc = self._enc_calc()



	def __str__(self):
		"""
		Override the __str__, now prints the amino acid, all codons that code for it and their
		count and RSCU values. Does this for each amino acid.
		"""
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

	def __eq__(self, other):
		same = True
		for codon in self.codon_table:
			if self.codon_table[codon] != other.codon_table[codon]:
				same = False
		if self.name != other.name:
			same = False
		return same














