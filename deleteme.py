


from Bio import SeqIO
from codon_table_class import Codon_Table

codon_db_order = ['CGA', 'CGC', 'CGG', 'CGU', 'AGA', 'AGG', 'CUA', 'CUC', 
'CUG', 'CUU', 'UUA', 'UUG', 'UCA', 'UCC', 'UCG', 'UCU', 'AGC', 'AGU', 'ACA', 
'ACC', 'ACG', 'ACU', 'CCA', 'CCC', 'CCG', 'CCU', 'GCA', 'GCC', 'GCG', 'GCU', 
'GGA', 'GGC', 'GGG', 'GGU', 'GUA', 'GUC', 'GUG', 'GUU', 'AAA', 'AAG', 'AAC', 
'AAU', 'CAA', 'CAG', 'CAC', 'CAU', 'GAA', 'GAG', 'GAC', 'GAU', 'UAC', 'UAU', 
'UGC', 'UGU', 'UUC', 'UUU', 'AUA', 'AUC', 'AUU', 'AUG', 'UGG', 'UAA', 'UAG', 'UGA']

def gen_CUDBlike_table(filename):

	outfile = open(filename[:-4]+"_CUDB.txt", "w")
	codon_count_list = [0]*len(codon_db_order)
	for record in SeqIO.parse(filename, "fasta"):

		sequence = ''
		revcomp = ''
		if 'U' in record.seq:
			sequence = record.seq
			revcomp = record.seq.reverse_complement()
		else:
			sequence = record.seq.transcribe()
			revcomp = record.seq.reverse_complement().transcribe()



		all_coding_seq = ""
		
		for strand, nuc in [(+1, sequence), (-1, revcomp)]:
		    for i in xrange(len(nuc)):
		    	if nuc[i:i+3] == "AUG":
		    		j = i+3
		    		gene = nuc[i:i+3]
		    		while nuc[j:j+3] != "UAG" and nuc[j:j+3] != "UGA" and nuc[j:j+3] != "UAA" and j < len(nuc):
		    			gene += nuc[j:j+3]
		    			j += 3
		    		gene += nuc[j:j+3]
		    		if len(gene) > 300:
		    			all_coding_seq += gene

		i = 0
		while i < len(all_coding_seq):
			codon_count_list[codon_db_order.index(all_coding_seq[i:i+3])] += 1
			i += 3
		

	ugly_string = '    '
	for count in codon_count_list:
		ugly_string += str(count) + ' '

	# bs = Codon_Table(record.id, ugly_string)
	outfile.write(record.id + "\n")
	outfile.write(ugly_string[0:-2]+"\n")

	outfile.close()


gen_CUDBlike_table("GMBV_MARU10962.fna")


	
    			

