
from codon_table_class import Codon_Table
import subprocess

def make_usage_table(filename):
	with open(filename, "r") as infile, open("rscu_score.txt", "w") as outfile:
		all_tables = []
		headers = []
		for line in infile:
			if ">" in line:
				outfile.write(line)
				headers.append(line)
			elif len(line) < 5:
				pass
			else:
				x = Codon_Table(headers[len(headers)-1], line)
				all_tables.append(x)
				outfile.write(str(x))
		return all_tables


def find_org(o):

	filename = "all_cu.txt"
	org = o

	with open(filename, "r") as infile, open("findings.txt", "w") as outfile, open("found_codons.txt", "w") as hiddenfile:

		append_next = False
		line_count = 0
		found_org = True

		for line in infile.readlines():
			if append_next:
				hiddenfile.write(line)
				append_next = False

			if org in line:
				line_count += 1
				outfile.write(str(line_count) + ".   >" + line)
				hiddenfile.write(">" + line)
				append_next = True

		if line_count == 0:
			found_org = False
			outfile.write("No match to " + org + " was found.")
		return found_org
				

def main():
	print "Welcome to the RSCU calculation tool!"
	search = raw_input("Please enter the keyword you would like to search: ")
	print 
	results = find_org(search)
	# subprocess.call(["open", "findings.txt"])

	userinput = "get it started in here!"
	while not results and userinput != 'exit':
		retry = raw_input("There were no sequences found that match your keyword, \
please try another search or type 'exit' to quit.\n")
		if retry == 'exit':
			userinput = 'exit'
		else:
			results = find_org(retry)

	tables = make_usage_table('found_codons.txt')

	if userinput != 'exit':
		print "Results were found for your search!"
	else:
		print "Thank you for using the RSCU calculation tool!"

	while userinput != "exit":
		userinput = raw_input("Press the 'f' key to see what sequences were found.\n \
Press the 'r' key to see codon count and RSCU values for each match \
to your keyword.\nPress the 'c' key to compare two sequences' RSCU values \
(input must be sequence's numeric coding, found by chosing the 'f' option).\nPress \
the 's' key to try another search.\nPress 'a' to calculate the CAI for an individual sequence.\n\
Type 'exit' to end the session.\n")

		print 
		if userinput == 'r':
			subprocess.call(["open", "rscu_score.txt"])
		elif userinput == 'f':
			subprocess.call(["open", "findings.txt"])
		elif userinput == 'c':
			one = input('Please input the first table\'s index: ')
			two = input('Please input the second table\'s index: ')
			try:
				tables[one-1].compare(tables[two-1])
			except IndexError:
				print "One of the sequences numbers you selected does not exist in the list!"
		elif userinput == 's':
			newsearch = raw_input("Please input the keyword you would like to search: ")
			results = find_org(newsearch)

			while not results and newsearch != 'exit':
				newsearch = raw_input("There were no sequences found that match your keyword, \
please try another search or type 'exit' to quit.\n")
				if newsearch != 'exit':
					results = find_org(newsearch)

			tables = make_usage_table('found_codons.txt')
			# subprocess.call(["open", "findings.txt"])
		elif userinput == 'a':
			indiv = input("Please select which sequence to calculate CAI: ")
			indiv -= 1
			cai = tables[indiv].CAI(tables)
			print cai
		elif userinput == 'e':
			indiv = input("Please select which sequence to calculate ENC: ")
			indiv -= 1
			print tables[indiv].enc
		elif userinput == 'ef':
			host = input("Please select which sequence is the host: ")
			virus = input("Please selece the virus sequence: ")
			print tables[virus-1].host_usage_effect(tables[host-1])
		elif userinput == 'exit':
			print "Thank you for using the RSCU calculation tool!"
		else:
			print "Command not recognized."
		print 

main()
