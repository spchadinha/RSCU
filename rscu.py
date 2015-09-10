
from codon_table_class import Codon_Table
import subprocess

def make_usage_table(filename):
	with open(filename, "r") as infile, open("rscu_score.txt", "w") as outfile:
		all_tables = []
		headers = []
		for line in infile:
			if ">" in line:
				outfile.write(line)
			elif len(line) < 5:
				pass
			else:
				x = Codon_Table(line)
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
the 's' key to try another search.\nType 'exit' to end the session.\n")

		if userinput == 'r':
			subprocess.call(["open", "rscu_score.txt"])
		elif userinput == 'f':
			subprocess.call(["open", "findings.txt"])
		elif userinput == 'c':
			one = input('Please input the first table\'s index: ')
			two = input('Please input the second table\'s index: ')
			tables[one-1].compare(tables[two-1])
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
		elif userinput == 'exit':
			print "Thank you for using the RSCU calculation tool!"
		else:
			print "Command not recognized."
		print 

main()
