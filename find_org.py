
def count_sum(line):
	line_list = line.split(" ")
	count = 0
	for num in line_list:
		if num != "\n":
			count += int(num)
	return count

def find_org():

	filename = "all_cu.txt" #raw_input("Please input the name of the file to search: ")
	org = "Anopheles" #raw_input("Please input the name of the organism to search for: ")

	with open(filename, "r") as infile, open("findings.txt", "w") as outfile:

		append_next = False
		line_count = 0
		current_headder = ""
		best_header = ""
		codon_counts = 0
		codon_list = ""

		for line in infile.readlines():
			if append_next:
				outfile.write(line)
				# counts = count_sum(line)
				# if counts > codon_counts:
				# 	codon_counts = counts
				# 	codon_list = line
				# 	best_header = current_headder
				append_next = False

			if org in line:
				outfile.write(">" + line)
				# current_headder = line
				append_next = True
				line_count += 2

		# if best_header:
		# 	outfile.write(best_header)
		# 	outfile.write(codon_list)
		# else:
		# 	outfile.write(org + " was not found in " + filename + "!")
			

find_org()

