

def detect_org():

	searchfile = "all_cu.txt" #raw_input("Please input the name of the file to search: ")
	listfile = "mosquito_list.txt"
	with open(searchfile, "r") as infile, open(listfile, "r") as findall, open("findings.txt", "w") as outfile:

		detect_list = []
		for line in findall.readlines():
			detect_list.append(line.replace("\n", ""))

		codon_usage = infile.readlines()
		line_count = 0
		used = []

		for org in detect_list:
			for line in codon_usage:
				if org in line and org not in used:
					outfile.write(org + " was found in the file!!\n")
					line_count += 1
					used.append(org)

			if line_count == 0:
				outfile.write(org + " not found.\n")
			line_count = 0
			
			
			

detect_org()