

with open("all_cu.txt", "r") as infile, open("sorted_orgs.txt", "w") as outfile:
	lst = []
	for line in infile.readlines():
		if ":" in line:
			x = line[line.find(":")+1:]
			y = x[:x.find(":")] + "\n"
			lst.append(y)
	slst = lst.sort()
	for line in lst:
		outfile.write(line)

