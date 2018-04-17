def readfile(filename):
	fin = open(filename,'r')
	catalog = []
	for line in fin:
		if line[0] == '#':
			continue
		if line[0] == '\n':
			continue
		star = [x for x in line.split(',')]
		if float(star[4]) <= 5.:
			catalog.append(star)
	fin.close()
	return catalog

print(readfile('star_catalog.txt'))