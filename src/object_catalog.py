import numpy as np
import matplotlib
matplotlib.use('Agg')

def readfile(filename):
	f = open(filename)
	catalog = []
	lines = f.readlines()
	print lines
	for line in lines:
		if line[0] == '#':
			continue
		if line[0] == '\n':
			continue
		star = [x.strip() for x in line.split(',')]
		if float(star[4]) <= 5.:
			catalog.append(star)
	f.close()
	return catalog

if __name__=="__main__":
	x= readfile('star_catalog.txt')
	from conversions import sex_to_deg
	y = np.copy(x)
	for kk, obj in enumerate(x):
		name, star_type, ra, dec, mag, epoch = obj
		ra_deg = sex_to_deg(ra,ra=True)
		dec_deg = sex_to_deg(dec,debug=True,ra=False)
		y[kk] = [name, star_type, ra_deg, dec_deg, float(mag),int(epoch)]

	print y

	np.save('star_catalog.npy',y)