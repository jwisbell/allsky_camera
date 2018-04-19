import sys
import time
import numpy as np
from subprocess import call
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import gridspec
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
from scipy.optimize import curve_fit, leastsq, least_squares

#-------------------
import conversions


# -------  constants ---------
HORIZON = 10
MAG_LIM = 2.5


#load the test image
def rawFits(filename):
    f = fits.open(filename)
    h = f[0].header
    d = f[0].data
    f.close()
    return h, d

#------------ THE ALLSKY PIPELINE ------------

def pipe(fname, verbose=True, grid=True, names=True, points=False, save=False):
	#open the image
	head, im = rawFits(fname)
	#im = np.rot90(im,1)

	#store the utc time of the obervation, iowa city lat and long
	time = head['date-obs']
	print time

	#open the catalog
	star_catalog = np.load('star_catalog.npy')
	#print star_catalog[0]

	#using time, lat, and long, convert ra, dec to alt,az
	#alt_az catalog
	sky_cat = []
	for kk,star in enumerate(star_catalog):
		name, star_type, ra, dec, mag, epoch = star
		ra = float(ra); dec = float(dec); mag = float(mag)

		#convert ra, dec to alt, az 
		alt, az = conversions.radec_to_altaz(ra, dec, time,debug=False)
		if alt > HORIZON and mag < MAG_LIM:
			#convert alt,az to x,y using Mason's fit
			x,y = conversions.altaz_to_xy(alt, az)
			sky_cat.append( [name, x,y, alt, az, ra, dec, mag] )


	fig, ax = plt.subplots()
	ax.imshow(-1*np.sqrt(im), origin='lower',cmap='gray')#,vmax=0, vmin=-0.5*np.max(np.sqrt(im)))
	#plot objects with names
	for kk, obj in enumerate(sky_cat):
		name = obj[0]; x = obj[1]; y = obj[2]
		if points:
			ax.scatter(y,x,marker='*',color='cyan',alpha=0.5)
		if names:
			ax.text(y,x,name, fontsize=5)

	
	#plot alt,az grid (flag for ra,dec grid)
	if grid:
		circs, lines = conversions.altaz_grid(im, spacing=15, verbose=False)

		for c in circs:
				ax.add_patch(c)

		for l in lines:
			ax.plot( l[0],l[1],alpha=0.2,color='gray')

	ax.set_xlim([0,480])
	ax.set_ylim([0,640])
	if verbose:
		plt.show()

	#compare to see where nearest stars are
	for kk, obj in enumerate(sky_cat):
		x = obj[1]; y = obj[2]
		#try:
		yc,xc = conversions.centroid(np.power(im,2.),x,y,s=20, verbose=True)

	#flag for image analysis -- edges, brightness, difference?

	#save as a new file
	if save:
		plt.savefig(outname, bbinches='tight')

def errfunc(p, *args):
	a,b,c,d = p
	im, star_catalog, time = args
	#star_catalog, time = args
	#using time, lat, and long, convert ra, dec to alt,az
	#alt_az catalog
	sky_cat = []
	for kk,star in enumerate(star_catalog):
		name, star_type, ra, dec, mag, epoch = star
		ra = float(ra); dec = float(dec); mag = float(mag)

		#convert ra, dec to alt, az 
		alt, az = conversions.radec_to_altaz(ra, dec, time,debug=False)
		if alt > HORIZON and mag < MAG_LIM:
			#convert alt,az to x,y using Mason's fit
			x,y = conversions.altaz_to_xy(alt, az, pxpdeg=a, b= b, imgcenter=(c,d))
			print x, y
			sky_cat.append( [name, x,y, alt, az, ra, dec, mag] )
	
	#compute the distance from each predicted star to the nearest centroid
	chi2 = 0
	for kk, obj in enumerate(sky_cat):
		x = obj[1]; y = obj[2]
		try:
			yc,xc = conversions.centroid(im,x,y,s=30, verbose=False)
			chi2 += np.power( x-xc,2. ) + np.power(y-yc,2.)
		except:
			print x,y
			return 30
	chi2 = chi2 / len(sky_cat)
	return np.ones(4)*chi2




def min_distances(fname):
	head, im = rawFits(fname)
	#im = np.rot90(im,1)

	#store the utc time of the obervation, iowa city lat and long
	time = head['date-obs']
	
	#open the catalog
	star_catalog = np.load('star_catalog.npy')

	p0 = [-3.3, 0, 320., 240.]
	x = np.arange(480)
	y = np.arange(640)
	xv,yv = np.meshgrid(x,y)
	#p, cov = curve_fit(errfunc, (xv,yv), im, p0, args=(star_catalog, time))
	res = leastsq( errfunc, p0, args=(im,star_catalog, time) )

	print res

	




if __name__ == '__main__':
	pipe('./assets/clear_with_shield.FIT', grid=True)
	#min_distances('./assets/clear_with_shield.FIT')





