import sys
import time
import numpy as np
from subprocess import call
from argparse import ArgumentParser
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
from matplotlib import gridspec
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
from scipy.optimize import curve_fit, leastsq, least_squares

#-------------------
import conversions
import planets
import weather


# -------  constants ---------
HORIZON = 15
MAG_LIM = 3.


#load the test image
def raw_fits(filename):
    #returns FITS header and data
    if filename.lower().endswith(('.fits', '.fit')):
        f = pyfits.open(filename)
        h = f[0].header
        d = f[0].data
        f.close()
        return h, d
    else:
        raise Exception('Invalid Filetype')


#---------- WEATHER STATS -----------------
def weather_graph(im, time,time_str):
	b = weather.brightness(-1*im)
	e = weather.sobel(-1*im,verbose=False)

	#load the previous points
	try:
		edges = np.load('./weather_e.npy')
		edges = np.append(edges, e)
		brights = np.load('./weather_b.npy')
		brights = np.append(brights, b)
		times = np.load('./weather_t.npy')
		times = np.append(times, time)
		time_labs = np.load('./weather_l.npy')
		time_labs = np.append(time_labs, time_str)
		#print 'brightness:', np.min(brights), np.max(brights)
		#print 'edges', np.min(edges), np.max(edges)
	except:
		edges = [e]
		brights = [b]
		times = [time]
		time_labs = [time_str]
	
	if len(edges) >= (12*20):
		edges = np.delete(edges,0)
		brights = np.delete(brights,0)
		times = np.delete(times,0)
		time_labs = np.delete(time_labs,0)
	
	fig,ax = plt.subplots(figsize=(6,8))
	ax.plot(times, brights, 'b--',label='Brightness')
	ax.plot(times, edges, 'r--',label='Edge Percentage (Cloud Cover)')
	ax.set_xticks(times[::5])
	ax.set_xticklabels(time_labs[::5],rotation='vertical')
	plt.savefig('./testing_weather.png')#,bbinches='tight')
	np.save('./weather_e.npy',edges)
	np.save('./weather_b.npy',brights)
	np.save('./weather_t.npy',times)
	np.save('./weather_l.npy',time_labs)
	plt.close()



#------------ THE ALLSKY PIPELINE ------------

def pipe(fname, verbose=True, grid=True, names=True, points=False, do_weather=True, save=False):
	#open the image
	head, im = raw_fits(fname)
	#im = np.rot90(im,1)

	#store the utc time of the obervation, iowa city lat and long
	time_str = head['date-obs']
	#print time_str
	time = Time(time_str) - 1*u.hour
	lst = conversions.utc_to_lst(time_str) - (1./24 * 360)

	#open the catalog
	star_catalog = np.load('star_catalog.npy')
	planet_catalog = planets.mk_cat(fname)

	#using time, lat, and long, convert ra, dec to alt,az
	#alt_az catalog
	sky_cat = []
	for kk,star in enumerate(star_catalog):
		name, star_type, ra, dec, mag, epoch = star
		ra = float(ra); dec = float(dec); mag = float(mag)
		#convert ra, dec to alt, az 
		alt, az = conversions.radec_to_altaz(ra, dec, lst,debug=False)
		if alt > HORIZON and mag < MAG_LIM:
			#convert alt,az to x,y using Mason's fit
			x,y = conversions.altaz_to_xy(alt, az)
			sky_cat.append( [name, x,y, alt, az, ra, dec, mag] )

	for kk,planet in enumerate(planet_catalog):
		name, star_type, ra, dec, mag, epoch = planet
		ra = conversions.sex_to_deg(ra,ra=True); dec = conversions.sex_to_deg(dec, ra=False) ; mag = float(mag)
		#convert ra, dec to alt, az 
		alt, az = conversions.radec_to_altaz(ra, dec, lst,debug=False)
		if alt > HORIZON: 
			#convert alt,az to x,y using Mason's fit
			x,y = conversions.altaz_to_xy(alt, az)
			sky_cat.append( [name, x,y, alt, az, ra, dec, mag] )

	fig, ax = plt.subplots(figsize=(10,12))
	ax.imshow(1*np.sqrt(im), origin='lower',cmap='gray',vmin=.4*np.max(np.sqrt(im)), vmax=1*np.max(np.sqrt(im)))
	#plot objects with names
	for kk, obj in enumerate(sky_cat):
		name = obj[0]; x = obj[1]; y = obj[2]
		if points:
			ax.scatter(y,x,marker='*',color='cyan',alpha=0.5)
		if names:
			if y < 460:
				ax.text(y,x,name, fontsize=8,alpha=0.95, color='gray')

	
	#plot alt,az grid (flag for ra,dec grid)
	do_altaz = False
	do_radec = True
	if grid:
		if do_altaz:
			circs, lines = conversions.altaz_grid(im, spacing=15, verbose=False)
			for c in circs:
				ax.add_patch(c)
			for l in lines:
				ax.plot( l[0],l[1],alpha=0.35,color='red')

		if do_radec:
			color = 'dodgerblue'
			ra_spacings, dec_spacings, ra_filler = conversions.radec_grid(im, lst)
			ralabels = []; dec_labels = []
			for c in ra_spacings:
				imcoords_x, imcoords_y,s = c
				ax.plot(imcoords_y, imcoords_x, c=color, alpha=0.35)
				ralabels.append(s)

				if imcoords_y[0]*1. < 460 and  imcoords_y[0]*1. > 5:
					if imcoords_x[0]*1.07 < 320:
						if imcoords_y[0] < 240:
							imcoords_y[0] -= 15
						if imcoords_y[0]*.99 > 210 and imcoords_y[0]*.99 < 300:
							ax.text(imcoords_y[0]*.99, imcoords_x[0]*1.1 - 50, '%s'%(s) ,color='gray')
						else:
							ax.text(imcoords_y[0]*.99, imcoords_x[0]*1.1 - 30, '%s'%(s) ,color='gray')
					else:
						ax.text(imcoords_y[0]*1., imcoords_x[0]*1.07, '%s'%(s) ,color='gray')

			for c in dec_spacings:
				imcoords_x, imcoords_y, s = c
				ax.plot(imcoords_y, imcoords_x, c=color, alpha=0.35)
				dec_labels.append(s)
				#ax.text(imcoords_y[0], imcoords_x[0], str(s))
			for c in ra_filler:
				imcoords_x, imcoords_y,s = c
				ax.plot(imcoords_y, imcoords_x, c=color, alpha=0.15)
				#ax.text(imcoords_y[-1], imcoords_x[-1], s)

			angles = np.linspace(0,350,len(ralabels))
			offset = 12
			r = 290
			for kk in range(len(ralabels)):
				a = angles[kk]+ offset
				#ax.text(r*np.cos( np.deg2rad(a)) + 480/2, r*np.sin(np.deg2rad(a)) + 640/2,'%.1f'%(ralabels[kk]),color='white',alpha=0.)
		

		

	ax.set_xlim([0,480])
	ax.set_ylim([0,640])
	plt.axis('off')
	fig.subplots_adjust(left=0,right=1,bottom=0,top=1)
	#save as a new file
	outname = fname.split('.FIT')[0] + '.jpg'
	if save:
		plt.savefig(outname, bbinches='tight')
		call('convert %s -crop 890x1200+50+0 +repage %s'%(outname, outname),shell=True) 
	if verbose:
		plt.show()

	#compare to see where nearest stars are
	for kk, obj in enumerate(sky_cat):
		x = obj[1]; y = obj[2]
		try:
			yc,xc = conversions.centroid(np.power(im,1.),x,y,s=20, verbose=False)
		except:
			continue

	#flag for image analysis -- edges, brightness, difference?
	if do_weather:
		weather_graph(im, lst, time_str.split('T')[-1])

	plt.close()
	return outname

def errfunc(p, *args):
	a,b,c,d,e = p
	im, star_catalog, time,ret_all = args
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
			x,y = conversions.altaz_to_xy(alt, az, pxpdeg=a, b= b, imgcenter=(c,d), az_rot=e)
			sky_cat.append( [name, x,y, alt, az, ra, dec, mag] )
	
	#compute the distance from each predicted star to the nearest centroid
	chi2 = 0
	if ret_all:
		chi2 = []
	for kk, obj in enumerate(sky_cat):
		x = obj[1]; y = obj[2]
		try:
			yc,xc = conversions.centroid(np.log(im),x,y,s=20, verbose=False)
			if ret_all:
				chi2.append(np.sqrt( np.power( x-xc,2. ) + np.power(y-yc,2.) )) 
			else:
				chi2 += np.sqrt( np.power( x-xc,2. ) + np.power(y-yc,2.) )
		except:
			continue
	if ret_all:
		return np.array(chi2) / len(sky_cat)
	else:
		chi2 = chi2 / len(sky_cat)
		return np.array([chi2,0,0,0,0])

def min_distances(fname):
	head, im = raw_fits(fname)
	#im = np.rot90(im,1)

	#store the utc time of the obervation, iowa city lat and long
	time_str = head['date-obs']
	lst = conversions.utc_to_lst(time_str) - (1./24 * 360)
	
	#open the catalog
	star_catalog = np.load('star_catalog.npy')

	p0 = [-3.29332406,  -9.5604571 , 316.51919006, 241.93552041, 4.2]
	x = np.arange(480)
	y = np.arange(640)
	xv,yv = np.meshgrid(x,y)
	res = leastsq( errfunc, p0, args=(im, star_catalog, lst,False),full_output=True )

	print res[0]#,res[1]/(640*480.)
	var = errfunc(res[0],im, star_catalog, lst,False)[0]
	print var
	#TODO revisit the mean of the distances once the hardware is complete
	

def make_gif(file_list,outname):
	import imageio
	from subprocess import call
	images = []
	for filename in file_list:
	    images.append(imageio.imread(filename))
	imageio.mimsave(outname, images,quality=5)


if __name__ == '__main__':
	parser = ArgumentParser(description='Generate overlay jpg for given date')
	parser.add_argument('-g', '--grid', type=bool, default=True,help='Insert grid', required=False)
	parser.add_argument('-p', '--points', type=bool, default=False, help='Label Objects', required=False)
	parser.add_argument('-n', '--names', type=bool, default=True, help='names', required=False)
	parser.add_argument('-v', '--verbose', type=bool, default=False, help='verbose output', required=False)
	parser.add_argument('-w', '--weather', type=bool, default=False, help='Weather plotting, default off', required=False)
	parser.add_argument('-m', '--many_files', type=bool, default=False, help='Given a folder, will generate a gif', required=False)
	parser.add_argument('-f', '--filename', type=str, default="", help='filename', required=True)


	args = parser.parse_args()
	filename = args.filename
	grid = args.grid
	points = args.points
	names = args.names
	verbose = args.verbose
	do_weather = args.weather
	many = args.many_files

	if not many:
		pipe(filename, grid=grid,points=points,names=names,verbose=verbose,do_weather=do_weather,save=True)
	else:
		import os
		ims = []
		files = os.listdir(filename)
		for f in files:
			if len(f.split('.FIT')) > 1:
				print (filename+f)
				ims.append( pipe(filename+f, grid=grid,points=points,names=names,verbose=verbose,do_weather=do_weather,save=True) )

		#make a gif from a list of jpgs
		make_gif(ims,filename+'/night.mp4')

	#fname = './assets/clear_with_shield.FIT'
	#min_distances('./assets/clear_with_shield.FIT')





