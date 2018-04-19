import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from astropy.io import fits
from astropy.time import Time
from matplotlib.patches import Circle
import matplotlib.image as mpimg
#from mpl_toolkits.basemap import Basemap

lat = 41.6611 #North
lon = 91.5302 #West
LAT = np.radians(lat)
LON = np.radians(lon)


def sex_to_deg(sex_string,ra=True,debug=False):
	#ra = True indicates that you are inputing ra
	#ra =  False indicates that you input dec
	if debug:
		print(sex_string)
	chunks = sex_string.split(':')
	h = float(chunks[0])
	m = float(chunks[1])
	s = float(chunks[2])
	if ra:
		ra = (h + m/60. + s/3600.) / 24. * 360
		if debug:
			print(ra)
		return ra
	else:
		dec = h + m/60. + s/3600.
		if debug:
			print(dec)
		return dec


def utc_to_lst(time,hours=False):
	#time is in format given by FITS header
	#need the julian date
	t = Time(time,format='isot',scale='utc')
	jd = t.jd
	#using navy conversion accurate to 0.1 seconds
	d = jd - 2451545.0
	GMST = 18.697374558 + 24.06570982441908*d
	gmst = GMST
	#convert to 24hr clock
	while gmst > 360:
		gmst -= 360

	if hours:
		#convert degrees to hours
		gmst = gmst / 15.
		local_sidereal_time = gmst - lon/15.
		return local_sidereal_time
	return gmst - lon


def radec_to_altaz(ra,dec,time,debug=False):
	if debug:
		print ra, dec,time
	#NOTE -- numpy trig is all in radians
	RA = np.radians(ra)
	DEC = np.radians(dec)
	#hour angle
	H = utc_to_lst(time) - ra #time must be LST

	#now get altitude
	alt = np.arcsin( np.sin(DEC)*np.sin(LAT) + np.cos(DEC)*np.cos(LAT)*np.cos(np.radians(H)) )

	#azumith
	az = np.arccos(  (np.sin(DEC) - np.sin(LAT)*np.sin(alt))/(np.cos(LAT) *np.cos(alt) ) )
	if debug:
		print alt, az

	return np.degrees(alt), np.degrees(az)


def altaz_to_radec(alt,az,time):
	ALT = np.radians(alt)
	AZ = np.radians(az)

	#first the declination
	DEC = np.arcsin( np.sin(ALT)*np.sin(LAT) + np.cos(ALT)*np.cos(LAT)*np.cos(AZ)  )

	#next the hour angle
	H = np.arcsin( -1 * np.sin(AZ)*np.cos(ALT) / np.cos(DEC)  )

	RA = utc_to_lst(time) - H

	return np.degrees(RA), np.degrees(DEC)


def altaz_to_xy(alt, az, imgcenter=(320., 240.), pxpdeg=-3.3, b=0., az_rot = 8.):
    #generalized transformation from Alt/Az (in degrees) to image x, y (in pixels)
    #pxperdeg found from altitude fitting, angleoffset from az fit (b/w x-axis and NORTH)
    x_star = imgcenter[0] + ((90. - alt)*pxpdeg + b)*np.cos(np.deg2rad((az + az_rot)*1.00)) + 10
    y_star = imgcenter[1] + ((90. - alt)*pxpdeg + b)*np.sin(np.deg2rad((az + az_rot)*1.00)) + 35
    return x_star, y_star


def centroid(image, x,y,s, verbose=False):
	#calculate the centroid of the peak within a radius
	x = int(x); y = int(y); s = int(s)
	print x, y, s
	region = image[x-s/2:x+s/2, y-s/2:y+s/2]
	I = []
	J = []
	for i in range(s):
		ival = np.sum(region[i,:])
		jval = np.sum(region[:,i])
		I.append(ival)
		J.append(jval)
	mean_i = np.sum(I) / s
	mean_j = np.sum(J) / s
	I = np.array(I)
	J = np.array(J)
	diff_i = (I - mean_i)[ (I - mean_i) > 0]
	diff_j = (J - mean_j)[ (J - mean_j) > 0]
	xvals = np.arange(0,s)[ (I - mean_i) > 0]
	yvals = np.arange(0,s)[ (J - mean_j) > 0]
	yc = np.sum( diff_i*xvals) / np.sum(diff_i) 
	xc = np.sum( diff_j*yvals) / np.sum(diff_j) 

	if verbose:
		print(xc,yc)
		fig = plt.figure()
		plt.imshow(region,cmap='viridis',origin='lower')
		plt.scatter(xc,yc)
		plt.show()

	return xc+y-s/2,yc+x-s/2



def altaz_grid(im, spacing=10, verbose=False):
	#spacing in degrees
	north = 8
	alts = np.linspace(0,100,spacing)
	azs = np.zeros(alts.shape)#np.linspace(0+north,360+north,spacing)

	#need a circle at each alt, and radial lines at each az
	x,y = altaz_to_xy(alts, azs)
	print x, y

	circs = []
	for radius in np.sqrt( np.power(x-320,2.) + np.power(y-240,2.)):
		c = Circle( (480/2,640/2),radius=radius,fc='none',ec='gray',alpha=0.2)
		circs.append(c)

	#make the lines
	length=np.max( np.sqrt( np.power(x-320,2.) + np.power(y-240,2.)) )
	lines = []
	azs = np.arange(0+north,360+north,15)
	for az in azs:
		lines.append(  ([480/2, 480/2+length*np.cos(np.radians(az))], [640/2,640/2+length*np.sin(np.radians(az))]) )

	#do the plotting
	if verbose:
		fig,ax = plt.subplots()
		ax.imshow(im, origin='upper',cmap='gray')
		for c in circs:
			ax.add_patch(c)

		for l in lines:
			ax.plot( l[0],l[1],alpha=0.2,color='gray')
		ax.set_xlim([0,480])
		ax.set_ylim([0,640])

		plt.show()
	return circs, lines

#def lambert()


def radec_grid(im,time, spacing=10,lat = 41.661):
	ncp = lat
	meridian = utc_to_lst(time)

	coords = []
	for ra in np.arange(0,15,15):
		for dec in np.arange(0,90,15):
			coords.append([meridian,dec])
	im_coords = []
	for c in coords:
		im_coords.append( altaz_to_xy( *radec_to_altaz(c[0],c[1],time)) )

	print im_coords

	fig, ax = plt.subplots()
	ax.imshow(im, origin='lower')

	for c in im_coords:
		x,y = c
		ax.scatter(640/2 + x, 480/2 + y)

	plt.show()

	


if __name__ == "__main__":
	im = np.zeros( (480,640) )


	altaz_grid(im)
	#radec_grid(im, '2018-03-26T08:18:32.410')









