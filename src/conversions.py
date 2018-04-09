import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from astropy.io import fits
from astropy.time import time

lat = 41.6611 #North
lon = 91.5302 #West
LAT = np.radians(lat)
LON = np.radians(lon)


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


def radec_to_altaz(ra,dec,time):
	#NOTE -- numpy trig is all in radians
	RA = np.radians(ra)
	DEC = np.radians(dec)
	#hour angle
	H = utc_to_lst(time) - ra #time must be LST

	#now get altitude
	alt = np.arcsin( np.sin(DEC)*np.sin(LAT) + np.cos(DEC)*np.cos(LAT)*np.cos(np.radians(H)) )

	#azumith
	az = np.arccos(  (np.sin(DEC) - np.sin(LAT)*np.sin(alt))/(np.cos(LAT) *np.cos(alt) ) )

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



