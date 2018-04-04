

import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy.visualization import astropy_mpl_style

#plt.style.use(astropy_mpl_style)

filepath = 'C://Users/Mason/Desktop/AstroLAB/allsky/'
filename1 = 'tesCurrentImage.fit'
#filename1 = 'testimg.fit'

VAO = EarthLocation(lat=41.6611*u.deg, lon=-91.5302*u.deg, height=225*u.m)

def rawFits(filename):
    f = pyfits.open(filename)
    h = f[0].header
    d = f[0].data
    f.close()
    return h, d

def getTime(header):
    if time.localtime().tm_isdst > 0:
        utcoffset = 5*u.hour
    else:
        utcoffset = 6*u.hour
    t = Time(np.array(header['DATE-OBS']), format='fits', scale='utc') - utcoffset
    return t

def altradius(coords):
    deltax = abs(300 - coords[0])
    deltay = abs(245 - coords[1])
    radius = np.sqrt(deltax**2. + deltay**2.)
    return radius


#####################################################################################
h1 = rawFits(filepath+filename1)[0]
d1 = rawFits(filepath+filename1)[1]
time = getTime(h1)
print(time)


starxy = [ ['polaris', (459, 246)],
           ['sirius', (102, 270)],
           ['pollux', (253, 219)],
           ['castor', (265, 228)],
           ['betelgeuse', (189, 313)]
           ]


projection = [ ['horizon', (300, 246), 290],
               ['zenith', (300, 246), 2]
               ]


colmap=plt.get_cmap('gray')
fig = plt.gcf()
fig.set_size_inches(9,6)
ax = fig.gca()
ax.cla()

altitudes = []
pixelalts = []

for star in starxy:
    #print(SkyCoord.from_name(star[0]))
    #print(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=time, location=VAO)))
    star.append('{0.alt:.3}'.format(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=time, location=VAO))))
    star.append('{0.az:.3}'.format(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=time, location=VAO))))
    projection.append([ star[0]+'Alt', star[1], altradius(star[1])])
    ax.add_artist(patches.Circle(star[1], 2, color='b', fill=False))
    altitudes.append(star[2])


for marker in projection:
    if marker[1][0] < 300:
        pixelalts.append(300 - marker[1][0])
    else:
        pixelalts.append(marker[2])
    ax.add_artist(patches.Circle((300, 246), marker[2], color='m', linestyle='--', fill=False, alpha=.25))


print(projection)
altitudes2 = [0.00, 90.0]

for string in altitudes:
    string=string[:-4]
    altitudes2.append(float(string))

print(starxy)
print(altitudes2)
print(pixelalts)

plt.imshow(d1, cmap = colmap, origin='lower')
#plt.show()

plt.figure(2)
plt.xlabel('pixels from zenith')
plt.ylabel('degrees')
plt.scatter(pixelalts, altitudes2)
plt.show()
