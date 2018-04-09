import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from astropy.visualization import astropy_mpl_style
from scipy.optimize import curve_fit

#plt.style.use(astropy_mpl_style)

filepath = 'C://Users/Mason/Desktop/AstroLAB/allsky/'
#filename1 = 'tesCurrentImage.fit'
filename1 = 'testimg.fit'

VAO = EarthLocation(lat=41.6611*u.deg, lon=-91.5302*u.deg, height=225*u.m)
xyimgcenter = (320, 240)

def rawFits(filename):
    f = pyfits.open(filename)
    h = f[0].header
    d = f[0].data
    f.close()
    return h, d

def getTime(header):
    """
    if time.localtime().tm_isdst > 0:
        utcoffset = 5*u.hour
    else:
        utcoffset = 6*u.hour
    """
    t = Time(np.array(header['DATE-OBS']), format='fits', scale='utc')# - utcoffset
    return t

def altradius(coords):
    deltax = abs(xyimgcenter[0] - coords[0])
    deltay = abs(xyimgcenter[1] - coords[1])
    radius = np.sqrt(deltax**2. + deltay**2.)
    return radius


def getangle(line2D):
    line2Dxy = line2D.get_xydata()
    angle = np.degrees(np.arctan2((line2Dxy[1][1] - line2Dxy[0][1]), float(line2Dxy[1][0] - line2Dxy[0][0])))
    return angle


def fitfunc(x, m, b):
    return m*x + b

#####################################################################################
h1 = rawFits(filepath+filename1)[0]
d1 = np.rot90((rawFits(filepath+filename1)[1]))
t1 = getTime(h1)
print(t1)


starxy = [ ['rigel', (518, 337)],
           ['sirius', (465, 421)],
           ['procyon', (450, 283)],
           ['betelgeuse', (423, 374)],
           ['polaris', (167, 251)],
           ['capella', (284, 349)],
           ['castor', (359, 280)],
           ['pollux', (373, 273)],
           ['aldebaran', (360, 426)],
           ['alnilam', (447, 397)]
           ]

d1 = np.flipud(d1)

projection = [ ['horizon', xyimgcenter, 296],
               ['zenith', xyimgcenter, 2]
               ]


colmap=plt.get_cmap('gray')
fig = plt.gcf()
fig.set_size_inches(9,6)
ax = fig.gca()
ax.cla()

altitudes = []
pixelalts = []

azimuths = []
pixelazs = []

test = []

for star in starxy:
    star.append('{0.alt:.3}'.format(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=t1, location=VAO))))
    star.append('{0.az:.3}'.format(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=t1, location=VAO))))
    projection.append([ star[0]+'Alt', star[1], altradius(star[1])])
    ax.add_artist(patches.Circle(star[1], 2, color='b', fill=False))
    meridian = lines.Line2D([xyimgcenter[0], star[1][0]],
                             [xyimgcenter[1], star[1][1]],
                             lw=1, linestyle='--', color='m', axes=ax)
    test.append(meridian)
    ax.add_line(meridian)
    altitudes.append(star[2])
    azimuths.append(star[3])


for marker in projection:
    pixelalts.append(marker[2])
    ax.add_artist(patches.Circle(xyimgcenter, marker[2], color='r', linestyle='--', fill=False, alpha=.25))



#print(projection)
altitudes2 = [19.65, 90.0]
azimuths2 = []

for string in altitudes:
    string=string[:-4]
    altitudes2.append(float(string))

for string in azimuths:
    string=string[:-4]
    azimuths2.append(float(string))


for i in test:
    #print(getangle(i))
    pixelazs.append(getangle(i))

pixelalts = np.asarray(pixelalts)
pixelazs = np.asarray(pixelazs)



#print(starxy)
print()
print('alt (deg): ', altitudes2)
print('alt (px): ', pixelalts)
print()
print('az (deg): ', azimuths2)
print('pixaz (deg): ', pixelazs)



######################### CURVE FITTING ##########################################

poptAlt, pcovAlt = curve_fit(fitfunc, pixelalts[2:], altitudes2[2:])

poptAz, pcovAz = curve_fit(fitfunc, pixelazs, azimuths2)







############################# PLOTTING ###########################################

plt.imshow(d1, cmap = colmap, origin='lower')
#plt.show()

plt.figure(2)
plt.title('Altitude as a function of pixel distance')
plt.xlabel('pixels from zenith')
plt.ylabel('altitude (deg')
plt.ylim(0, 95)
plt.scatter(pixelalts, altitudes2)
plt.plot(pixelalts, fitfunc(pixelalts, *poptAlt), 'r-', ms=8, label='fit: m=%5.3f, b=%5.3f' % tuple(poptAlt))
plt.legend(loc='best')

plt.figure(3)
plt.title("Azimuth as a function of angle")
plt.xlabel('angle from zenith || x-axis (deg)')
plt.ylabel('azimuth (deg)')
plt.xlim(0, 360)
plt.scatter(pixelazs, azimuths2)
plt.plot(pixelazs, fitfunc(pixelazs, *poptAz), 'r-', ms=8, label='fit: m=%5.3f, b=%5.3f' % tuple(poptAz))
plt.legend(loc='best')

plt.show()