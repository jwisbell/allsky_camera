import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as lines
import astropy.io.fits as pyfits
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from scipy.optimize import curve_fit

#from astropy.visualization import astropy_mpl_style
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

def altaz_to_xy(alt, az, imgcenter=(320., 240.), pxpdeg=-3.3, b=0., az_rot = 0.):    #pxperdeg found from altitude fitting, angleoffset from az fit (b/w x-axis and NORTH)
    x_star = imgcenter[0] + ((90. - alt)*pxpdeg + b)*np.cos(np.deg2rad((az + az_rot)*1.00))
    y_star = imgcenter[1] + ((90. - alt)*pxpdeg + b)*np.sin(np.deg2rad((az + az_rot)*1.00))
    return x_star, y_star

#####################################################################################
h1 = rawFits(filepath+filename1)[0]
d1 = np.rot90((rawFits(filepath+filename1)[1]))
t1 = getTime(h1)
print(t1)


starxy = [ ['sirius', (518, 337)],
           ['rigel', (465, 421)],
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

projection = []


colmap=plt.get_cmap('gray')
fig = plt.gcf()
fig.set_size_inches(9,6)
ax = fig.gca()
ax.cla()
ax.set_clip_on(True)


altitudes = []
pixelalts = []

azimuths = []
pixelazs = []

anglefromline = []

for star in starxy:
    star.append('{0.alt:.3}'.format(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=t1, location=VAO))))
    star.append('{0.az:.3}'.format(SkyCoord.from_name(star[0]).transform_to(AltAz(obstime=t1, location=VAO))))
    projection.append([ star[0]+'Alt', star[1], altradius(star[1])])
    ax.add_artist(patches.Circle(star[1], 2, color='b', fill=False))
    meridian = lines.Line2D([xyimgcenter[0], star[1][0]],
                             [xyimgcenter[1], star[1][1]],
                             lw=1, linestyle='--', color='m', axes=ax)
    anglefromline.append(meridian)
    ax.add_line(meridian)
    altitudes.append(star[2])
    azimuths.append(star[3])


for marker in projection:
    pixelalts.append(marker[2])
    ax.add_artist(patches.Circle(xyimgcenter, marker[2], color='r', linestyle='--', fill=False, alpha=.25))


##### positive x axis & angle annotation ####################################
ax.add_artist(patches.Circle(xyimgcenter, 2, color='y', linestyle='--', fill=False, alpha=.8))    #Zenith
ax.add_artist(patches.Circle(xyimgcenter, 307.36, color='y', linestyle='--', fill=False, alpha=.8))   #Horizon
xaxis = lines.Line2D(xdata=[xyimgcenter[0],640], ydata=[xyimgcenter[1], xyimgcenter[1]+1], color='w', alpha=0.3)
ax.add_artist(xaxis)

#angle_plot = get_angle_plot(xaxis, test[0], 1)
#angle_text = get_angle_text(angle_plot)

#ax.add_patch(angle_plot) # To display the angle arc
#ax.text(*angle_text) # To display the angle value
##############################################################################



#print(projection)
altitudes2 = []
azimuths2 = []

for string in altitudes:
    string=string[:-4]
    altitudes2.append(float(string))

for string in azimuths:
    string=string[:-4]
    azimuths2.append(float(string))

for i in anglefromline:
    #print(getangle(i))
    pixelazs.append(getangle(i))

pixelalts = np.asarray(pixelalts)
pixelazs = np.asarray(pixelazs)
altitudes2 = np.asarray(altitudes2)
azimuths2 = np.asarray(azimuths2)



#print(starxy)
print()
print('alt (deg): ', altitudes2)
print('alt (px): ', pixelalts)
print()
print('az (deg): ', azimuths2)
print('pixaz (deg): ', pixelazs)



######################### CURVE FITTING ##########################################

poptAlt, pcovAlt = curve_fit(fitfunc, altitudes2,  pixelalts)

poptAz, pcovAz = curve_fit(fitfunc, azimuths2, pixelazs)

print(poptAlt, pcovAlt)
print(poptAz, pcovAz)
print()


#### TEST
compDATA = []
forwardDATA = []
testparam = 5.

for i in range(0,10,1):

    tau = altaz_to_xy(alt=altitudes2[i], az=azimuths2[i], imgcenter=xyimgcenter, pxpdeg=np.abs(poptAlt[0]), b=5.0, az_rot=np.abs(poptAz[1]+testparam))
    #print(tau)
    forwardDATA.append(tau)
    compDATA.append(starxy[i][1])
    ax.add_artist(patches.Circle(tau, 2, color='c', linestyle='-', fill=True))

forwardDATA = np.asarray(forwardDATA)
compDATA = np.asarray(compDATA)


############################# PLOTTING ###########################################

plt.figure(1)
plt.title('Preliminary Data Set with Visual Aid')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
#plt.ylim([0, 480])
plt.imshow(d1, cmap = colmap, origin='lower')

plt.figure(2)
plt.title('Altitude as a Function of Pixel Distance')
plt.ylabel('Pixels from zenith')
plt.xlabel('Altitude (deg)')
#plt.ylim(0, 95)
plt.scatter(altitudes2, pixelalts,)
plt.plot(altitudes2, fitfunc(altitudes2, *poptAlt), 'r-', ms=8, label='fit: m=%5.3f, b=%5.3f' % tuple(poptAlt))
plt.legend(loc='best')

plt.figure(3)
plt.title("Azimuth as a Function of Image Angle")
plt.ylabel('Angle from zenith || +x-axis (deg)')
plt.xlabel('Azimuth (deg)')
#plt.xlim(0, 360)
#plt.ylim(0, 360)
plt.scatter(azimuths2, pixelazs)
plt.plot(azimuths2, fitfunc(azimuths2, *poptAz), 'r-', ms=8, label='fit: m=%5.3f, b=%5.3f' % tuple(poptAz))
plt.legend(loc='best')

plt.figure(4)
plt.title('Image XY as a Function of AltAz Coords')
plt.xlabel('AltAz Coords')
plt.ylabel('Image Pixel Coords')
plt.scatter(forwardDATA[:,0], forwardDATA[:,1], label='tau')
plt.scatter(compDATA[:,0], compDATA[:,1], label='true')
plt.legend(loc='best')

plt.show()