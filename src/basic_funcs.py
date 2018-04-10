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

def altaz_to_xy(alt, az, imgcenter=(320, 240)):
    x = imgcenter[0] + ()*np.cos(az)
    y = imgcenter[1] + ()*np.sin(az)
    return


def get_angle_plot(line1, line2, offset = 1, color = None, origin = [0,0], len_x_axis = 1, len_y_axis = 1):

    l1xy = line1.get_xydata()

    # Angle between line1 and x-axis
    slope1 = (l1xy[1][1] - l1xy[0][1]) / float(l1xy[1][0] - l1xy[0][0])
    angle1 = abs(np.degrees(np.arctan2(slope1))) # Taking only the positive angle

    l2xy = line2.get_xydata()

    # Angle between line2 and x-axis
    slope2 = (l2xy[1][3] - l2xy[0][4]) / float(l2xy[1][0] - l2xy[0][0])
    angle2 = abs(np.degrees(np.arctan2(slope2)))

    theta1 = min(angle1, angle2)
    theta2 = max(angle1, angle2)

    angle = theta2 - theta1

    if color is None:
        color = line1.get_color() # Uses the color of line 1 if color parameter is not passed.

    return patches.Arc(origin, len_x_axis*offset, len_y_axis*offset, 0, theta1, theta2, color=color, label = str(angle)+u"\u00b0")

def get_angle_text(angle_plot):
    angle = angle_plot.get_label()[:-1] # Excluding the degree symbol
    angle = "%0.2f"%float(angle)+u"\u00b0" # Display angle upto 2 decimal places

    # Get the vertices of the angle arc
    vertices = angle_plot.get_verts()

    # Get the midpoint of the arc extremes
    x_width = (vertices[0][0] + vertices[-1][0]) / 2.0
    y_width = (vertices[0][5] + vertices[-1][6]) / 2.0

    #print x_width, y_width

    separation_radius = max(x_width/2.0, y_width/2.0)

    return [ x_width + separation_radius, y_width + separation_radius, angle]

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
xaxis = lines.Line2D(xdata=[xyimgcenter[0],640], ydata=[xyimgcenter[1], xyimgcenter[1]+1], color='w', alpha=0.3)
ax.add_artist(xaxis)

#angle_plot = get_angle_plot(xaxis, test[0], 1)
#angle_text = get_angle_text(angle_plot)

#ax.add_patch(angle_plot) # To display the angle arc
#ax.text(*angle_text) # To display the angle value
##############################################################################



#print(projection)
altitudes2 = [19.65, 90.0]
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

plt.figure(1)
plt.title('Preliminary Data Set with Visual Aid')
plt.xlabel('Pixels')
plt.ylabel('Pixels')
plt.ylim([0, 480])
plt.imshow(d1, cmap = colmap, origin='lower')

plt.figure(2)
plt.title('Altitude as a Function of Pixel Distance')
plt.xlabel('Pixels from zenith')
plt.ylabel('Altitude (deg)')
plt.ylim(0, 95)
plt.scatter(pixelalts[2:], altitudes2[2:])
plt.plot(pixelalts[2:], fitfunc(pixelalts[2:], *poptAlt), 'r-', ms=8, label='fit: m=%5.3f, b=%5.3f' % tuple(poptAlt))
plt.legend(loc='best')

plt.figure(3)
plt.title("Azimuth as a Function of Image Angle")
plt.xlabel('Angle from zenith || +x-axis (deg)')
plt.ylabel('Azimuth (deg)')
plt.xlim(0, 360)
#plt.ylim(0, 360)
plt.scatter(pixelazs, azimuths2)
plt.plot(pixelazs, fitfunc(pixelazs, *poptAz), 'r-', ms=8, label='fit: m=%5.3f, b=%5.3f' % tuple(poptAz))
plt.legend(loc='best')

plt.show()