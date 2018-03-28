import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import matplotlib.cm as cm
import astropy.io.fits as pyfits
from astropy import units as u

from astropy.time import Time


def rawFits(filename):
    f = pyfits.open(filename)
    h = f[0].header
    d = f[0].data
    f.close
    return h, d

def dispFits(data, cmap='gray'):
    colmap=plt.get_cmap(cmap)
    plt.figure(1)
    plt.imshow(d, cmap = colmap, origin='lower')
    plt.show()
    return

def getLST(header): #Local Sidereal Time from Fits Observation Date: Units of Hour Angle
    t = Time(np.array(header['DATE-OBS']), format='fits', scale='utc')
    t = t.sidereal_time('apparent', longitude='-91.53')
    return t 
