import numpy as np
from astropy.io import fits as pyfits
import ephem
import matplotlib
matplotlib.use('Agg')

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

def mk_cat(fname):
    head, im          = raw_fits(fname)
    time_str          = head['date-obs']
    DateOfObservation = (time_str.replace('-','/')).replace('T',' ')
    icity = ephem.Observer()
    icity.lon = '-91.5302'
    icity.lat = '41.6611'
    icity.elevation = 100
    icity.date = DateOfObservation

    ListOfPlanets = [ephem.Mercury(icity), ephem.Venus(icity), ephem.Mars(icity), ephem.Jupiter(icity), ephem.Saturn(icity), ephem.Uranus(icity), ephem.Neptune(icity)]
    PlanetCatalog = []

    for i in ListOfPlanets:
    	Planet = i
    	#Planet.compute(str(DateOfObservation))
        #p = Planet(icity)
    	PlanetInfo = [str(Planet.name), 'Planet', str(Planet.ra), str(Planet.dec), str(Planet.mag), '2000']
    	PlanetCatalog.append(PlanetInfo)

    PlanetCatalog = np.array(PlanetCatalog)

    return PlanetCatalog