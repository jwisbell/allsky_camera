import numpy as np
from astropy.io import fits as pyfits
import ephem

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

head, im          = raw_fits('025007.fit')
time_str          = head['date-obs']
DateOfObservation = (time_str.replace('-','/')).replace('T',' ')
print(DateOfObservation)

ListOfPlanets = [ephem.Mercury(), ephem.Venus(), ephem.Mars(), ephem.Jupiter(), ephem.Saturn(), ephem.Uranus(), ephem.Neptune()]
PlanetCatalog = []

for i in ListOfPlanets:
	Planet = i
	Planet.compute(str(DateOfObservation))
	PlanetInfo = [str(Planet.name), 'Planet', str(Planet.ra), str(Planet.dec), str(Planet.mag), '2000']
	PlanetCatalog.append(PlanetInfo)

PlanetCatalog = np.array(PlanetCatalog)

print(PlanetCatalog)