from astropy.coordinates import EarthLocation, AltAz, ICRS, SkyCoord
from astropy.time import Time
import astropy.units as u
import numpy as np

loc = EarthLocation.from_geodetic(lon=14.974609, lat=37.693267, height=1750*u.m) #ASTRI Site Coordinates                                                                                                   

mainsource = SkyCoord(ra='2 31 49.09',dec='+89 15 50.8',unit=(u.hourangle, u.deg),frame=ICRS) #ASTRO DRACO field position               
starttime=Time('2019-05-09T01:37:54.72806')


print('Source',mainsource.ra.hms,mainsource.dec.hms)
print('Source RA and DEC',mainsource.ra.deg,mainsource.dec.deg)

mainaltaz=mainsource.transform_to(AltAz(location=loc,obstime=starttime))
print('Source alt and az',mainaltaz.alt.deg,mainaltaz.az.deg)

obsalt = 37.3*u.degree
obsaz = 0.2*u.degree

aa=AltAz(alt=obsalt,az=obsaz,location=loc,obstime=starttime)
coords=aa.transform_to(ICRS())
print('Coords',coords.ra.hms,coords.dec.hms)
print(coords.ra.deg,coords.dec.deg)
