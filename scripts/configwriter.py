from astropy.coordinates import EarthLocation, AltAz, ICRS
from astropy.time import Time, TimeDelta
import astropy.units as u
import numpy as np

def generateconfig(observing_time,observing_location,altaz,version='hess_basic',image_size=200,healpix_level=11,gauss=0.0,moon_above_horizon=0.0,sun_below_horizon=-18.0):
    # Generates a dictionary for an astropy observing time and Earth location, alt and az (in astropy altaz form), with other parameters defined.
    
    configuration={'version':version,'image_size':image_size,'healpix_level':healpix_level,'gauss':gauss,'moon_above_horizon':moon_above_horizon,'sun_below_horizon':sun_below_horizon}
    return configuration

def writeconfig(filename,configuration):
    #Actually writes all the elements of a configuration dictionary to a give filename 
    return 0

#paranal=EarthLocation(lat='-24.681546',lon='-24.681546',height=2161.25*u.m)
paranal=EarthLocation.from_geodetic(lon=14.974609, lat=37.693267, height=1750*u.m) #ASTRI Site Coordinates

obsalt=73.21*u.degree
obsaz=0.5*u.degree
print(paranal)

starttime=Time('2019-05-09T01:37:54.72806')
aa=AltAz(alt=obsalt,az=obsaz,location=paranal,obstime=starttime)
coords=aa.transform_to(ICRS())
print(coords.ra.hms,coords.dec.hms)
raise KeyboardInterrupt
dt=TimeDelta(20*60.0*u.s)*np.linspace(0,100)

observing_times=starttime+dt

for obstime in observing_times:
    print(obstime)
    aa=AltAz(alt=obsalt,az=obsaz,location=paranal,obstime=obstime)
    coords=aa.transform_to(ICRS())
    print(coords)
