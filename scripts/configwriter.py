from astropy.coordinates import EarthLocation, AltAz, ICRS, SkyCoord
from astropy.time import Time, TimeDelta
import astropy.units as u
import numpy as np

def generateconfig(observing_time,observing_location,altaz,version='hess_basic',image_size=200,healpix_level=11,gauss=0.0,moon_above_horizon=0.0,sun_below_horizon=-18.0):
    # Generates a dictionary for an astropy observing time and Earth location, alt and az (in astropy altaz form), with other parameters defined.
    obstime=observing_time.strftime('%Y/%m/%d %H:%M:%S')
    Lon=observing_location.lon.deg
    Lat=observing_location.lat.deg
    elevation=observing_location.height.value
    alt=altaz.alt.deg
    az=altaz.az.deg
    configuration={'time':obstime,
                   'Lon':Lon,
                   'Lat':Lat,
                   'elevation':elevation,
                   'version':version,
                   'image_size':image_size,
                   'alt':alt,
                   'az':az,
                   'healpix_level':healpix_level,
                   'gauss':gauss,
                   'moon_above_horizon':moon_above_horizon,
                   'sun_below_horizon':sun_below_horizon}
    return configuration

def writeconfig(filename,configuration):
    #Actually writes all the elements of a configuration dictionary to a give filename
    f=open(filename,'w')
    for key in configuration.keys():
        f.write(str(key)+' = '+str(configuration[key])+'\n')
    f.close()
    return 0

if __name__ == "__main__":

    #log=EarthLocation(lat='-24.681546',lon='-24.681546',height=2161.25*u.m) # Paranal
    loc=EarthLocation.from_geodetic(lon=14.974609, lat=37.693267, height=1750*u.m) #ASTRI Site Coordinates
    
    obsalt=73.21*u.degree
    obsaz=0.5*u.degree
    
    
    starttime=Time('2019-05-09T02:37:54.72806')
    
    aa=AltAz(alt=obsalt,az=obsaz,location=loc,obstime=starttime)
    conf=generateconfig(starttime,loc,aa)
    writeconfig('testconf.cfg',conf)
    coords=aa.transform_to(ICRS())
    print(coords.ra.hms,coords.dec.hms)
    print(coords.ra.deg,coords.dec.deg)
