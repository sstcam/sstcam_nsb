CROTA2 Variable Investigation, Create a NSB fits file for the ASTRI DRACO field using the following parameters:

starttime = Time('2019-05-08T23:37:54.72806') #This absolutely must be in UTC, \
by virtue of NSB failing if you enter in local time.                            
loc = EarthLocation.from_geodetic(lon=14.974609, lat=37.693267, height=1750*u.m\
) #ASTRI Site Coordinates                                                       
obsalt = 73.74*u.deg
obsaz =0.465*u.deg

And for mypycat.txt:
ASTRI DRACO;17 44 36.834;+53 57 16.60;STAR

Then change fits CROTA2 (with HEALPIX option) in ten degree increments until star positions line up with hotspots. 
