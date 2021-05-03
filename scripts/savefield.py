import matplotlib as mpl
mpl.use('Agg')
from astropy.time import Time
from astropy.coordinates import EarthLocation,SkyCoord,AltAz
import astropy.units as u
import copy
import numpy as np
import matplotlib.pyplot as plt
from ctapipe.io import EventSource
from ctapipe.calib import CameraCalibrator
from ctapipe.utils import get_dataset_path
from ctapipe.visualization import ArrayDisplay
from ctapipe.instrument import CameraGeometry
from ctapipe.instrument.optics import OpticsDescription
from ctapipe.coordinates import (GroundFrame,TiltedGroundFrame,NominalFrame,TelescopeFrame,CameraFrame,)
from nsb import config
from nsb.model import nsbModel
from nsb.mypycat import mypycat
from nsb.gaia import Gaia
from nsb.nsbtools import makeDateString, plotMaps
import ephem
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from matplotlib import cm
import sys
from configwriter import generateconfig, writeconfig

#np.set_printoptions(threshold=sys.maxsize)
# Make plots and fonts larger                                                                                                                                                                              
runname=sys.argv[1] # Specify runname at launch

plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16

fovval = 12.0
npixels = 10000 #Number of pixels on fov map
nsky = 2000 #Number of pixels on all sky map
mversion = 'hess_basic'

print('Begin analysis')
con = config.TheConfiguration()

starttime = Time('2019-05-08T23:37:54.72806')
loc = EarthLocation.from_geodetic(lon=14.974609, lat=37.693267, height=1750*u.m) #ASTRI Site Coordinates 
obsalt = 64.1*u.degree
obsaz = 41.65*u.degree
sourcename = 'ASTRI DRACO' # Mypycat Source Name

aa = AltAz(alt=obsalt,az=obsaz,location=loc,obstime=starttime)


conf = generateconfig(starttime,loc,aa)
writeconfig('/home/spencers/sstcam_nsb/configs/'+runname+'.cfg',conf)
con.readConfig('/home/spencers/sstcam_nsb/configs/'+runname+'.cfg') #Hacky solution to avoid config files

dstring = starttime.strftime('%Y/%m/%d %H:%M:%S')
print(dstring)

time_and_date= ephem.Date(dstring)

mpc = mypycat()
source = mpc.get(sourcename)

gaiamap = Gaia(level=11, verbose=True)
model = nsbModel(con, gaiamap, time_and_date, version=mversion, verbose=True)

# draw what you want

model.drawAllSky(size=nsky)
# show the results on screen
plotMaps(model.skymap.data, 'Allskymaps for %s' % (model.observer_source.date))
plt.savefig('/home/spencers/allsky'+runname+'.png',dpi=300)

# draw something else:
fig=plt.figure()
model.drawFOV_source(source=source, fov=fovval, size=npixels)
# show the results on screen
im1=plotMaps(model.skymap.data, 'FOV for %s' % (model.observer_source.date))
# finally call show
#plt.show()
print(model.skymap.data)
# save fits files

hdul = pyfits.HDUList([model.skymap])
hdul.writeto("NSB_of_"+ makeDateString(time_and_date) + "_"+runname+".fits",'exception',True)

plt.savefig('/home/spencers/fov_'+runname+'.png',dpi=300)
