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

#np.set_printoptions(threshold=sys.maxsize)
# Make plots and fonts larger                                                                                                                                                                              
runname=sys.argv[1]
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16

fovval=12.0

print('Begin analysis')
con = config.TheConfiguration()
#con.readStandardConfig()
con.readConfig("/home/spencers/astri_config.cfg")

time_and_date = ephem.Date("2019/05/09 00:37:54.728026")

mpc = mypycat()
source = mpc.get("ASTRI DRACO")

gaiamap = Gaia(level=11, verbose=True)
model = nsbModel(con, gaiamap, time_and_date, version="hess_basic", verbose=True)

# draw what you want

model.drawAllSky(size=800)
# show the results on screen
plotMaps(model.skymap.data, 'Allskymaps for %s' % (model.observer_source.date))
plt.savefig('/home/spencers/allsky'+runname+'.png')

# draw something else:
fig=plt.figure()
model.drawFOV_source(source=source, fov=fovval, size=10000)
# show the results on screen
im1=plotMaps(model.skymap.data, 'FOV for %s' % (model.observer_source.date))
# finally call show
#plt.show()
print(model.skymap.data)
# save fits files

hdul = pyfits.HDUList([model.skymap])
hdul.writeto("NSB_of_"+ makeDateString(time_and_date) + "_"+runname+".fits",'exception',True)

plt.savefig('/home/spencers/fov_'+runname+'.png',dpi=300)
