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

# Make plots and fonts larger                                                                                                                                                                               
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16

fovval=10.0

def customMap(image, title,fovval):
    plt.figure()
    vmin, vmax = np.nanmin(image), np.minimum(4 * np.nanmedian(image), np.nanmax(image))
    imshape=image.shape[0]
    #plt.imshow(image, cmap=cm.viridis, vmin=vmin, vmax=vmax,extent=[-image.shape[1]/2., image.shape[1]/2., -image.shape[0]/2., image.shape[0]/2. ])
    im1=plt.imshow(image, cmap=cm.viridis, vmin=vmin, vmax=vmax,extent=[-fovval/2.0, fovval/2.0, -fovval/2.0, fovval/2.0])
    plt.colorbar(label='Brightness (nLb)')
    #plt.gca().get_xaxis().set_visible(False)
    #plt.gca().axes.get_yaxis().set_visible(False)
    plt.title(title)
    plt.draw()
    return im1


con = config.TheConfiguration()
#con.readStandardConfig()
con.readConfig("/home/spencers/my_config.cfg")

time_and_date = ephem.Date("2022/04/02 00:30:00")

mpc = mypycat()
source = mpc.get("Vela Pulsar")

gaiamap = Gaia(level=8, verbose=True)
model = nsbModel(con, gaiamap, time_and_date, version="hess_basic", verbose=True)

# draw what you want
'''
model.drawAllSky(size=800)
# show the results on screen
plotMaps(model.skymap.data, 'Allskymaps for %s' % (model.observer_source.date))
print(model.skymap.data)
plt.savefig('/home/spencers/allsky.png')
'''
# draw something else:
fig=plt.figure()
model.drawFOV_source(source=source, fov=fovval, size=2048)
# show the results on screen
im1=customMap(model.skymap.data, 'FOV for %s' % (model.observer_source.date),fovval)
# finally call show
#plt.show()

# save fits files
#hdul = pyfits.HDUList([model.skymap])
#hdul.writeto("NSB_of_"+ makeDateString(time_and_date) + "_.fits", 'exception', True)


# Set observatory location and time
location = EarthLocation.from_geodetic(-70.317876,-24.681546)
obstime = Time('2022-04-02T00:30')

vela = SkyCoord.from_name("vela pulsar")

altaz = AltAz(location=location, obstime=obstime)

pointing=vela.transform_to(altaz)

cam = CameraGeometry.from_name('CHEC') #Note this is very old, needs to be updated for prod 5.
pix_x = cam.pix_x
pix_y = cam.pix_y

optics = OpticsDescription.from_name('SST-ASTRI') #As is this
focal_length = optics.equivalent_focal_length

camera_frame = CameraFrame(focal_length=focal_length,rotation=0*u.deg,telescope_pointing=pointing)

cam_coords = SkyCoord(pix_x,pix_y,frame=camera_frame)

telescope_frame = TelescopeFrame(
    telescope_pointing=pointing,
    obstime=pointing.obstime,
    location=pointing.location,
)
telescope_coords = cam_coords.transform_to(telescope_frame)
plt.savefig('/home/spencers/fov.png',dpi=300)
fig=plt.figure()
#wrap_angle = telescope_pointing.az + 180* u.deg
print(np.shape(telescope_coords.fov_lon.deg))
print(np.shape(model.skymap.data))
plt.axis('equal')
plt.scatter(
    telescope_coords.fov_lon.deg,
    telescope_coords.fov_lat.deg,
    alpha=0.4,
    color='gray'
)

print('Coords',telescope_coords.fov_lon.deg,telescope_coords.fov_lat.deg)

for i, name in enumerate(['Vela Pulsar']):
    star = SkyCoord.from_name(name)
    star_tel = star.transform_to(telescope_frame)
    print(star)

    plt.plot(star_tel.fov_lon.deg, star_tel.fov_lat.deg, '*')
    plt.annotate(
        s=name, xy=(star_tel.fov_lon.deg, star_tel.fov_lat.deg), xytext=(5, 5),
        textcoords='offset points', color=f'C{i}',
    )


plt.xlabel('fov_lon / {}'.format(telescope_coords.altaz.az.unit))
plt.ylabel('fov_lat / {}'.format(telescope_coords.altaz.alt.unit))
#plt.savefig('/home/spencers/fov.png')
plt.savefig('/home/spencers/geom.png',dpi=300)

#Combined plot

fig=plt.figure()

im1=customMap(model.skymap.data, 'FOV for %s' % (model.observer_source.date),fovval)

plt.axis('equal')
plt.scatter(
    telescope_coords.fov_lon.deg,
    telescope_coords.fov_lat.deg,
    alpha=0.4,
    color='gray'
)

print('Coords',telescope_coords.fov_lon.deg,telescope_coords.fov_lat.deg)

for i, name in enumerate(['Vela Pulsar']):
    star = SkyCoord.from_name(name)
    star_tel = star.transform_to(telescope_frame)
    print(star)

    plt.plot(star_tel.fov_lon.deg, star_tel.fov_lat.deg, '*')
    plt.annotate(
        s=name, xy=(star_tel.fov_lon.deg, star_tel.fov_lat.deg), xytext=(5, 5),
        textcoords='offset points', color='white',
    )


plt.xlabel('fov_lon / {}'.format(telescope_coords.altaz.az.unit))
plt.ylabel('fov_lat / {}'.format(telescope_coords.altaz.alt.unit))
plt.savefig('/home/spencers/combined.png',dpi=300)
