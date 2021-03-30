import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.time import Time
import copy
from astropy.coordinates import SkyCoord,EarthLocation,AltAz,ICRS,Galactic
from photutils.aperture import SkyRectangularAperture
from astropy.io import fits
from astropy.wcs import WCS
import sys
from matplotlib import cm
from photutils.aperture import RectangularAperture, aperture_photometry, SkyCircularAperture
from ctapipe.io import EventSource
from ctapipe.calib import CameraCalibrator
from ctapipe.utils import get_dataset_path
from ctapipe.visualization import ArrayDisplay
from ctapipe.instrument import CameraGeometry
from ctapipe.instrument.optics import OpticsDescription
from ctapipe.coordinates import (GroundFrame,TiltedGroundFrame,NominalFrame,TelescopeFrame,CameraFrame,EngineeringCameraFrame)

np.set_printoptions(threshold=2000)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16

filename='/store/spencers/NSBmaps/NSB_of_20220402003000_.fits'

hdu1=fits.open(filename)
#print(hdu1)
#print(hdu1.info())
#print(dir(hdu1))
#print(hdu1[0].header)
print(hdu1.info())
fov=hdu1[0].data
fov=np.nan_to_num(fov)

wcs_input_dict=wcs_input_dict = {'CTYPE1': 'RA','CUNIT1': 'deg','CDELT1': -0.0005,'CRPIX1': 5000,'CRVAL1':128.83606354,'NAXIS1': 10000,'CTYPE2': 'DEC','CUNIT2': 'deg','CDELT2': 0.0005,'CRPIX2': 5000,'CRVAL2':-45.17643181,'NAXIS2': 10000}

wcs=WCS(wcs_input_dict)
funit=u.dimensionless_unscaled
print('funit',funit)
hdu1.close()

fig=plt.figure()
print(fov.shape)
plt.subplot(projection=wcs)
vmin, vmax = np.nanmin(fov), np.minimum(4 * np.nanmedian(fov), np.nanmax(fov))
plt.imshow(fov, cmap=cm.viridis, vmin=vmin, vmax=vmax,origin='lower')

plt.grid(color='white', ls='solid')
plt.xlabel('Galactic Longitude')
plt.ylabel('Galactic Latitude')
plt.colorbar(label='Brightness (nLb)')
plt.savefig('fovd.png',dpi=300)

data=u.Quantity(fov,unit=funit)
print(data,np.shape(data))
#data=w.pixel_to_world(data)

location = EarthLocation.from_geodetic(-70.317876,-24.681546)
obstime = Time('2022-04-02T00:30')

vela = SkyCoord.from_name("vela pulsar")
print(vela)
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
sky_coords=telescope_coords.transform_to(ICRS())
print(sky_coords)
#print('Telescope coordinates:',telescope_coords,dir(telescope_coords))
print(wcs)
aper=SkyRectangularAperture(sky_coords,w=0.19*u.deg,h=0.19*u.deg) #ASTRI Values for 7mm silicon
#aper=SkyCircularAperture(sky_coords,0.1*u.arcsec)
print(aper,dir(aper))
print(data,np.shape(data))
phot_table=aperture_photometry(data,aper,wcs=wcs)

print(phot_table)
phot_table.write('fluxes.csv')
fig=plt.figure()
np.set_printoptions(sys.maxsize)

sums=phot_table['aperture_sum']
print(sums)
plt.hist(sums)
plt.xlabel('Pixel')
plt.ylabel('Integrated Sky Brightness (nLb)')
plt.savefig('hist.png')
