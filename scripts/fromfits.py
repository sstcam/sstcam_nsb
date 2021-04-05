import numpy as np
import matplotlib as mpl
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

filename='/store/spencers/NSBmaps/NSB_of_20220322013000_.fits'

hdu1=fits.open(filename)
#print(hdu1)
#print(hdu1.info())
#print(dir(hdu1))

print(hdu1.info())
fov=hdu1[0].data
print(hdu1[0].header)
fov=np.nan_to_num(fov)
wcs_fits=WCS(hdu1[0].header)
print('WCS from fits file:',wcs_fits)
wcs_input_dict={'CTYPE1': 'RA---HPX','CUNIT1': 'deg','CDELT1': -0.001,'CRPIX1': 5000,'CRVAL1':161.26477294,'NAXIS1': 10000,'CTYPE2': 'DEC--HPX','CUNIT2': 'deg','CDELT2': 0.001,'CRPIX2': 5000,'CRVAL2':-59.68443085,'NAXIS2': 10000,'CROTA1':0,'CROTA2':0}

wcs_dict=WCS(wcs_input_dict)
print('wcs_dict',wcs_dict)
funit=u.dimensionless_unscaled
print('funit',funit)
hdu1.close()

fig=plt.figure()
print(fov.shape)
ax=plt.subplot(projection=wcs_dict,label='overlays')
vmin, vmax = np.nanmin(fov), np.minimum(4 * np.nanmedian(fov), np.nanmax(fov))
f1=ax.imshow(fov, cmap=cm.viridis, vmin=vmin, vmax=vmax,origin='lower')


ax.coords.grid(True, color='white', ls='solid')
ax.coords[0].set_axislabel('Right Ascension')
ax.coords[1].set_axislabel('Declination')
ax.coords[0].set_format_unit('deg')

overlay = ax.get_coords_overlay('galactic')
overlay.grid(color='white', ls='dotted')
overlay[0].set_axislabel('Galactic Longitude')
overlay[1].set_axislabel('Galactic Latitude')

fig.colorbar(f1,ax=ax,label='Brightness (nLb)',pad=0.2)
plt.savefig('fovd.png',dpi=300)

data=u.Quantity(fov,unit=funit)
print(data,np.shape(data))
#data=w.pixel_to_world(data)

location = EarthLocation.from_geodetic(-70.317876,-24.681546)
obstime = Time('2022-03-22T01:30')

vela = SkyCoord.from_name("Eta Carinae")
print(vela)
altaz = AltAz(location=location, obstime=obstime)

pointing=vela.transform_to(altaz)

cam = CameraGeometry.from_name('CHEC') #Note this is very old, needs to be updated for prod 5.
pix_x = cam.pix_x
pix_y = cam.pix_y

optics = OpticsDescription.from_name('SST-ASTRI') #As is this
print(optics)

focal_length = optics.equivalent_focal_length

camera_frame = CameraFrame(focal_length=focal_length,rotation=0*u.deg,telescope_pointing=pointing)

cam_coords = SkyCoord(pix_x,pix_y,frame=camera_frame)

telescope_frame = TelescopeFrame(
    telescope_pointing=pointing,
    obstime=pointing.obstime,
    location=pointing.location,
)
telescope_coords = cam_coords.transform_to(telescope_frame)
print(telescope_coords.fk5)
sky_coords=telescope_coords.transform_to(ICRS())
print(sky_coords)
#print('Telescope coordinates:',telescope_coords,dir(telescope_coords))
print('WCS used:',wcs_dict)

aper=SkyRectangularAperture(sky_coords,w=0.19*u.deg,h=0.19*u.deg) #ASTRI Values for 7mm silicon
#aper=SkyCircularAperture(sky_coords,0.1*u.arcsec)
print(data,np.shape(data))
phot_table=aperture_photometry(data,aper,wcs=wcs_dict)

print(phot_table)
phot_table.write('fluxes.csv',overwrite=True)
fig=plt.figure()

sums=phot_table['aperture_sum']
sums=np.nan_to_num(sums)
print(sums)
plt.hist(sums.value)
plt.semilogy()
plt.xlabel('Pixel')
plt.ylabel('Integrated Sky Brightness (nLb / pixel)')
plt.savefig('hist.png')

fig=plt.figure()

plt.axis('equal')
plt.scatter(
    telescope_coords.fov_lon.deg,
    telescope_coords.fov_lat.deg,
    c=sums.value,
    norm=mpl.colors.LogNorm(),
    cmap=cm.viridis
)

plt.colorbar(label='Integrated brightness (nLb/pixel)')
plt.xlabel('fov_lon / {}'.format(telescope_coords.altaz.az.unit))
plt.ylabel('fov_lat / {}'.format(telescope_coords.altaz.alt.unit))
plt.savefig('integrated.png',dpi=300)
