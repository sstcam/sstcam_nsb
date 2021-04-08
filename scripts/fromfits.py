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
from ctapipe.visualization import ArrayDisplay, CameraDisplay
from ctapipe.instrument import CameraGeometry
from ctapipe.instrument.optics import OpticsDescription
from ctapipe.coordinates import (GroundFrame,TiltedGroundFrame,NominalFrame,TelescopeFrame,CameraFrame,EngineeringCameraFrame)

np.set_printoptions(threshold=2000)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16

filename='/store/spencers/NSBmaps/NSB_of_20220322015000_delta20mins.fits'

hdu1=fits.open(filename)
#print(hdu1)
#print(hdu1.info())
#print(dir(hdu1))

print(hdu1.info())
fov=hdu1[0].data
print(hdu1[0].header)
fov=np.nan_to_num(fov)
wcs_fits=WCS(hdu1[0].header)

funit=1e-5*u.cd*np.pi**(-1)*u.m**(-2) # Photutils doesn't support nanolamberts as a unit natively, so convert to candela/m^2, which it does support.

print('WCS from fits file:',wcs_fits)
'''                                                                                                                                                                                                         
Nsb doesn't write fits headers correctly, so you have to define a wcs input dictionary yourself here that depends on the fov size used, the number of pixels, and the RA/DEC of your source. Essentially C
RVAL1 and 2 are RA/DEC in degrees, CRPIX1 defines the centre of the field as the 5000th pixel (assuming you're using 10000x10000 pixels), CRDELT defines the number of degrees per fov map pixel. CTYPE def
ines the co-ordinate system used, there are variations that could be valid but using the healpix (HPX) versions seems to work well (variations in this are likely to be too small to matter for our purpose
s anyway).                                                                                                                                                                                                  
CDELT1 and CDELT2 for 10 degree fov: -0.001, 0.001                                                                                                                                                          
CDELT1 and CDELT2 for 12 degree fov: -0.0012, 0.0012                                                                                                                                                        
'''

wcs_input_dict={'CTYPE1': 'RA---HPX','CUNIT1': 'deg','CDELT1': -0.0012,'CRPIX1': 5000,'CRVAL1':161.26477294,'NAXIS1': 10000,'CTYPE2': 'DEC--HPX','CUNIT2': 'deg','CDELT2': 0.0012,'CRPIX2': 5000,'CRVAL2':-59.68443085,'NAXIS2': 10000,'CROTA1':0,'CROTA2':0}

wcs_dict=WCS(wcs_input_dict)
print('wcs_dict',wcs_dict)
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
obstime = Time('2022-03-22T01:50')

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

camera_frame = CameraFrame(focal_length=focal_length,rotation=0*u.deg,telescope_pointing=pointing) #Hacky solution to camera rotation problem

cam_coords = SkyCoord(pix_x,pix_y,frame=camera_frame)

telescope_frame = TelescopeFrame(
    telescope_pointing=pointing,
    obstime=pointing.obstime,
    location=pointing.location,
)
engineering_frame=EngineeringCameraFrame(n_mirrors=2,location=pointing.location,obstime=pointing.obstime,focal_length=focal_length,telescope_pointing=pointing)
telescope_coords = cam_coords.transform_to(engineering_frame)
tc2=cam_coords.transform_to(telescope_frame)

print(telescope_coords.fk5)
sky_coords=telescope_coords.transform_to(ICRS())
sc2=tc2.transform_to(ICRS())

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


engineering_cam=cam.transform_to(engineering_frame)

fig, axs = plt.subplots(1, 2, constrained_layout=True, figsize=(12, 6))
display_camera = CameraDisplay(
    engineering_cam, 
    ax=axs[0], 
    image=sums.value,
    title="EngineeringCameraFrame",                                                                                                                                                                        
    norm=mpl.colors.LogNorm(),
    cmap='viridis'
)

#display_camera.set_limits_minmax(vmin,vmax)
display_camera.add_colorbar(label='Brightness (nLb/pixel)')

display_engineering = CameraDisplay(
    engineering_cam.transform_to(camera_frame),
    ax=axs[1],
    image=sums.value,
    title="CameraFrame",                                                                                                                                                                                   
    norm=mpl.colors.LogNorm(),
    cmap='viridis'
)
#display_engineering.set_limits_minmax(vmin,vmax)
display_engineering.add_colorbar(label='Brightness (nLb/pixel)')

plt.savefig('integrated.png',dpi=300)


fig=plt.figure()
plt.axis('equal')

angle=-90.0
theta=np.pi*angle/180
rotation=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
pos=np.array([tc2.fov_lat.deg,tc2.fov_lon.deg])
pos=np.dot(rotation,pos)

plt.scatter(pos[0],pos[1],c=sums.value,norm=mpl.colors.LogNorm(),cmap=cm.viridis,marker='s') #Hacky solution to rotate sky field to be the same rotation as the engineering camera frame
plt.colorbar(label='Integrated brightness (nLb/pixel)')                                                                                                                                                    
plt.xlabel('fov_lon / {}'.format(tc2.altaz.az.unit))
plt.ylabel('fov_lat / {}'.format(tc2.altaz.alt.unit))                                                                                                                                               
plt.savefig('skyframe.png',dpi=300)

