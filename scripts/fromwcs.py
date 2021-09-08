import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy import constants as const
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
from astroquery.vizier import Vizier
from nsb.mypycat import mypycat

np.set_printoptions(threshold=2000)
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 16

# Set run options here

filename='/more/spencers/mwpan2_RGB.fits'

location = EarthLocation.from_geodetic(-70.317876,-24.681546,height=2176.6) #Paranal SST-1
#location = EarthLocation.from_geodetic(lon=14.974609, lat=37.693267,height=1750*u.m) #ASTRI
obstime = Time('2022-02-07T04:54:0')
#raval = 266.1836287564894*u.deg # Source right ascension
#decval = 54.49192701147445*u.deg # Source declination
#raval=266.15347499999996*u.deg
#decval=53.95461111111111*u.deg #Draco coords for astri obs
raval=161.3*u.deg
decval=59.7*u.deg

# Width and height of pixels on sky
pixelw=0.19*u.deg
pixelh=0.19*u.deg
angle=-90.0 #Angle to rotate NSB fits images to match Engineering Camera Frame, set to -90 for consistant results with SSTCAM sandbox.
a2=0 #Rotation angle for stars in the Engineering Camera Frame plots.
a3=90.0 #Rotation angle for stars in the Skyframe (camera on sky) plots.
crota=0 #Rotation angle fed to WCS dictionary CROTA2 variable.

plotstars=True # Whether or not to plot star position overlays
nstars=5 # Number of stars to plot if selected
searchradius=5*u.deg # Radius to search for stars in
UTCtimecorrection=-3 #UTC time correction in hours
dt=TimeDelta(3600*UTCtimecorrection,format='sec') #Convert to hours (not currently used)
print('Rotation angle', crota)
#mainsource = SkyCoord.from_name("ASTRI DRACO")
mainsource = SkyCoord(ICRS(ra=raval,dec=decval)) #ASTRO DRACO field position

#Values for NSB Hz calculation
PDEaverage=0.4 #Average PDE
sapixel= 8*10**-6 #Solid angle of pixel in steradians
mirrorarea = 7.3 #Mirror area in m^2
transmission = 0.85 #Telescope transmission

cam = CameraGeometry.from_name('CHEC') # Note this is very old, needs to be updated for prod 5.
optics = OpticsDescription.from_name('SST-ASTRI') # As is this  

# Output filenames
tablefilename = 'fluxes_'+'crota_'+str(crota)+'.csv'
histname = 'hist'+'crota_'+str(crota)+'.png'
skyframename = 'skyframe'+'crota_'+str(crota)+'.png'
camframename = 'integrated'+'crota_'+str(crota)+'.png'
fovname='fovd'+'crota_'+str(crota)+'_TAN.png'

modestr='nLb/pixel'

hdu1=fits.open(filename)
print(hdu1)
hdu1.verify('silentfix')
print('Fits header',hdu1[0].header)
print('Fits WCS',WCS(hdu1[0].header))
#print(hdu1.info())
#print(dir(hdu1))

fov=hdu1[0].data
fov=np.nan_to_num(fov)

funit=1e-5*u.cd*np.pi**(-1)*u.m**(-2) # Photutils doesn't support nanolamberts as a unit natively, so convert to candela/m^2, which it does support.

wcs_dict=WCS(hdu1[0].header,naxis=2)

print('wcs_dict',wcs_dict)
print('funit',funit)
hdu1.close()

data=u.Quantity(fov,unit=funit)
print(data,np.shape(data))
#data=wcs_dict.pixel_to_world(data)

vmin, vmax = np.nanmin(fov), np.maximum(4 * np.nanmedian(fov), np.nanmax(fov)) #replace fov with data
print('VVals',vmin,vmax)

print(mainsource)
o2=obstime+dt
print(o2)

altaz = AltAz(location=location, obstime=obstime)
altaz2=AltAz(location=location,obstime=o2)
pointing=mainsource.transform_to(altaz)
p2=mainsource.transform_to(altaz2)


pix_x = cam.pix_x
pix_y = cam.pix_y


print(optics)

#focal_length = optics.equivalent_focal_length
focal_length = 2.15191*u.m #ASTRI value from sstcam-sandbox

camera_frame = CameraFrame(focal_length=focal_length,telescope_pointing=p2)

cam_coords = SkyCoord(pix_x,pix_y,frame=camera_frame)

telescope_frame = TelescopeFrame(
    telescope_pointing=p2,
    obstime=o2,
    location=pointing.location,
)
telframe2 = TelescopeFrame(telescope_pointing=p2,obstime=o2,location=p2.location)
engineering_frame=EngineeringCameraFrame(n_mirrors=2,location=p2.location,obstime=p2.obstime,focal_length=focal_length,telescope_pointing=p2)

print(engineering_frame, engineering_frame.rotation)

ef2=EngineeringCameraFrame(n_mirrors=2,location=p2.location,obstime=p2.obstime,focal_length=focal_length,telescope_pointing=p2)
telescope_coords = cam_coords.transform_to(engineering_frame)
tc2=cam_coords.transform_to(telframe2)

print(telescope_coords.fk5)
sky_coords=telescope_coords.transform_to(ICRS())
sc2=tc2.transform_to(ICRS())
print(sky_coords)
#print('Telescope coordinates:',telescope_coords,dir(telescope_coords))
print('WCS used:',wcs_dict)
data=np.mean(data,axis=0)
aper=SkyRectangularAperture(sky_coords,w=pixelw,h=pixelh) #ASTRI Values for 7mm silicon
#aper=SkyCircularAperture(sky_coords,0.1*u.arcsec)
print(data,np.shape(data),aper,wcs_dict,)
phot_table=aperture_photometry(data[:,:],aper,wcs=wcs_dict)

print(phot_table)

fig=plt.figure()

sums=phot_table['aperture_sum'] #Values in nLb

phot_table.remove_column('aperture_sum')

phot_table['aperture_sum_nLb']=sums
sums=np.nan_to_num(sums)


print('Umodified sum values',sums.value)
sv=sums.value/10**9
div=(1*10**4)*(1.0/np.pi)*((680.0-330.0)/505.0)
div=1.0/div
sv=sv/div
hv=sv.copy()
sv=sv*u.photon/(u.ns*u.sr*u.nm*u.m**2)

print('Mean value photons',np.mean(sv))

phot_table['aperture_sum_phot']=sv

hv = hv * 250 #NSB range in nm

print('Total photons: ',hv)
hv = hv*PDEaverage*sapixel*mirrorarea*transmission*10**9
hv = hv * u.Hz

print('Mean value Hz',np.mean(hv))

phot_table['aperture_sum_Hz']=hv

print(phot_table)
phot_table.write(str(tablefilename),overwrite=True)

sums=sums.value # Continue calculation in nLb

plt.hist(sums)
plt.semilogy()
plt.xlabel('Pixel')
plt.ylabel('Integrated Sky Brightness ('+modestr+')')
plt.savefig(str(histname))


engineering_cam=cam.transform_to(engineering_frame)
ec2=cam.transform_to(ef2)
if plotstars==True:
    # Get stars
    vizier = Vizier(catalog='V/50',columns=['_RAJ2000', '_DEJ2000', 'Vmag', 'Name','GLON','GLAT'])
    t = vizier.query_region(p2, radius=searchradius, catalog='V/50')[0]
    t.sort("Vmag")
    t = t[:nstars]
    stars = SkyCoord(ra=t['_RAJ2000'], dec=t['_DEJ2000'], frame='icrs')
    print(t)
    #stars_cam = stars.transform_to(camera_frame)
    stars_eng = stars.transform_to(ef2)
    

fig, axs = plt.subplots(1, 2, constrained_layout=True, figsize=(12, 6))
display_camera = CameraDisplay(
    engineering_cam, 
    ax=axs[0], 
    image=sums,
    norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax),
    title="EngineeringCameraFrame",
    cmap='viridis'
)
#norm=mpl.colors.LogNorm(),

#display_camera.set_limits_minmax(vmin,vmax)

try:
    display_camera.add_colorbar(label='Brightness ('+modestr+')')
except Exception:
    print('Colorbar exception')

display_engineering = CameraDisplay(
    engineering_cam.transform_to(camera_frame),
    ax=axs[1],
    image=sums,
    norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax),
    title="CameraFrame",
    cmap='viridis'
)
#display_engineering.set_limits_minmax(vmin,vmax)
display_engineering.add_colorbar(label='Brightness ('+modestr+')',pad=0.2)

theta=np.pi*angle/180
rotation=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]) #Hacky rotation matrix solution to engineeringcameraframe problem. 

theta2=np.pi*a2/180
rot2=np.array([[np.cos(theta2),-np.sin(theta2)],[np.sin(theta2),np.cos(theta2)]]) #Hacky rotation matrix solution to engineeringcameraframe problem.          

if plotstars==True:                                            
    engpos=np.asarray([stars_eng.x.value,stars_eng.y.value])
    engpos=np.dot(rot2,engpos)

theta3=np.pi*a3/180
rot3=np.asarray([[np.cos(theta3),-np.sin(theta3)],[np.sin(theta3),np.cos(theta3)]])
if plotstars==True:
    axs[0].plot(-1*engpos[0],engpos[1], 'wo', mfc='none', ms=25, mew=2)
    #axs[1].plot(stars_cam.x.value, stars_cam.y.value,'wo',mfc='none',ms=25,mew=2)
plt.savefig(str(camframename),dpi=300)


fig=plt.figure()
plt.axis('equal')

pos=np.array([tc2.fov_lat.deg,tc2.fov_lon.deg])
pos=np.dot(rotation,pos)

if plotstars==True:
    stars_tframe=stars.transform_to(tc2)
    starpos=np.array([stars_tframe.fov_lat.deg,stars_tframe.fov_lon.deg])
    starpos=np.dot(rot3,starpos)
    plt.plot(-1*starpos[0],starpos[1], 'wo', mfc='none', ms=25, mew=2)

plt.scatter(pos[0],pos[1],c=sums,cmap=cm.viridis,norm=mpl.colors.LogNorm(),marker='s')
#norm=mpl.colors.LogNorm(),

try:
    plt.colorbar(label='Integrated brightness ('+modestr+')')                                                                                                                                                
except Exception:
    print('Colorbar exception')

plt.xlabel('fov_lon / {}'.format(tc2.altaz.az.unit))
plt.ylabel('fov_lat / {}'.format(tc2.altaz.alt.unit))                                                                                                                                               
plt.savefig(str(skyframename),dpi=300)

fig=plt.figure()
fov=np.mean(fov,axis=0)
ax=plt.subplot(projection=wcs_dict,label='overlays')
f1=ax.imshow(fov,cmap=cm.viridis,norm=mpl.colors.LogNorm(),origin='lower')


ax.coords.grid(True, color='white', ls='solid')
ax.coords[0].set_axislabel('Right Ascension')
ax.coords[1].set_axislabel('Declination')
ax.coords[0].set_format_unit('deg')

overlay = ax.get_coords_overlay('galactic')
overlay.grid(color='white', ls='dotted')
overlay[0].set_axislabel('Galactic Longitude')
overlay[1].set_axislabel('Galactic Latitude')

if plotstars==True:
    ax.scatter(t['_RAJ2000'],t['_DEJ2000'], transform=ax.get_transform('icrs'),s=300,edgecolor='white', facecolor='none')

fig.colorbar(f1,ax=ax,label='Brightness ('+modestr+')',pad=0.2)

plt.savefig(str(fovname),dpi=300)
