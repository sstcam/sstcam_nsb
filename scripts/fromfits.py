import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
from matplotlib import cm

np.set_printoptions(threshold=2000)

filename='/home/spencers/NSB_of_20220402003000_.fits'
#hdu1=fits.open('/home/spencers/NSB_of_20220402003000_.fits')
#print(hdu1)
#print(hdu1.info())
#print(dir(hdu1))
#print(hdu1[0].header)
#hdu1.close()

fov=fits.getdata(filename,ext=0)
fig=plt.figure()
print(fov.shape)
vmin, vmax = np.nanmin(fov), np.minimum(4 * np.nanmedian(fov), np.nanmax(fov))
plt.imshow(fov, cmap=cm.viridis, vmin=vmin, vmax=vmax)
plt.colorbar(label='Brightness (nLb)')
plt.savefig('fovd.png',dpi=300)
