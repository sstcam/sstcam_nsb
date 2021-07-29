import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column, MaskedColumn
from matplotlib import cm
import matplotlib as mpl
from scipy.ndimage import rotate
inputtable=Table.read('/home/spencers/sstcam_nsb/results/etacardarkgauss3/fluxes_crota_0.csv')

print(inputtable.keys())
posx=inputtable['xcenter'].data
posy=inputtable['ycenter'].data
sums=inputtable['aperture_sum_nLb'].data
sHz=inputtable['aperture_sum_Hz'].data/10**6 #Values in MHz
windowsize_r = 8
windowsize_c = 8
wsr2 = 2
wsc2 = 2
mean=np.mean(sums)
meanhz=np.mean(sHz)
sig=np.std(sums)
sighz=np.std(sHz)
maxval=np.amax(sums)
minval=np.amin(sums)
maxvalhz=np.amax(sHz)
minvalhz=np.amin(sHz)

def cam_squaremaker(data):
    '''Function to translate CHEC-S integrated images into square arrays for
    analysis purposes.'''
    square = np.zeros(2304)
    i = 0
    while i < 48:
        if i < 8:
            xmin = 48 * i + 8
            xmax = 48 * (i + 1) - 8
            square[xmin:xmax] = data[i * 32:(i + 1) * 32]
            i = i + 1
        elif i > 7 and i < 40:
            square[384:1920] = data[256:1792]
            i = 40
        else:
            xmin = 48 * i + 8
            xmax = 48 * (i + 1) - 8
            square[xmin:xmax] = data[512 + i * 32:544 + i * 32]
            i = i + 1

    square.resize((48, 48))
    square = np.flip(square, 0)
    return square

def im_maxdiff(image):
    '''Takes a 2D array and calculates the maximum differences between x and y pixels. Extremely stupid implementation.'''

    xshape=np.shape(image)[0]
    yshape=np.shape(image)[1]
    newim=np.zeros((xshape,yshape))

    for i in np.arange(xshape):
        for j in np.arange(yshape):
            if image[i,j]>1e-6:
                nearvals=[]
                try:
                    left=abs(image[i,j]-image[i-1,j])
                    nearvals.append(left)
                except IndexError:
                    pass
                try:
                    right=abs(image[i,j]-image[i+1,j])
                    nearvals.append(right)
                except IndexError:
                    pass
                try:
                    up=abs(image[i,j]-image[i,j+1])
                    nearvals.append(up)
                except IndexError:
                    pass
                try:
                    down=abs(image[i,j]-image[i,j-1])
                    nearvals.append(down)
                except IndexError:
                    pass
                newim[i,j]=max(nearvals)
            else:
                continue
    
    return newim
            


sums=np.asarray(cam_squaremaker(sums))
sums=np.rot90(sums,3)
fig=plt.Figure()
plt.imshow(sums,cmap=cm.viridis,norm=mpl.colors.LogNorm())
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%mean)+', $\sigma$='+str('%.1f'%sig)+'\n Max='+str('%.1f'%maxval)+', Min='+str('%.1f'%minval))
plt.colorbar(label='Relative Brightness (nLb/Pixel)')
plt.xlabel('x (Pixels)')
plt.ylabel('y (Pixels)')
plt.tight_layout()
plt.savefig('nLb_pixel.png')

tms=[]

# Crop out the window and calculate the histogram
for r in range(0,sums.shape[0], windowsize_r):
    for c in range(0,sums.shape[1], windowsize_c):
        window = sums[r:r+windowsize_r,c:c+windowsize_c]
        print(window)
        tms.append(np.mean(window))

tms=np.asarray(tms)
print(np.shape(tms))
tms=tms.reshape((6,6))
flat=tms[np.where(tms>0)]
m2=np.mean(flat)
s2=np.std(flat)
maxval2=np.amax(flat)
minval2=np.amin(flat)

plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(tms)
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%m2)+', $\sigma$='+str('%.1f'%s2)+'\n Max='+str('%.1f'%maxval2)+', Min='+str('%.1f'%minval2))
plt.colorbar(label='Mean Relative Brightness (nLb/TM)')
plt.xlabel('x (TMs)')
plt.ylabel('y (TMs)')
plt.tight_layout()
plt.savefig('nLb_TM.png')
print(tms.shape)

# Crop out the window and calculate the histogram                                                                                                                                                           
tm_superpixel=[]
for r in range(0,sums.shape[0], wsr2):
    for c in range(0,sums.shape[1], wsc2):
        window = sums[r:r+wsr2,c:c+wsc2]
        print(window)
        tm_superpixel.append(np.mean(window))

tms=np.asarray(tm_superpixel)
print(np.shape(tm_superpixel))
tms=tms.reshape((24,24))
flat=tms[np.where(tms>0)]
m2=np.mean(flat)
s2=np.std(flat)
maxval2=np.amax(flat)
minval2=np.amin(flat)

plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(tms)
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%m2)+', $\sigma$='+str('%.1f'%s2)+'\n Max='+str('%.1f'%maxval2)+', Min='+str('%.1f'%minval2))
plt.colorbar(label='Mean Relative Brightness (nLb/Superpixel)')
plt.xlabel('x (Superpixels)')
plt.ylabel('y (Superpixels)')
plt.tight_layout()
plt.savefig('nLb_Superpixel.png')
print(tms.shape)





sHz=np.asarray(cam_squaremaker(sHz))
sHz=np.rot90(sHz,3)
plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(sHz,cmap=cm.viridis,norm=mpl.colors.LogNorm())
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%meanhz)+' MHz, $\sigma$='+str('%.1f'%sighz)+' MHz\n Max='+str('%.1f'%maxvalhz)+' MHz, Min='+str('%.1f'%minvalhz)+' MHz')
plt.colorbar(label='Relative Brightness (MHz/Pixel)')
plt.xlabel('x (Pixels)')
plt.ylabel('y (Pixels)')
plt.tight_layout()
plt.savefig('Hz_pixel.png')


plt.clf()
plt.cla()
fig=plt.Figure()
diffim=im_maxdiff(sHz)
plt.imshow(diffim)
plt.title('Engineering Camera Frame\nMaximum Change in NSB Per Pixel\n Mean='+str('%.1f'%np.mean(diffim[np.where(diffim>0)]))+' MHz, $\sigma$='+str('%.1f'%np.std(diffim[np.where(diffim>0)]))+' MHz\n Max='+str('%.1f'%np.nanmax(diffim[np.where(diffim>0)]))+' MHz, Min='+str('%.1f'%np.nanmin(diffim[np.where(diffim>0)]))+' MHz')
plt.colorbar(label='Max Difference in Mean Relative Brightness (MHz/Pixel)')
plt.xlabel('x (Pixels)')
plt.ylabel('y (Pixels)')
plt.tight_layout()
plt.savefig('Hz_pixel_diff.png')


tmshz=[]

# Crop out the window and calculate the histogram
for r in range(0,sHz.shape[0], windowsize_r):
    for c in range(0,sHz.shape[1], windowsize_c):
        window = sHz[r:r+windowsize_r,c:c+windowsize_c]
        print(window)
        tmshz.append(np.mean(window))

tmshz=np.asarray(tmshz)
print(np.shape(tmshz))
tmshz=tmshz.reshape((6,6))
flathz=tmshz[np.where(tmshz>0)]
m2hz=np.mean(flathz)
s2hz=np.std(flathz)
maxval2hz=np.amax(flathz)
minval2hz=np.amin(flathz)

plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(tmshz)
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%m2hz)+' MHz, $\sigma$='+str('%.1f'%s2hz)+' MHz\n Max='+str('%.1f'%maxval2hz)+' MHz, Min='+str('%.1f'%minval2hz)+' MHz')
plt.colorbar(label='Mean Relative Brightness (MHz/TM)')
plt.xlabel('x (TMs)')
plt.ylabel('y (TMs)')
plt.tight_layout()
plt.savefig('Hz_TM.png')

plt.clf()
plt.cla()
fig=plt.Figure()
diffim=im_maxdiff(tmshz)
plt.imshow(diffim)
plt.title('Engineering Camera Frame\nMaximum Change in NSB Per TM\n Mean='+str('%.1f'%np.mean(diffim[np.where(diffim>0)]))+' MHz, $\sigma$='+str('%.1f'%np.std(diffim[np.where(diffim>0)]))+' MHz\n Max='+str('%.1f'%np.nanmax(diffim[np.where(diffim>0)]))+' MHz, Min='+str('%.1f'%np.nanmin(diffim[np.where(diffim>0)]))+' MHz')
plt.colorbar(label='Max Difference in Mean Relative Brightness (MHz/TM)')
plt.xlabel('x (TMs)')
plt.ylabel('y (TMs)')
plt.tight_layout()
plt.savefig('Hz_TM_diff.png')

tmshz=[]

# Crop out the window and calculate the histogram                                                                                                                                                           
for r in range(0,sHz.shape[0], wsr2):
    for c in range(0,sHz.shape[1], wsc2):
        window = sHz[r:r+wsr2,c:c+wsc2]
        print(window)
        tmshz.append(np.mean(window))

tmshz=np.asarray(tmshz)
print(np.shape(tmshz))
tmshz=tmshz.reshape((24,24))
flathz=tmshz[np.where(tmshz>0)]
m2hz=np.mean(flathz)
s2hz=np.std(flathz)
maxval2hz=np.amax(flathz)
minval2hz=np.amin(flathz)

plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(tmshz)
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%m2hz)+' MHz, $\sigma$='+str('%.1f'%s2hz)+' MHz\n Max='+str('%.1f'%maxval2hz)+' MHz, Min='+str('%.1f'%minval2hz)+' MHz')
plt.colorbar(label='Mean Relative Brightness (MHz/Superpixel)')
plt.xlabel('x (Superpixels)')
plt.ylabel('y (Superpixels)')
plt.tight_layout()
plt.savefig('Hz_Superpixel.png')

plt.clf()
plt.cla()
fig=plt.Figure()
diffim=im_maxdiff(tmshz)
plt.imshow(diffim)
plt.title('Engineering Camera Frame\nMaximum Change in NSB Per Superpixel\n Mean='+str('%.1f'%np.mean(diffim[np.where(diffim>0)]))+' MHz, $\sigma$='+str('%.1f'%np.std(diffim[np.where(diffim>0)]))+' MHz\n Max='+str('%.1f'%np.nanmax(diffim[np.where(diffim>0)]))+' MHz, Min='+str('%.1f'%np.nanmin(diffim[np.where(diffim>0)]))+' MHz')
plt.colorbar(label='Max Difference in Mean Relative Brightness (MHz/Superpixel)')

plt.xlabel('x (Superpixels)')
plt.ylabel('y (Superpixels)')
plt.tight_layout()
plt.savefig('Hz_Superpixel_diff.png')
