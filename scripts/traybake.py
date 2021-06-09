import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column, MaskedColumn
from matplotlib import cm
import matplotlib as mpl
from scipy.ndimage import rotate
inputtable=Table.read('/home/spencers/sstcam_nsb/results/etacardarkgauss3/fluxes_crota_0.csv')

print(inputtable)
posx=inputtable['xcenter'].data
posy=inputtable['ycenter'].data
sums=inputtable['aperture_sum'].data
windowsize_r = 8
windowsize_c = 8
mean=np.mean(sums)
sig=np.std(sums)
maxval=np.amax(sums)
minval=np.amin(sums)

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

sums=np.asarray(cam_squaremaker(sums))
sums=np.rot90(sums,3)
fig=plt.Figure()
plt.imshow(sums,cmap=cm.viridis,norm=mpl.colors.LogNorm())
plt.title('Engineering Camera Frame\n Mean='+str('%.1f'%mean)+', $\sigma$='+str('%.1f'%sig)+'\n Max='+str('%.1f'%maxval)+', Min='+str('%.1f'%minval))
plt.colorbar(label='Relative Brightness (nLb/Pixel)')
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
plt.tight_layout()
plt.savefig('nLb_TM.png')
print(tms.shape)
