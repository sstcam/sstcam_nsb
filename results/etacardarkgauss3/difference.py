import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column, MaskedColumn
from matplotlib import cm
import matplotlib as mpl
from scipy.ndimage import rotate

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


inputtable1=Table.read('/home/spencers/sstcam_nsb/results/etacardarkminus30mins/fluxes_crota_0.csv')
inputtable2=Table.read('/home/spencers/sstcam_nsb/results/etacardarkgauss3/fluxes_crota_0.csv')

posx1=inputtable1['xcenter'].data
posy1=inputtable1['ycenter'].data
sums1=inputtable1['aperture_sum_nLb'].data
sHz1=inputtable1['aperture_sum_Hz'].data/10**6 #Values in MHz                                                                                                                                              


posx2=inputtable2['xcenter'].data
posy2=inputtable2['ycenter'].data
sums2=inputtable2['aperture_sum_nLb'].data
sHz2=inputtable2['aperture_sum_Hz'].data/10**6 #Values in MHz   

windowsize_r = 8
windowsize_c = 8
wsr2 = 2
wsc2 = 2

diff_hz=sHz2-sHz1
print(diff_hz,np.shape(diff_hz))
mean=np.mean(diff_hz)
meanhz=np.mean(diff_hz)
sig=np.std(diff_hz)
sighz=np.std(diff_hz)
maxval=np.amax(diff_hz)
minval=np.amin(diff_hz)
maxvalhz=np.amax(diff_hz)
minvalhz=np.amin(diff_hz)

diff_hz=np.asarray(cam_squaremaker(diff_hz))
diff_hz=np.rot90(diff_hz,3)
fig=plt.Figure()
plt.imshow(diff_hz,cmap=cm.viridis,norm=mpl.colors.LogNorm())
plt.title('Difference Over 30 Minutes \n Engineering Camera Frame\n Mean='+str('%.1f'%mean)+' MHz, $\sigma$='+str('%.1f'%sig)+' MHz \n Max='+str('%.1f'%maxval)+' MHz , Min='+str('%.1f'%minval)+' MHz')
plt.colorbar(label='Change In Relative Brightness (Hz/Pixel)')
plt.xlabel('x (Pixels)')
plt.ylabel('y (Pixels)')
plt.tight_layout()
plt.savefig('diff_Hz_pixel.png')

tmshz=[]

# Crop out the window and calculate the histogram                                                                                                                                                          \
                                                                                                                                                                                                            
for r in range(0,diff_hz.shape[0], wsr2):
    for c in range(0,diff_hz.shape[1], wsc2):
        window = diff_hz[r:r+wsr2,c:c+wsc2]
        print(window)
        tmshz.append(np.mean(window))

tmshz=np.asarray(tmshz)
print(np.shape(tmshz))
tmshz=tmshz.reshape((24,24))
#flathz=tmshz[np.where(tmshz>0)]
flathz=tmshz
m2hz=np.mean(flathz)
s2hz=np.std(flathz)
maxval2hz=np.amax(flathz)
minval2hz=np.amin(flathz)

plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(tmshz)
plt.title('Difference Over 30 Minutes \n Engineering Camera Frame\n Mean='+str('%.1f'%m2hz)+' MHz, $\sigma$='+str('%.1f'%s2hz)+' MHz\n Max='+str('%.1f'%maxval2hz)+' MHz, Min='+str('%.1f'%minval2hz)+' MHz')
plt.colorbar(label='Change in Mean Relative Brightness (MHz/Superpixel)')
plt.xlabel('x (Superpixels)')
plt.ylabel('y (Superpixels)')
plt.tight_layout()
plt.savefig('diff_Hz_Superpixel.png')


tms=[]

# Crop out the window and calculate the histogram                                                                                                                                                           
for r in range(0,diff_hz.shape[0], windowsize_r):
    for c in range(0,diff_hz.shape[1], windowsize_c):
        window = diff_hz[r:r+windowsize_r,c:c+windowsize_c]
        print(window)
        tms.append(np.mean(window))

tms=np.asarray(tms)
tms=tms.reshape((6,6))
#flat=tms[np.where(tms>0)]
flat=tms
m2=np.mean(flat)
s2=np.std(flat)
maxval2=np.amax(flat)
minval2=np.amin(flat)

plt.clf()
plt.cla()
fig=plt.Figure()
plt.imshow(tms)
plt.title('Difference Over 30 Minutes \n Engineering Camera Frame\n Mean='+str('%.1f'%m2)+' MHz, $\sigma$='+str('%.1f'%s2)+' MHz \n Max='+str('%.1f'%maxval2)+' MHz , Min='+str('%.1f'%minval2)+' MHz')
plt.colorbar(label='Change in Mean Relative Brightness (Hz/TM)')
plt.xlabel('x (TMs)')
plt.ylabel('y (TMs)')
plt.tight_layout()
plt.savefig('diff_Hz_TM.png')

