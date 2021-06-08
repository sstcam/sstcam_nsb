import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, Column, MaskedColumn
from matplotlib import cm
import matplotlib as mpl

inputtable=Table.read('fluxes_crota_0.csv')
print(inputtable)
posx=inputtable['xcenter'].data
posy=inputtable['ycenter'].data
sums=inputtable['aperture_sum'].data
windowsize_r = 8
windowsize_c = 8

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

sums=cam_squaremaker(sums)

plt.imshow(sums,cmap=cm.viridis,norm=mpl.colors.LogNorm())

plt.colorbar()
plt.show()

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
plt.imshow(tms)
plt.colorbar(label='Relative Brightness (nLb/TM)')
plt.show()
print(tms.shape)
