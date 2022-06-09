import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'large',
                   'axes.labelsize': 'large',
                   'axes.titlesize':'x-large'}

pylab.rcParams.update(params)

nsb=np.genfromtxt('nsb.txt')

wavelength=nsb[:,0]
intensity=nsb[:,1]
plt.plot(wavelength,intensity)

plt.xlabel('$\mathrm{Wavelength\ (nm)}$')
plt.ylabel('$\mathrm{Intensity\ I_{\\nu}\ (photons\ cm^{-2} ns^{-1} sr^{-1})}$')
plt.title('NSB Spectrum From La Palma')
plt.xlim(200,950)
plt.savefig('bandeplot.png')
plt.show()
