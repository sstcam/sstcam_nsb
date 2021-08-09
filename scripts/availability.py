import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D


mpl.use('Agg')

def calc_field_rotation(Lat,Azs,Alts):
    K=15.04/3600 #Work in seconds
    Azs=Azs*np.pi/180.0
    Alts=Alts*np.pi/180.0
    Lat=Lat*np.pi/180.0 #Must be in rad for numpy
    rots=K*np.cos(Lat)*np.cos(Azs)/np.cos(Alts)
    return rots

if __name__=='__main__':
    
    Lat=-24.681546
    Maxtrack=0.19/40 #Assume 0.19 degree pixel size and 40s calibration time
    print(Maxtrack)
    Azvals=np.linspace(0,90,1000)
    Altvals=np.linspace(0,360,1000)
    Azs,Alts=np.meshgrid(Azvals,Altvals)
    rots=calc_field_rotation(Lat,Azs,Alts)
    print(rots)
    fig=plt.figure()
    sigma=np.std(rots)
    vmin=-Maxtrack
    vmax=Maxtrack
    cmap = plt.get_cmap('PiYG')
    print(rots.min(),rots.max())
    levels = MaxNLocator(nbins=1000).tick_values(vmin, vmax)
    cf=plt.contourf(Azs,Alts,rots,levels=levels)
    cs=plt.contour(cf,levels=[-Maxtrack,Maxtrack],colors='r',linestyles='dashed')
    plt.colorbar(cf,label='Degrees/Second')
    colors = ['red']
    lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='dashed') for c in colors]
    labels=['Calibration Limit']
    plt.legend(lines,labels)
    plt.xlabel('Altitude (Degrees)')
    plt.ylabel('Azimuth (Degrees)')
    plt.title('Calibration Limit '+'%.4f'%Maxtrack+ ' Degrees/Second')
    plt.savefig('rot.png')
