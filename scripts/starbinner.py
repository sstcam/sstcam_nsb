

import numpy as np
import astropy.units as u
from astropy.time import Time
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import MaxNLocator
from matplotlib.lines import Line2D
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astroquery.vizier import Vizier
mpl.use('Agg')
from skyfield.api import Star, load
from skyfield.data import hipparcos
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd

if __name__=='__main__':
    df=pd.read_csv('/home/spencers/gaiacat8.txt')
    highmaglimit=4
    lowmaglimit=8
    covlimit=10

    df=df[df['Magnitude'] <= lowmaglimit]
    df=df[df['Magnitude'] > highmaglimit]
    ra=df['RA']
    dec=df['DEC']
    #nbins=(np.linspace(0,360,18),np.linspace(-90,90,36)) #10x10 degree bins
    nbins=(np.linspace(0,360,25),np.linspace(-90,90,51)) #7x7 degree bins
    hist,yedges,xedges=np.histogram2d(ra,dec,nbins)
    mstars=np.mean(hist)
    lenhist=np.shape(hist)[0]*np.shape(hist)[1]
    covlocs=hist[hist>covlimit]
    skycov=(np.float(len(covlocs))/lenhist)*100
    extent = [0,360,-90,90]
    fig=plt.figure(figsize=(12,6))
    ax = plt.gca()
    plt.xlabel('RA (deg)')
    plt.ylabel('DEC (deg)')
    plt.title('Hipparcos Stars with a Magnitude Brighter or Equal to '+str(lowmaglimit)+' mag,\n But Dimmer Than '+str(highmaglimit)+' mag.\n Mean Stars per Bin: '+str(mstars)+'\n Sky Coverage Assuming '+str(covlimit)+' Stars Needed: '+str(skycov)+'%')
    im=ax.imshow(hist,extent=extent)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im,label='Number of Stars Present',cax=cax)
    plt.tight_layout()
    plt.savefig('starhist'+str(highmaglimit)+'_'+str(lowmaglimit)+'_'+str(covlimit)+'_7deg.png',dpi=300)
