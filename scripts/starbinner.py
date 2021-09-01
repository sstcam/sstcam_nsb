

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


if __name__=='__main__':
    print(hipparcos.URL)
    obstime = Time('2022-03-24T07:32:00.000')
    loc = EarthLocation.from_geodetic(lon=-70.317876,lat=-24.681546,height=2176.6*u.m) #SST1 Paranal Position                                                                                                                                
    highmaglimit=4
    lowmaglimit=8
    covlimit=30
    Lat=loc.lat

    with load.open('hip_main.dat') as f:
        df = hipparcos.load_dataframe(f)

    df=df[df['magnitude'] <= lowmaglimit]
    df=df[df['magnitude'] > highmaglimit]
    ts=load.timescale()
    t=ts.from_astropy(obstime)
    planets = load('de421.bsp')
    earth = planets['earth']
    stars= Star.from_dataframe(df)
    astrometric = earth.at(t).observe(stars)
    ra, dec, distance = astrometric.radec()
    #nbins=(np.linspace(0,360,18),np.linspace(-90,90,36)) #10x10 degree bins
    nbins=(np.linspace(0,360,25),np.linspace(-90,90,51)) #7x7 degree bins
    hist,yedges,xedges=np.histogram2d(ra._degrees,dec.degrees,nbins)
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
