import numpy as np
import sys
import glob
import fnmatch
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
from starbinner import *

if __name__=='__main__':
    specdic={}
    print(hipparcos.URL)
    obstime = Time('2022-03-24T07:32:00.000')
    loc = EarthLocation.from_geodetic(lon=-70.317876,lat=-24.681546,height=2176.6*u.m) #SST1 Paranal Position                                                                                                                                
    highmaglimit=4
    lowmaglimit=7
    covlimit=30
    NSIDE=8
    Lat=loc.lat

    with load.open('hip_main.dat') as f:
        df = hipparcos.load_dataframe(f)
    print(df)
    df=df[df['magnitude'] <= lowmaglimit]
    df=df[df['magnitude'] > highmaglimit]
    ra=df['ra_degrees'].to_numpy()
    dec=df['dec_degrees'].to_numpy()
    spectype=df['SpType'].to_numpy()
    for myfile in sorted(glob.glob('/home/spencers/dat/*')):
        fstring=myfile[19:-4]
        data=np.genfromtxt(myfile)[:,:2]
        specdic[fstring]=np.asarray(data)
        
        
    colors=['*o*','*b*','*a*','*f*','*w*','*g*','*k*','*m*']
    colordict={}
        
    for color in colors:
        keys=fnmatch.filter(specdic,color)
        sums=np.zeros((len(specdic['uko5v']),2))
        for key in keys:
            sums=sums+specdic[key]
        sums=sums/len(keys)
        colordict[color]=sums

    for key in colordict.keys():
        plt.plot(colordict[key][:,0],colordict[key][:,1],label=key)

    plt.legend(title='Spectral Class')
    plt.xlabel('Wavelength (Angstroms)')
    plt.ylabel('Mean Flux')
    plt.savefig('spectrum.png',dpi=300)
    
    fig=plt.figure(figsize=(10,6))
    randopoint=calcobserve()
    avmap=cat2hpx(randopoint.ra,randopoint.dec,NSIDE,radec=True)
    coverage=np.zeros(np.shape(avmap))
    covlocs=np.where(avmap>0)
    coverage[covlocs]=1
    im=hp.mollview(coverage,title='CTA-South Sky Coverage',hold=True,xsize=1200,unit='Source Observability')
    plt.savefig('ctasouthskycoverage.png')
    plt.clf()
    plt.cla()
    goodra=np.where(np.logical_and(dec<90,dec>-90))
    ra=ra[goodra]
    dec=dec[goodra]
    spectype=spectype[goodra]
    hpx_map=cat2hpx(ra,dec,nside=NSIDE,radec=True)
    sumcov=np.sum(coverage)
    skycov=len(np.where(np.logical_and(hpx_map>covlimit,coverage>0))[0])/float(sumcov)*100
    hpx_map[np.where(coverage==0)]=0
    mstars=np.mean(hpx_map)
    gs = gridspec.GridSpec(2, 1,height_ratios=[1,200])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    #plt.axes(ax1)                                                                                                                                                                                                                                                                                                                                                          
    ax1.axis('off')
    plt.suptitle('Hiparcos Stars with a Magnitude Brighter or Equal to '+str(lowmaglimit)+' mag,\n but Dimmer Than '+str(highmaglimit)+' mag.\n Bin Size: '+'%.3f'%hp.nside2pixarea(NSIDE)+' Steradians, Mean Stars per Bin: '+str('%.3f'%mstars)+'\n Sky Coverage Given '+str(covlimit)+' Stars Needed and CTA-South Observability: '+str('%.3f'%skycov)+'%')
    plt.axes(ax2)
    im=hp.mollview(hpx_map,title='',unit='Number of Stars',hold=True,xsize=1200)
    plt.savefig('starhist'+str(highmaglimit)+'_'+str(lowmaglimit)+'_'+str(covlimit)+'_CTAS.png',dpi=300)



