import numpy as np
import astropy.units as u
from astropy.time import Time, TimeDelta
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
import healpy as hp
from astropy.coordinates import EarthLocation,SkyCoord,AltAz,Galactic, ICRS
from astropy.time import Time
import healpy as hp
import numpy as np
import matplotlib.gridspec as gridspec

def cat2hpx(lon, lat, nside, radec=True):
    """
    Convert a catalogue to a HEALPix map of number counts per resolution
    element.

    Parameters
    ----------
    lon, lat : (ndarray, ndarray)
        Coordinates of the sources in degree. If radec=True, assume input is in the icrs
        coordinate system. Otherwise assume input is glon, glat

    nside : int
        HEALPix nside of the target map

    radec : bool
        Switch between R.A./Dec and glon/glat as input coordinate system.

    Return
    ------
    hpx_map : ndarray
        HEALPix map of the catalogue number counts in Galactic coordinates

    """

    npix = hp.nside2npix(nside)

    if radec:
        eq = SkyCoord(lon, lat, frame='icrs', unit='deg')
        l, b = eq.galactic.l.value, eq.galactic.b.value
    else:
        l, b = lon, lat

    # conver to theta, phi
    theta = np.radians(90. - b)
    phi = np.radians(l)

    # convert to HEALPix indices
    indices = hp.ang2pix(nside, theta, phi)

    idx, counts = np.unique(indices, return_counts=True)

    # fill the fullsky map
    hpx_map = np.zeros(npix, dtype=int)
    hpx_map[idx] = counts

    return hpx_map

def pixarea2nside(area):
    """Closest HEALPix nside for a given area, shamelessly stolen from desitarget.
    Parameters
    ----------
    area : :class:`float`
        area in square degrees.
    Returns
    -------
    :class:`int`
        HEALPix nside that corresponds to passed area.
    Notes
    -----
        - Only considers 2**x nside values (1, 2, 4, 8 etc.)
    """
    # ADM loop through nsides until we cross the passed area.
    nside = 1
    while (hp.nside2pixarea(nside, degrees=True) > area):
        nside *= 2

    # ADM there is no nside of 0 so bail if nside is still 1.
    if nside > 1:
        # ADM is the nside with the area that is smaller or larger
        # ADM than the passed area "closest" to the passed area?
        smaller = hp.nside2pixarea(nside, degrees=True)
        larger = hp.nside2pixarea(nside//2, degrees=True)
        if larger/area < area/smaller:
            return nside//2

    return nside

def calcobserve():
    '''Random Pointing generator for Paranal'''
    noyears=20
    azs=np.random.randint(0,360,size=8760*noyears)*u.deg
    alts=np.random.randint(20,90,size=8760*noyears)*u.deg
    loc = EarthLocation.from_geodetic(lon=-70.317876,lat=-24.681546,height=2176.6*u.m) #SST1 Paranal Position
    dt=TimeDelta(3600,format='sec')
    obstimes=Time("2021-02-07T04:54:54.0")+np.random.randint(0,8760*noyears,size=8760*noyears)*dt
    skycoords=AltAz(az=azs,alt=alts,location=loc,obstime=obstimes)
    skycoords=skycoords.transform_to(ICRS())
    print(skycoords)
    return skycoords

if __name__=='__main__':
    NSIDE=8
    print('Appropriate NSIDE to use',pixarea2nside(52.524))
    df=pd.read_csv('/home/spencers/gaiacat8.txt')
    highmaglimit=4
    lowmaglimit=7
    covlimit=30
    fig=plt.figure(figsize=(10,6))
    randopoint=calcobserve()
    avmap=cat2hpx(randopoint.ra,randopoint.dec,NSIDE,radec=True)
    coverage=np.zeros(np.shape(avmap))
    covlocs=np.where(avmap>0)
    coverage[covlocs]=1
    print(avmap)
    im=hp.mollview(coverage,title='CTA-South Sky Coverage',hold=True,xsize=1200)
    plt.savefig('ctasouthskycoverage.png')
    plt.clf()
    plt.cla()
    df=df[df['Magnitude'] <= lowmaglimit]
    df=df[df['Magnitude'] > highmaglimit]
    ra=df['RA']
    dec=df['DEC']
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
    plt.suptitle('Gaia Stars with a Magnitude Brighter or Equal to '+str(lowmaglimit)+' mag,\n but Dimmer Than '+str(highmaglimit)+' mag.\n Bin Size: '+'%.3f'%hp.nside2pixarea(NSIDE)+' Steradians, Mean Stars per Bin: '+str('%.3f'%mstars)+'\n Sky Coverage Given '+str(covlimit)+' Stars Needed and CTA-South Observability: '+str('%.3f'%skycov)+'%')
    plt.axes(ax2)
    im=hp.mollview(hpx_map,title='',unit='Number of Stars',hold=True,xsize=1200)
    plt.savefig('starhist'+str(highmaglimit)+'_'+str(lowmaglimit)+'_'+str(covlimit)+'_CTAS.png',dpi=300)
