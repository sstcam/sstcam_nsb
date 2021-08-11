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

def calc_field_rotation(Lat,Azs,Alts):
    K=15.04/3600 #Work in seconds
    Azs=Azs*np.pi/180.0
    Alts=Alts*np.pi/180.0
    Lat=Lat*np.pi/180.0 #Must be in rad for numpy
    rots=K*np.cos(Lat)*np.cos(Azs)/np.cos(Alts)
    return rots

if __name__=='__main__':
    print(hipparcos.URL)
    obstime = Time('2022-03-24T07:32:00.000')
    loc = EarthLocation.from_geodetic(lon=-70.317876,lat=-24.681546,height=2176.6*u.m) #SST1 Paranal Position                                                                                                                                
    maglimit=4
    plotstars=True
    Maxtrack=0.19/10
    Lat=loc.lat
    with load.open('hip_main.dat') as f:
        df = hipparcos.load_dataframe(f)
    df=df[df['magnitude']<=maglimit]

    ts=load.timescale()
    t=ts.from_astropy(obstime)
    planets = load('de421.bsp')
    earth = planets['earth']
    stars= Star.from_dataframe(df)
    astrometric = earth.at(t).observe(stars)
    ra, dec, distance = astrometric.radec()
    sc=SkyCoord(ra=ra.radians*u.rad,dec=dec.radians*u.rad)
    AltAz=AltAz(location=loc,obstime=obstime)
    sc=sc.transform_to(AltAz)
    alts=sc.alt.deg
    azs=sc.az.deg
    Azvals=np.linspace(20,89.2,1000)
    Altvals=np.linspace(0,360,1000)
    Azs,Alts=np.meshgrid(Azvals,Altvals)
    rots=calc_field_rotation(Lat,Azs,Alts)
    print(rots)
    vmin=-Maxtrack
    vmax=Maxtrack
    cmap = plt.get_cmap('PiYG')
    print(rots.min(),rots.max())
    levels = MaxNLocator(nbins=1000).tick_values(vmin, vmax)
    cf=plt.contourf(Azs,Alts,rots,levels=levels)
    if plotstars==True:
        plt.scatter(alts,azs,marker='+',color='orange')
    cs=plt.contour(cf,levels=[-Maxtrack,Maxtrack],colors='r',linestyles='dashed')
    plt.colorbar(cf,label='Degrees/Second')
    colors = ['red']
    lines = [Line2D([0], [0], color=c, linewidth=3, linestyle='dashed') for c in colors]
    if plotstars==True:
        lines.append(Line2D([],[],color='orange',linestyle='None',marker='+'))
        labels=['Calibration Limit','Bright Star']
    else:
        lables=['Calibration Limit']
    plt.legend(lines,labels,loc='upper right')
    plt.xlim(20,89.2)
    plt.ylim(0,360)
    plt.xlabel('Altitude (Degrees)')
    plt.ylabel('Azimuth (Degrees)')
    plt.title('UTC Time '+str(obstime.utc)+'\n Field Rotation With Stars Brighter Than Mag '+str(maglimit)+'\n Calibration Limit '+'%.4f'%Maxtrack+ ' Degrees/Second')
    plt.tight_layout()
    plt.savefig('rot.png')
