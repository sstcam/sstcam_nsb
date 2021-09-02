import numpy as np
from tqdm import trange
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
import pandas as pd

fov = 8.8

r=pd.read_csv('/home/spencers/gaiacat8.txt')

for j in range(9):
    M_max = 4
    M_min = 5 + j * 0.5
    N = 5000
    ra_list = 360 * np.random.rand(N)
    dec_list = (180 / np.pi) * np.arcsin(2.0 * (np.random.rand(N) - 0.5))
    n_array = []
    f = open(f"gaia_results_M{M_max}_M{M_min}.csv", "w")
    t = trange(N, desc="Bar desc", leave=True)
    rp=r[r['Magnitude'] > M_max]
    rp=r[r['Magnitude'] < M_min ]
    for i in t:
        ra = ra_list[i]
        dec = dec_list[i]
        delta_dec = fov / 2.0
        delta_ra = (fov / 2.0) / np.cos(dec * np.pi / 180.0)
        ra1 = ra - delta_ra
        ra2 = ra + delta_ra
        dec1 = dec - delta_dec
        dec2 = dec + delta_dec
        mask1 = np.logical_and(rp["RA"] > ra1, rp["RA"] < ra2)
        mask2 = np.logical_and(rp["DEC"] > dec1, rp["DEC"] < dec2)
        mask = np.logical_and(mask1, mask2)
        n = len(rp[mask])
        n_array.append(n)
        f.write(f"{i}, {ra:.2f}, {dec:.2f}, {n}\n")
        t.set_description(f"{i} of {N} RA={ra:.1f}, DEC={dec:.1f}, N={n}")
        t.refresh()
    mn = np.min(n_array)
    mx = np.max(n_array)
    NBINS = int((mx - mn) / 2)
    plt.hist(n_array, NBINS, alpha=0.8,label=str(M_min))

plt.title('FOV: '+str(fov)+' Degrees')
plt.xlabel('Number of Sources per FOV')
plt.ylabel('Frequency')
plt.legend(title='Minimum Magnitude')
plt.savefig('richfig.png')

