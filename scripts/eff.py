import matplotlib.pyplot as plt
import numpy as np
from sstcam_simulation.utils.efficiency import CameraEfficiency, NSB_FLUX_UNIT
from sstcam_simulation.utils.sipm import SiPMSpecification
from sstcam_simulation.data import get_data
from astropy import units as u
from astropy.table import Table

tab=Table.read('/home/spencers/sstcam_nsb/results/etacarmoonpoint53gauss3/fluxes_crota_0.csv')
print(tab)
print('NSB rate Rich calc:',np.mean(tab['aperture_sum_Hz'].data))
nsb_fluxes=tab['aperture_sum_phot'].data*1.0/(u.ns*u.sr*u.nm*u.m**2)
print(nsb_fluxes)
meanval=np.mean(nsb_fluxes)
print('Mean NSB',meanval)

prod4_cam_eff = CameraEfficiency.from_prod4()

nsb_fluxes=prod4_cam_eff._integrate_nsb(nsb_fluxes,300*u.nm,550*u.nm)

print('Scaled NSB rate based on input:',prod4_cam_eff.get_scaled_nsb_rate(nsb_fluxes).to(u.Hz))

print('Nominal NSB rate:',prod4_cam_eff.nominal_nsb_rate.to(u.Hz))

