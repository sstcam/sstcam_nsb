from nsb import config
from nsb.model import nsbModel
from nsb.mypycat import mypycat
from nsb.gaia import Gaia
from nsb.nsbtools import plotTimespan
import ephem
import matplotlib.pyplot as plt

con = config.TheConfiguration()
con.readConfig('/home/spencers/sstcam_nsb/configs/etacardarkgauss3.cfg')
t1=ephem.Date("2022/02/07 04:54:00")
t2= t1 + 30.0

mpc=mypycat()
source=mpc.get("Eta Carinae")

gaiamap=Gaia(level=11)
model = nsbModel(con, gaiamap, t1, t2, version="hess_basic", threshold=400, timeresolution=15, verbose=False)
model.setSource(source=source)
model.calculateTimespan()
plotTimespan(model)
plt.savefig('timespan_monthlong_etacar.png')
