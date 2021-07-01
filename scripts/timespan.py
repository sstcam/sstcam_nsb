from nsb import config
from nsb.model import nsbModel
from nsb.mypycat import mypycat
from nsb.gaia import Gaia
from nsb.nsbtools import plotTimespan
import ephem
import matplotlib.pyplot as plt
import numpy as np

runname='timespan_Hz_monthlong_etacar'
con = config.TheConfiguration()
con.readConfig('/home/spencers/sstcam_nsb/configs/etacardarkgauss3.cfg')
t1=ephem.Date("2022/02/07 04:54:00")
t2= t1 + 30.0

PDEaverage=0.4 #Average PDE                                                          
sapixel= 8*10**-6 #Solid angle of pixel in steradians                                
mirrorarea = 7.3 #Mirror area in m^2                                                 
transmission = 0.85 #Telescope transmission                                          

mpc=mypycat()
source=mpc.get("Eta Carinae")

def plotTimespan_Hz(model,PDEaverage,sapixel,mirrorarea,transmission,norm=75.96*10**6):
    # verbose timing output
    if model.verbose:
        all_times = model.sunset + model.sunrise + model.sourceset + model.sourcerise + model.moonset + model.moonrise
        all_times = sorted(all_times, key=itemgetter(0))
        sun_down = False
        for i in range(len(all_times)):
            if "Sunset" in all_times[i][1]:
                sun_down = True
            if sun_down:
                if i > 0 and not ("Sunset" in all_times[i][1]):
                    print("\t\t\t%.2f hours" % ((all_times[i][0] - all_times[i - 1][0]) * 24))
                print(all_times[i][1])
            if "Sunrise" in all_times[i][1]:
                sun_down = False
                print()

    plt.rcParams['figure.figsize'] = 16, 12
    # fig, (sub1, sub2, sub3) = plt.subplots(3, 1, sharex=True)
    fig, (sub1, sub2) = plt.subplots(2, 1, sharex=True)
    plt.subplots_adjust(bottom=0.5)

    # source alt colorcoded
    co = []
    for value in model.sourcealt:
        co.append(np.interp(value, [0.0, 90.0], [0.0, 1.0]))

    if len(model.timestamps) <= 50:
        labels = ['      %s' % ephem.Date(t + model.time_start) for t in model.timestamps]
        plt.xticks(model.timestamps, labels, rotation=90)
    elif len(model.timestamps) > 50 and (model.time_end - model.time_start) < 2:
        ticks = []
        labels = []
        t = 0
        while t < (model.time_end - model.time_start):
            t += 1/24
            ticks.append(t)
            labels.append("      " + str(ephem.Date(t + model.time_start)))
        plt.xticks(ticks, labels, rotation=90)

    elif len(model.timestamps) > 50 and (model.time_end - model.time_start) < 15:
        ticks = []
        labels = []
        t = 0
        while t < (model.time_end - model.time_start):
            t += 0.25
            ticks.append(t)
            labels.append("      " + str(ephem.Date(t + model.time_start)))
        plt.xticks(ticks, labels, rotation=90)

    elif (model.time_end - model.time_start) < 50:
        ticks = []
        labels = []
        t = 0
        while t < (model.time_end - model.time_start):
            t += 1
            ticks.append(t)
            labels.append("      " + str(ephem.Date(t + model.time_start)))
        plt.xticks(ticks, labels, rotation=90)

    else:
        plt.xlabel('Days since %s' % model.time_start)

    sub1.minorticks_on()
    sub2.minorticks_on()
    bv=np.asarray(model.bright)
    sv = bv/10**9
    div = (1*10**4)*(1.0/np.pi)*((680.0-330.0)/505.0)
    div = 1.0/div
    sv = sv/div
    hv = sv*250 #NSB range in nm
    bv = hv*PDEaverage*sapixel*mirrorarea*transmission*10**9
    bv = bv * norm/np.nanmean(bv)
    bv = bv/10**6 #Value in MHz
    sub1.set_title('Estimated Sky-Brightness for Position of Source '+str(model.source.name)+'\n Mean Normalised to: '+str('%.1f'%(norm/10**6))+' Mhz,'+' Mean Dark Observation Rate'+'\n Mean = '+str('%.1f'%np.nanmean(bv))+' MHz, Median = '+str('%.1f'%np.nanmedian(bv))+' MHz \n $\sigma$ = '+str('%.1f'%np.nanstd(bv))+ ' MHz, Min = '+str('%.1f'%np.nanmin(bv))+' MHz, Max = '+str('%.1f'%np.nanmax(bv))+' MHz')
    # sub1.scatter(model.timestamps, model.bright, s=1, c=co, label='brightness')
    #+'\n Normalised to: '+str('%.1f'%(norm/10**6))+' Mhz,'
    sub1.scatter(model.timestamps, bv, s=1, label='Brightness (MHz)')
    #sub1.set_ylim((10, 10000))
    #sub1.set_yscale("log", nonposy='clip')
    sub1.set_xlim(left=0, right=model.time_end - model.time_start)
    #sub1.axhline(y=120, alpha=0.7, color='green', ls='--', linewidth=1)
    #sub1.axhline(y=420, alpha=0.7, color='orange', ls='--', linewidth=1)
    #sub1.axhline(y=550, alpha=0.7, color='red', ls='--', linewidth=1)
    sub1.set_ylabel('Brightness (MHz)')

    sub2.plot(model.timestamps, model.moonphase, label='Moon Phase')
    sub2.plot(model.timestamps, model.separation, label='Separation')
    sub2.plot(model.timestamps, model.moonalt, label='Moon Alt')
    sub2.plot(model.timestamps, model.sourcealt, label='Source Alt')
    sub2.set_xlim(left=0, right=model.time_end - model.time_start)
    sub2.axhspan(80, 100, facecolor='grey', alpha=0.2)
    sub2.axhline(y=0, linewidth=0.2, color='black')
    sub2.set_ylabel('Alt, Sep [deg]\nPhase [%]')

    # put the orange and grey bands for dayligt and moonlight
    for i in range(0, len(model.sunrise)):
        sub1.axvspan(model.sunrise[i][0], model.sunset[i][0], facecolor='orange', alpha=0.1)
        sub2.axvspan(model.sunrise[i][0], model.sunset[i][0], facecolor='orange', alpha=0.1)

    for i in range(0, len(model.moonrise)):
        sub1.axvspan(model.moonrise[i][0], model.moonset[i][0], facecolor='grey', alpha=0.1)
        sub2.axvspan(model.moonrise[i][0], model.moonset[i][0], facecolor='grey', alpha=0.1)

    sub1.legend(loc='upper right')
    sub2.legend(loc='upper right')

    # Pad margins so that markers don't get clipped by the axes
    # plt.margins(0.2)
    fig.tight_layout()

    plt.draw()
    return bv

gaiamap=Gaia(level=11)
model = nsbModel(con, gaiamap, t1, t2, version="hess_basic", threshold=400, timeresolution=15, verbose=False)
model.setSource(source=source)
model.calculateTimespan()
bv=plotTimespan_Hz(model,PDEaverage,sapixel,mirrorarea,transmission)
plt.savefig(str(runname)+'.png')
np.save(str(runname)+'_brightness.npy',bv)
