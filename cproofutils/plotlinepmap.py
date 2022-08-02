
from pyglider.slocum import parse_logfiles
import pyglider.utils as pgutils

import datetime
import glob
import io
import numpy as np
import seawater
import xarray as xr

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def get_line_info(linename):

    if linename == 'linep':
        inst = """#Tofino -12602.988  4856.826
        #P4 -12640.0000   4839.0000
        #P5 -12710.0000   4841.5000
        #P6 -12740.0000   4844.6000
        #P7 -12810.0000   4846.6000
        #P8 -12840.0000   4849.0000
        #P9 -12910.0000   4851.4000
        #P10 -12940.0000   4853.6000
        #P11 -13010.0000   4856.0000
        #P12 -13040.0000   4858.2000
        #P13 -13140.0000   4902.6000
        #P14 -13240.0000   4907.4000
        #P15 -13340.0000   4912.0000
        #P16 -13440.0000   4917.0000
        #P17 -13540.0000   4921.0000
        #P18 -13640.0000   4926.0000
        #P19 -13740.0000   4930.0000
        #P20 -13840.0000   4934.0000
        #P21 -13940.0000   4938.0000
        #P22 -14040.0000   4942.0000
        #P23 -14140.0000   4946.0000
        #P24 -14240.0000   4950.2000
        #P25 -14336.3000   5000.0000
        #P35 -14418.2000   5000.0000
        #P26 -14500.0000   5000.0000"""
    elif linename == 'LCshort':
        inst = """#Tofino -12602.988  4856.826
        #LC00 -12640.0000 4815.01
        #LC01 -12735.9200 4750.9500"""


    f = io.StringIO(inst)
    lineP = {}
    lineP['wps'] = None
    lineP['nm']  = []
    for line in f:
        st = line.strip().split(' ')
        print(st)
        lon = pgutils.nmea2deg(float(st[1]))
        lat = pgutils.nmea2deg(float(st[-1]))
        if lineP['wps'] is None:
            lineP['wps'] = np.array([[lat, lon]])
        else:
            lineP['wps'] = np.append(lineP['wps'], np.array([[lat, lon]]), axis=0)
        lineP['nm'] += [st[0][1:]]

    # get distances
    dist = np.zeros(len(lineP['wps']))
    dist, ang  = seawater.extras.dist(lineP['wps'][:, 0], lineP['wps'][:, 1])
    dist = np.cumsum(np.append([0], dist))
    dist = dist - dist[1]
    lineP['dist'] = -dist  # make P4 = 0, and negative to west
    return lineP


def plotLinePMissionMap(figdir='./figs/', linename='linep',
                        outname='LinePMissionMap.png', dpi=200, logdir='./logs',
                        lonlim=[-145.2, -124.4], latlim=[46.1, 52.2],
                        topofile = '~cproof/Sites/gliderdata/smith_sandwell_topo_v8_2.nc'):

    utcnow = datetime.datetime.utcnow()

    lineP = get_line_info(linename=linename)

    # get positions, times, and ampH from logfiles:
    fns = glob.glob(f'{logdir}/*.log')
    fns.sort()
    print(fns)
    glider = parse_logfiles(fns)
    glider = glider.dropna(dim='surfacing')
    # get a distance along line.  Simple interp in lon which is prob OK here
    glider['LinePdist'] = ('surfacing', np.interp(glider['lon'],
                            lineP['wps'][::-1, -1], lineP['dist'][::-1]))

    print(glider.time[0].values)
    print(glider.time[-1].values)
    # make a 6-h interp of this...
    dt = np.timedelta64(6*3600, 's')
    time = np.arange(glider.time[-1].values, glider.time[0].values, -dt)[::-1]

    glider6h = xr.Dataset(coords={'time': time,
                                  'timemid': time[:-1] + dt / 2})

    for todo in ['lon', 'lat', 'ampH']:
        glider6h[todo] = ('time', np.interp(glider6h.time, glider.time, glider[todo]))
    dist, ang = seawater.extras.dist(glider6h['lat'], glider6h['lon'])
    glider6h['dist'] = ('timemid', dist)
    glider6h['head'] = ('timemid', np.mod(ang-270, 360))
    glider6h['speed'] = glider6h['dist'] / dt.astype(float) * 24 * 3600  # km/d
    glider6h['ampHperday'] = ('timemid', glider6h.ampH.diff(dim='time').values / dt.astype(float) * 24 * 3600)
    glider6h['distLineP'] = ('time', np.interp(glider6h['lon'], lineP['wps'][::-1, -1], lineP['dist'][::-1]))
    glider6h['speedLineP'] = ('timemid', glider6h.distLineP.diff(dim='time').values / dt.astype(float) * 24 * 3600)

    inds = min(5, len(glider6h.speedLineP))
    latestLinePspeed = np.mean(glider6h.speedLineP[-inds:].values)
    latestAmpHperday = np.mean(glider6h.ampHperday[-inds:].values)
    latestAmpHper100km = np.abs(latestAmpHperday / latestLinePspeed * 100)

    # string with current time:
    timest = f'{glider6h.time.values[-1]}'[:16]
    timest = timest.replace('T', ' ')

    # time since last:
    print(glider6h.time.values[-1].astype('datetime64[s]').tolist())
    sincelast = utcnow - (glider6h.time.values[-1].astype('datetime64[s]')).tolist()
    sincelast = str(sincelast)[:-10]

    # figure out how far we have to go:
    total = (lineP['dist'][0] - lineP['dist'][-1]) * 2
    dist_to_lineP = glider6h.distLineP[-1] - lineP['dist'][-1]
    print(total, dist_to_lineP.values)
    if glider6h.lon.min() < -144.5 and glider6h.speedLineP[-6] > 0:
        # we are returning
        dist_to_go = lineP['dist'][0] - glider6h.distLineP[-1]
    else:
        dist_to_go = total / 2 + dist_to_lineP

    # arrival time:
    try:
        predTime = glider6h.time[-1].values + np.timedelta64(int(dist_to_go / np.abs(latestLinePspeed) * 24 * 3600), 's')
        predTime = f'{predTime}'[:10]
    except (OverflowError):
        predTime ='Bad data'

    # Plots:

    fig, axs = plt.subplots(4, 1, constrained_layout=True, figsize=(8, 10),
                            gridspec_kw={'height_ratios':[1, 1, 1, 2]})

    axs[0].plot(glider6h.timemid, glider6h.speed, label='Absolute Speed')
    axs[0].plot(glider6h.timemid, -glider6h.speed, color='C0')
    axs[0].plot(glider6h.timemid, glider6h.speedLineP, label='Speed along Line')
    axs[0].plot(glider6h.timemid[-inds:], latestLinePspeed * np.ones(inds), color='C1', alpha=0.5, lw=4)
    text = f'Latest speed (30 h): {latestLinePspeed:1.1f} km/d\n'
    text += f'Predicted return time: {predTime}'
    axs[0].text(0.1, 0.5, text,
                transform=axs[0].transAxes)
    axs[0].grid('True')
    axs[0].legend()
    axs[0].set_ylabel('Speed [km/d]')
    axs[0].set_ylim([-40, 40])
    axs[0].set_title(f'{timest}; Processed: {str(utcnow)[:16]}\n{sincelast} since last update')

    axs[1].plot(glider6h.timemid, glider6h.ampHperday)
    axs[1].plot(glider6h.timemid[-inds:], latestAmpHperday * np.ones(inds), color='C0', alpha=0.5, lw=4)
    axs[1].grid('True')
    axs[1].set_ylabel('AmpH/d')
    axs[1].set_ylim(bottom=0, top=10)
    text =  f'Latest [30 h]:  {latestAmpHperday:4.1f} AmpH/d\n'
    text += f'Used [500 Amph] {glider6h.ampH.values[-1]:4.1f} Amph\n'
    text += f'Percent used:   {glider6h.ampH.values[-1] / 5:4.1f}%'
    axs[1].text(0.05, 0.1, text, family='monospace', fontsize='medium',
                        transform=axs[1].transAxes)

    axs[2].plot(glider6h.timemid, 100 * glider6h.ampHperday/glider6h.speed)
    axs[2].plot(glider6h.timemid, np.abs(100 * glider6h.ampHperday/glider6h.speedLineP))
    axs[2].plot(glider6h.timemid[-inds:], latestAmpHper100km * np.ones(inds), color='C1', alpha=0.5, lw=4)
    axs[2].set_ylim([0, 40])
    axs[2].grid('True')
    axs[2].set_ylabel('AmpH / 100 km')

    predAH = glider6h.ampH[-1].values + latestAmpHper100km*dist_to_go.values / 100
    print(dist_to_go / np.abs(latestLinePspeed))
    text = f'Latest [30 h]:  {latestAmpHper100km:4.1f} AmpH/100km\n'
    text += f'Predicted total AmpH: {predAH:1.0f} Amph'
    axs[2].text(0.05, 0.1, text, family='monospace', fontsize='medium',
                            transform=axs[2].transAxes)

    # plot map:
    ax = axs[3]
    with xr.open_dataset(topofile) as ds:
        # longitude, latitiude, ROSE
        ds = ds.sel(longitude=slice(-145.2 + 360, -124+360))
        ind = np.where((ds.latitude < 52.2) & (ds.latitude > 46.0))[0]
        ds = ds.sel(latitude=slice(46, 52.5))
        ds['longitude'] = ds.longitude - 360

        ax.contourf(ds.longitude, ds.latitude, -ds.ROSE, linewidths=0, cmap='Blues',
                levels=[0, 100, 200, 500, 1000, 2000, 3000, 4000, 5000], vmax=1000, vmin=-1000, alpha=0.9, zorder=-100)
    ax.set_ylim(latlim)
    ax.set_xlim(lonlim)
    ax.set_aspect(1/np.cos(49 * np.pi / 180))

    td = list(range(1, len(lineP['wps']), 4))
    td = [0] + td + [-1]
    for wpind in td:
        ax.plot(lineP['wps'][wpind, 1], lineP['wps'][wpind, 0], 'd', color='C1')

        ax.text(lineP['wps'][wpind, 1],
                lineP['wps'][wpind, 0],
                f" {lineP['nm'][wpind]} {lineP['dist'][wpind]:1.0f} km",
                color='C1', rotation=60, fontsize='smaller')
    ax.plot(glider.lon, glider.lat, 'm.', markersize=1)
    ax.plot(glider.lon[-1], glider.lat[-1], 'go', markersize=6)

    text = f'Distance travelled: {(total - dist_to_go.values):1.0f} km\n'
    text += f'Distance to go: {dist_to_go.values:1.0f} km'
    ax.text(0.05, 0.1, text, family='monospace', fontsize='medium',
                            transform=ax.transAxes, bbox={'facecolor':[1, 1, 1]})

    fig.savefig(f'{figdir}/{outname}', dpi=dpi)
