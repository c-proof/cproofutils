
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

    if linename == 'calvert':
        inst = """#Koeye -12753.3  5145.5
        #QCS01 -12814.3000   5142.3300
        #MPA4 -12839.6300   5124.5600
        #XSM -12907.5200   5118.8100
        #Shelf -12951.3300   5104.8700
        #Turn -13303.578   4956.457"""

    f = io.StringIO(inst)
    Calvert = {}
    Calvert['wps'] = None
    Calvert['nm']  = []
    for line in f:
        st = line.strip().split(' ')
        print(st)
        lon = pgutils.nmea2deg(float(st[1]))
        lat = pgutils.nmea2deg(float(st[-1]))
        if Calvert['wps'] is None:
            Calvert['wps'] = np.array([[lat, lon]])
        else:
            Calvert['wps'] = np.append(Calvert['wps'], np.array([[lat, lon]]), axis=0)
        Calvert['nm'] += [st[0][1:]]

    # get distances
    dist = np.zeros(len(Calvert['wps']))
    dist, ang  = seawater.extras.dist(Calvert['wps'][:, 0], Calvert['wps'][:, 1])
    dist = np.cumsum(np.append([0], dist))
    dist = dist - dist[1]
    Calvert['dist'] = -dist  # make QCS01 = 0, and negative to west
    return Calvert


def plotCalvertMissionMap(figdir='./figs/', linename='calvert',
                        outname='CalvertMissionMap.png', dpi=200, logdir='./logs',
                        lonlim=[-132, -126.5], latlim=[50.5, 52],
                        topofile = '~cproof/Sites/gliderdata/smith_sandwell_topo_v8_2.nc',
                        start=np.datetime64('1970-01-01')):

    utcnow = datetime.datetime.utcnow()

    Calvert = get_line_info(linename=linename)

    # get positions, times, and ampH from logfiles:
    fns = glob.glob(f'{logdir}/*.log')
    fns.sort()
    print(fns)
    glider = parse_logfiles(fns)
    glider = glider.dropna(dim='surfacing')
    glider = glider.sel(surfacing=(glider.time>start))
    # get a distance along line.  Simple interp in lon which is prob OK here
    glider['Calvertdist'] = ('surfacing', np.interp(glider['lon'],
                            Calvert['wps'][::-1, -1], Calvert['dist'][::-1]))

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
    glider6h['distCalvert'] = ('time', np.interp(glider6h['lon'], Calvert['wps'][::-1, -1], Calvert['dist'][::-1]))
    glider6h['speedCalvert'] = ('timemid', glider6h.distCalvert.diff(dim='time').values / dt.astype(float) * 24 * 3600)

    inds = min(5, len(glider6h.speedCalvert))
    latestCalvertspeed = np.mean(glider6h.speedCalvert[-inds:].values)
    latestAmpHperday = np.mean(glider6h.ampHperday[-inds:].values)
    latestAmpHper100km = np.abs(latestAmpHperday / latestCalvertspeed * 100)

    # string with current time:
    timest = f'{glider6h.time.values[-1]}'[:16]
    timest = timest.replace('T', ' ')

    # time since last:
    print(glider6h.time.values[-1].astype('datetime64[s]').tolist())
    sincelast = utcnow - (glider6h.time.values[-1].astype('datetime64[s]')).tolist()
    sincelast = str(sincelast)[:-10]

    # figure out how far we have to go:
    total = (Calvert['dist'][0] - Calvert['dist'][-1]) * 2
    dist_to_Calvert = glider6h.distCalvert[-1] - Calvert['dist'][-1]
    print(total, dist_to_Calvert.values)
    if glider6h.speedCalvert[-6] > 0:
        # we are returning
        dist_to_go = Calvert['dist'][0] - glider6h.distCalvert[-1]
    else:
        dist_to_go = total / 2 + dist_to_Calvert

    # arrival time:
    try:
        predTime = glider6h.time[-1].values + np.timedelta64(int(dist_to_go / np.abs(latestCalvertspeed) * 24 * 3600), 's')
        predTime = f'{predTime}'[:10]
    except (OverflowError):
        predTime ='Bad data'


    # Plots:

    fig, axs = plt.subplots(4, 1, constrained_layout=True, figsize=(8, 10),
                            gridspec_kw={'height_ratios':[1, 1, 1, 2]})

    axs[0].plot(glider6h.timemid, glider6h.speed, label='Absolute Speed')
    axs[0].plot(glider6h.timemid, -glider6h.speed, color='C0')
    axs[0].plot(glider6h.timemid, glider6h.speedCalvert, label='Speed along Line')
    axs[0].plot(glider6h.timemid[-inds:], latestCalvertspeed * np.ones(inds), color='C1', alpha=0.5, lw=4)
    text = f'Latest speed (30 h): {latestCalvertspeed:1.1f} km/d\n'
    text += f'Predicted return time: {predTime}'
    axs[0].text(0.1, 0.5, text,
                transform=axs[0].transAxes)
    axs[0].grid('True')
    axs[0].legend()
    axs[0].set_ylabel('Speed [km/d]')
    axs[0].set_ylim([-50, 50])
    axs[0].set_title(f'{timest}; Processed: {str(utcnow)[:16]}\n{sincelast} since last update')

    axs[1].plot(glider6h.timemid, glider6h.ampHperday)
    axs[1].plot(glider6h.timemid[-inds:], latestAmpHperday * np.ones(inds), color='C0', alpha=0.5, lw=4)
    axs[1].grid('True')
    axs[1].set_ylabel('AmpH/d')
    axs[1].set_ylim(bottom=0, top=15)
    text =  f'Latest [30 h]:  {latestAmpHperday:4.1f} AmpH/d\n'
    text += f'Used [230 Amph] {glider6h.ampH.values[-1]:4.1f} Amph\n'
    text += f'Percent used:   {glider6h.ampH.values[-1] / 2.3:4.1f}%'
    axs[1].text(0.05, 0.1, text, family='monospace', fontsize='medium',
                        transform=axs[1].transAxes)

    axs[2].plot(glider6h.timemid, 100 * glider6h.ampHperday/glider6h.speed)
    axs[2].plot(glider6h.timemid, np.abs(100 * glider6h.ampHperday/glider6h.speedCalvert))
    axs[2].plot(glider6h.timemid[-inds:], latestAmpHper100km * np.ones(inds), color='C1', alpha=0.5, lw=4)
    axs[2].set_ylim([0, 100])
    axs[2].grid('True')
    axs[2].set_ylabel('AmpH / 100 km')

    predAH = glider6h.ampH[-1].values + latestAmpHper100km*dist_to_go.values / 100
    print(dist_to_go / np.abs(latestCalvertspeed))
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

    td = list(range(1, len(Calvert['wps']), 1))
    td = [0] + td + [-1]
    for wpind in td:
        ax.plot(Calvert['wps'][wpind, 1], Calvert['wps'][wpind, 0], 'd', color='C1')

        ax.text(Calvert['wps'][wpind, 1],
                Calvert['wps'][wpind, 0],
                f" {Calvert['nm'][wpind]} {Calvert['dist'][wpind]:1.0f} km",
                color='C1', rotation=60, fontsize='smaller')
    ax.plot(glider.lon, glider.lat, 'm.', markersize=1)
    ax.plot(glider.lon[-1], glider.lat[-1], 'go', markersize=6)

    text = f'Distance travelled: {(total - dist_to_go.values):1.0f} km\n'
    text += f'Distance to go: {dist_to_go.values:1.0f} km'
    ax.text(0.05, 0.1, text, family='monospace', fontsize='medium',
                            transform=ax.transAxes, bbox={'facecolor':[1, 1, 1]})

    fig.savefig(f'{figdir}/{outname}', dpi=dpi)
