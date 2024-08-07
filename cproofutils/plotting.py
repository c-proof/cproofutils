import xarray as xr
import logging
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import yaml
import matplotlib.dates as mdates
import matplotlib.units as munits
import gsw
from datetime import datetime
import glob

import cmocean

_log = logging.getLogger(__name__)

try:
    converter = mdates.ConciseDateConverter()
    munits.registry[np.datetime64] = converter
except:
    # older matplotlib...
    pass

plt.rcParams['figure.constrained_layout.use']=True

def _autoclim(vals):
    vals = vals.values.flatten()
    vals = vals[~np.isnan(vals)]
    min = np.min(vals)
    max = np.max(vals)

    m15 = np.quantile(vals, 0.02 / 2)
    m85 = np.quantile(vals, 1 - 0.02 / 2)
    d = m85 - m15
    m15 = m15 - d / 10
    m85 = m85 + d / 10
    _log.info(f'm15, m85 {m15}, {m85}, {np.max((min, m15))}')
    return np.max((min, m15)), np.min((max, m85))

def timeseries_plots(fname, plottingyaml):

    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)
        if 'starttime' in config.keys():
            starttime = np.datetime64(config['starttime'], 'ns')
        else:
            starttime = None
        _log.info(f'starttime: {starttime} {starttime.dtype}')

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    with xr.open_dataset(fname, decode_times=True) as ds0:
        ds = ds0.where(ds0.time>starttime, drop=True).dropna(dim='time')

        # map!
        fig1, axs1 = plt.subplots(figsize=(7.5, 7.5))
        ax = axs1
        good = (ds.longitude < -120)
        ax.plot(ds.longitude[good], ds.latitude[good], '.',markersize=1)

        #ax.plot(ds.longitude, ds.latitude, '.') #jp commented
       # ax.set_aspect(1 / np.cos(np.deg2rad(ds.latitude.mean())))
        ax.set_ylabel('Lat [degrees north]')
        ax.set_xlabel('Lon [degrees west]')
        ax.grid()
        fig1.savefig(config['figdir'] + '/jp1_map%s.png'%ds.attrs['deployment_name'], dpi=200)

        fig, axs = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 1]})
        ax = axs[0]
        ax.plot(ds.time[good], ds.longitude[good], '.',markersize=1)
        ax.set_ylabel('Lon [degrees west]')
        ax.grid()

        ax = axs[1]
        ax.plot(ds.time[good], ds.latitude[good], '.',markersize=1)
        ax.set_ylabel('Lat [degrees north]')
        ax.grid()
        fig.savefig(config['figdir'] + '/jp2_map%s.png'%ds.attrs['deployment_name'], dpi=200)

        # timeseries of things....
        _log.info('Plotting timeseries data')
        keys = config['timeseries'].keys()
        N = len(keys)
        if 1:
            fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                    sharex=True, sharey=False)
            axs = axs.flat
            for n, k in enumerate(keys):
                a = [0,1,2,3]
                b = [4,5]
                _log.debug(f'key {k}')
                if (config['timeseries'][k] == 'True') and (k in ds):
                    ax = axs[n]
                    #convert 9999 to nan
                    ds[k] = ds[k].where(ds[k] != 9999, np.nan)
                    good = np.where(~np.isnan(ds[k]))[0]
                    if good.size != 0:
                        pc = ax.plot(ds.time[good], ds[k][good], '.',markersize=1)
                        min, max = _autoclim(ds[k][good])
                        ax.set_ylim(min, max)
                        ax.grid()
                        ax.set_ylabel(ds[k].attrs['long_name'] + ' [' +
                            ds[k].attrs['units'] + ']')
                    if n in a:
                        ax.xaxis.set_tick_params(labelbottom=True)
                    ax.set_title(ds[k].attrs['long_name'] + ' over time', loc='left', fontsize=9)
            fig.savefig(config['figdir'] + '/ts_%s.png'%ds.attrs['deployment_name'], dpi=200)

        _log.info('Plotting timeseries data')
        keys = config['timeseries'].keys()
        N = len(keys)
        if 1:
            fig, axs = plt.subplots(int(N / 2), 2, figsize=(7, 12),
                                    sharex=False, sharey=False)
            axs = axs.flat
            for n, k in enumerate(keys):
                _log.debug(f'key {k}')
                print(f'key {k}')
                if (config['timeseries'][k] == 'True') and (k in ds):
                    ax = axs[n]
                    if ax == axs[1]:
                        ax.set_xlim(30,35) #jp hardcode?
                       # ds[k] = ds[k].where(ds[k] != 0, np.nan)
                       # good = np.where(~np.isnan(ds[k]))[0]
                    if ax == axs[3]:
                        ax.set_xlim(0,0.004) #jp hardcode
                    if ax == axs[4]:
                        ax.set_xlim(0,7.5)
                    if ax == axs[5]:
                        ax.set_xlim(0,4)

                    good = np.where(~np.isnan(ds[k]))[0]
                    pc = ax.plot(ds[k][good],ds.depth[good], '.',markersize=1)
                    ax.grid()
                    ax.set_ylabel('DEPTH [m]')
                    ax.xaxis.set_tick_params(labelbottom=True)
                    ax.invert_yaxis()

                    ax.set_xlabel(ds[k].attrs['long_name'] + ' [' +
                            ds[k].attrs['units'] + ']')
                    ax.set_title(ds[k].attrs['long_name'] + ' over depth', loc='left', fontsize=9)

            fig.savefig(config['figdir'] + '/vv_%s.png'%ds.attrs['deployment_name'], dpi=200)


       # _log.info('Plotting timeseries data')
        keys = config['timeseries'].keys()
       # N = len(keys)
       ##water temp v depth plot
        if 1:
            #fig, axs = plt.subplots(int(N / 2), 2,figsize=(7.5, 7),
                             #       sharex=True, sharey=False)
            fig, axs = plt.subplots(1, 1,figsize=(7.5, 7),
                                    sharex=True, sharey=False)
          #  axs = axs.flat
            for n, k in enumerate(keys):
                if n == 0:
                    a = [0,1,2,3]
                    b = [4,5]
                    #_log.debug(f'key {k}')
                    #print(f'key {k}')
                    if config['timeseries'][k] == 'True':
                        ax = axs

                        good = np.where(~np.isnan(ds[k]))[0]
                        pc = ax.scatter(ds.temperature[good], ds.depth[good], c=ds.longitude[good], cmap='gray', s=1)

                        ax.grid()
                        ax.set_ylabel('DEPTH [m]')
                        x = np.random.randint(low=0, high=10, size=13)
                        plt.xticks(np.arange(4, len(x)+1, 1))

                        ax.set_title('water temp over depth')

                        ax.set_xlabel(ds[k].attrs['long_name'] + ' [' +
                            ds[k].attrs['units'] + ']')
                        plt.gca().invert_yaxis()

            cbar = fig.colorbar(pc, ax=axs)
            cbar.ax.set_yticklabels(["{:.4}".format(i) for i in cbar.get_ticks()])
            cbar.ax.set_ylabel('Longitude [degrees west]')
            fig.savefig(config['figdir'] + '/jp_ts_%s.png'%ds.attrs['deployment_name'], dpi=200)

        # colorline?
        if 0:
            fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                    sharex=True, sharey=True)
            axs = axs.flat
            _log.info('Plotting colorline data')

            for n, k in enumerate(keys):
                _log.debug(f'key, {k}')
                ax = axs[n]
                locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
                formatter = mdates.ConciseDateFormatter(locator)
                ax.xaxis.set_major_locator(locator)
                ax.xaxis.set_major_formatter(formatter)
                good = np.where(~np.isnan(ds[k]))[0]
                min, max = _autoclim(ds[k][good])
                pc = ax.scatter(ds.time[good], ds.depth[good], s=3, c=ds[k][good],
                    rasterized=True, vmin=min, vmax=max)

                fig.colorbar(pc, ax=ax, extend='both', shrink=0.6)
                ax.set_title(ds[k].attrs['long_name'] + ' [' +
                             ds[k].attrs['units'] + ']', loc='left', fontsize=9)
                t0 = ds.time[good][0]
                t1 = ds.time[good][-1]
                ax.set_xlim([t0, t1])
                ax.set_ylim([ds['depth'].max(), ds['depth'].min()])
                if n == 0:
                    ax.set_ylabel('DEPTH [m]')
            fig.savefig(config['figdir'] + '/cl_%s.png'%ds.attrs['deployment_name'], dpi=200)

            # depth profiles...
            fig, axs = plt.subplots(int(N / 2), 2, figsize=(7.5, 7),
                                    sharex=True, sharey=True)
            axs = axs.flat
            _log.info('Plotting colorline data')

            for n, k in enumerate(keys):
                _log.info('cl %s', k)
                if config['timeseries'][k] == 'True':
                    ax = axs[n]
                    locator = mdates.AutoDateLocator(minticks=3, maxticks=7)
                    formatter = mdates.ConciseDateFormatter(locator)
                    ax.xaxis.set_major_locator(locator)
                    ax.xaxis.set_major_formatter(formatter)
                    good = np.where(~np.isnan(ds[k] + ds.depth))[0]
                    min, max = _autoclim(ds[k][good])
                    pc = ax.scatter(ds.time[good], ds.depth[good], s=3, c=ds[k][good],
                        rasterized=True)
                    fig.colorbar(pc, ax=ax, extend='both', shrink=0.6)
                    ax.set_title(ds[k].attrs['long_name'] + ' [' +
                                 ds[k].attrs['units'] + ']', loc='left', fontsize=9)
                    t0 = ds.time[good][0]
                    t1 = ds.time[good][-1]
                    ax.set_xlim([t0, t1])
                    ax.set_ylim([ds['depth'].max(), ds['depth'].min()])
                    if n == 0:
                        ax.set_ylabel('DEPTH [m]')
            fig.savefig(config['figdir'] +
                        '/cl_%s.png'%ds.attrs['deployment_name'], dpi=200)

        # prop_v_prop:
        _log.info('property vs property')

        fig, axs = plt.subplots(1, 2, constrained_layout=True, sharey=True)
        ax = axs[0]

        good = ~np.isnan(ds.salinity)
        if len(good) < 2:
            return
        smin, smax = _autoclim(ds.salinity[good])
        good = ~np.isnan(ds.temperature)
        tmin, tmax = _autoclim(ds.temperature[good])

        s = np.linspace(ds.salinity.min(), ds.salinity.max(), 100)
        t = np.linspace(ds.temperature.min(), ds.temperature.max(), 100)
        S, T = np.meshgrid(s, t)
        sa = gsw.SA_from_SP(S, T*0, ds['longitude'].mean().values, ds['latitude'].mean().values)
        ct = gsw.CT_from_t(sa, T, T*0)
        pd = gsw.density.sigma0(sa, ct)
        levels = np.arange(2, 30, 0.5)
        c = ax.contour(s, t, pd, colors='0.5', levels=levels)
        ax.clabel(c, levels[::2], fontsize=8, fmt='%1.1f')
        ax.plot(ds['salinity'], ds['temperature'], '.', markersize=1)
        ax.set_xlabel('$S\ [psu]$')
        ax.set_ylabel('$T\ [^oC]$')
        ax.set_xlim(smin, smax)
        ax.set_ylim(tmin, tmax)
        ax.grid()
        try:
            ax = axs[1]
            ds['oxygen_concentration'][ ds['oxygen_concentration']==9999]=np.nan

            ax.plot(ds['oxygen_concentration'], ds['temperature'],
                    '.', markersize=1)
            ax.set_xlabel('$O^2\ [mmol/L]$')
            ax.grid()
        except:
            pass
        add_suptitle(fig, ds)
        fig.savefig(config['figdir'] +
            '/pvp_%s.png'%ds.attrs['deployment_name'], dpi=200)

def add_suptitle(fig, ds):
    sst = 'Glider/deployment: %s; ' % ds.attrs['deployment_name']
    sst += 'lat = %1.3f to %1.3f N; ' % (ds.attrs['geospatial_lat_min'],
                                  ds.attrs['geospatial_lat_max'])
    sst += 'lon = %1.3f to %1.3f E ' % (ds.attrs['geospatial_lon_min'],
                                  ds.attrs['geospatial_lon_max'])
    sst += '\n CPROOF: http://cproof.uvic.ca/; Preliminary data'

    fig.suptitle(sst, fontsize=8, fontfamily='courier')


def grid_plots(fname, plottingyaml):
    _log.info('Grid plots!')
    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)
        if 'starttime' in config.keys():
            starttime = np.datetime64(config['starttime'], 'ns')
        else:
            starttime = None

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    with xr.open_dataset(fname, decode_times=True) as ds:
        # ds = ds0.sel(time=slice(starttime, None))
        ds = ds.where(ds.time>starttime, drop=True)
        _log.debug(str(ds))

        keys = config['pcolor']['vars'].keys()
        N = len(keys)
        # get the max depth that data is at:
        tmean = ds.temperature.mean(axis=1)
        indmax = np.where(~np.isnan(tmean))[0][-1]
        depmax = ds.depth[indmax]
        fig, axs = plt.subplots(int(math.ceil(N / 2)), 2, figsize=(7.5, 7),
                                sharex=True, sharey=True, constrained_layout=True)
        axs = axs.flat
        for n, k in enumerate(keys):
            if not k in ds.keys():
                _log.warning(f'{k} not in dataset')
                continue
            _log.debug(f'key {k}')
            pconf = config['pcolor']['vars'][k]
            _log.debug(pconf)

            cmap = pconf.get('cmap','viridis')
            vmin = pconf.get('vmin', None)
            vmax = pconf.get('vmax', None)

            ax = axs[n]

            min, max = _autoclim(ds[k])
            if vmin is not None:
                min = vmin
            if vmax is not None:
                max = vmax
            # get good profiles.  i.e. those w data
            ind = np.where(np.sum(np.isfinite(ds[k].values), axis=0)>2)[0]
            _log.debug(ind)
            _log.debug(len(ds.time))
            if len(ind) > 2:
                time = ds.time[ind[1:]] + np.diff(ds.time[ind]) / 2
                print('TIME', time)
                time = np.hstack((time[0] - (time[1]-time[0]) / 2, time))
                depth = ds.depth[:-1] - np.diff(ds.depth)
                depth = np.hstack((depth, depth[-1] + np.diff(ds.depth)[-1]))
                pc = ax.pcolormesh(time, depth, ds[k][:, ind],
                    rasterized=True, vmin=min, vmax=max, cmap=cmap, shading='auto')
                ax.contour(ds.time[ind], ds.depth, ds.potential_density[:, ind], colors='0.5',
                       levels=np.arange(22, 28, 0.5)+1000, linewidths=0.5, alpha=0.7)
                fig.colorbar(pc, ax=ax, extend='both', shrink=0.6, pad=0.01)
            ax.set_title(ds[k].attrs['long_name'] + ' [' +
                         ds[k].attrs['units'] + ']', loc='left', fontsize=9)
            t0 = ds.time[0]
            t1 = ds.time[-1]
            ax.set_xlim([t0, t1])
            ax.set_ylim([depmax, 0])
            if n == 0:
                ax.set_ylabel('DEPTH [m]')
            ax.set_facecolor('0.8')
        now = str(datetime.utcnow())[:-10]
        lastdata = str(ds.time[-1].values)[:-13]
        fig.suptitle(f'Processed: {now}, Lastdata: {lastdata} ')
        fig.savefig((config['figdir'] + '/' +
                     'pcolor_%s.png'%ds.attrs['deployment_name']), dpi=200)
        ax.set_ylim([400, 0])
        fig.savefig((config['figdir'] + '/' +
                     'pcolor400_%s.png'%ds.attrs['deployment_name']), dpi=200)
        ax.set_ylim([depmax, 0])

        t1 = ds.time[-1]
        t0 = t1 - np.timedelta64(10, 'D')
        for ax in axs:
            ax.set_xlim([t0, t1])
        fig.savefig(config['figdir'] + '/pcolor_%s_last10d.png'%ds.attrs['deployment_name'], dpi=200)


def anomaly_plots(fname, plottingyaml):
    """
    Temperature anomalies from LineP means and standard deviations.
    """
    _log.info('Anomaly plots!')
    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)
        if 'starttime' in config.keys():
            starttime = np.datetime64(config['starttime'], 'ns')
        else:
            starttime = None

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    fig, ax = plt.subplots(layout='constrained', figsize=(4.5, 2.7))
    cmap = plt.cm.RdBu_r

    with (xr.open_dataset(fname, decode_times=True) as ds,
          xr.open_dataset('../../../meanProfilesLineP.nc', decode_times=True) as dsmean):
        ds = ds.where(ds.time>starttime, drop=True)

        ind = np.where(np.sum(np.isfinite(ds['temperature'].values), axis=0)>3)[0]
        ds = ds.isel(time=ind)
        tmean = ds.temperature.mean(dim='time', skipna=True)
        indmax = np.where(~np.isnan(tmean))[0][-1]
        depmax = ds.depth[indmax]

        ds['Tmean'] = ds.temperature * 0
        ds['Tstd'] = ds.temperature * 0
        for j in range(len(ds.time)):
            ds['Tmean'][:, j] = np.interp(ds.potential_density[:, j],
                                          dsmean.potential_density, dsmean.temperature)
            ds['Tstd'][:, j] = np.interp(ds.potential_density[:, j],
                                         dsmean.potential_density, dsmean.temperature_std)

        pc = ax.pcolormesh(ds.time, ds.depth, (ds.temperature - ds.Tmean)/ds.Tstd,
                           rasterized=True, vmin=-3, vmax=3, cmap=cmap, shading='nearest')
        ax.contour(ds.time, ds.depth, ds.potential_density, colors='0.5',
                   levels=np.arange(22, 28, 0.5)+1000, linewidths=0.5, alpha=0.7)

        fig.colorbar(pc, ax=ax, extend='both', shrink=0.6, pad=0.01, label='normalized temperature anom.')
        ax.set_title(ds.attrs['deployment_name'], loc='left', fontsize=9, fontstyle='italic')
        t0 = ds.time[0]
        t1 = ds.time[-1]
        ax.set_xlim([t0, t1])
        ax.set_ylim([depmax, 0])
        ax.set_facecolor('0.6')
        ax.set_ylabel('DEPTH [m]')
        fig.savefig(config['figdir'] + '/pcolor_%s_tempanom.png'%ds.attrs['deployment_name'], dpi=200)


def overview_plot(fname, plottingyaml, location='LineP',
                  to_try=None,
                  topofile='~cproof/Sites/gliderdata/smith_sandwell_topo_v8_2.nc'):
    """
    Make a plot w/ Temp, O2, and Chl, and a map in fourth place.

    Map location is set by "location" variable.
    """
    lonlim=[-145.2, -124.4]
    latlim=[46.1, 52.2]
    _log.info('Overview plots!')
    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)
        if 'starttime' in config.keys():
            starttime = np.datetime64(config['starttime'], 'ns')
        else:
            starttime = None

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    with xr.open_dataset(fname, decode_times=True) as ds:
        # ds = ds0.sel(time=slice(starttime, None))
        ds = ds.where(ds.time>starttime, drop=True)
        _log.debug(str(ds))

        keys = config['pcolor']['vars'].keys()
        # choose the first three of these
        if to_try is None:
            to_try = ['temperature','oxygen_concentration', 'chlorophyll', 'backscatter_700', 'cdom', 'salinity']
        new_keys = []
        for key in to_try:
            if (key in keys) and (key in ds.keys()):
                new_keys += [key]
        keys = new_keys[:3]
        # get the max depth that data is at:
        tmean = ds.temperature.mean(axis=1)
        indmax = np.where(~np.isnan(tmean))[0][-1]
        depmax = ds.depth[indmax]
        fig, axs = plt.subplots(2, 2, figsize=(5*8/4, 5),
                                constrained_layout=True)
        axs = axs.flat
        for n, k in enumerate(keys):

            _log.debug(f'key {k}')
            pconf = config['pcolor']['vars'][k]
            _log.debug(pconf)

            cmap = pconf.get('cmap','viridis')
            vmin = pconf.get('vmin', None)
            vmax = pconf.get('vmax', None)

            ax = axs[n]

            min, max = _autoclim(ds[k])
            if vmin is not None:
                min = vmin
            if vmax is not None:
                max = vmax
            # get good profiles.  i.e. those w data
            ind = np.where(np.sum(np.isfinite(ds[k].values), axis=0)>2)[0]
            _log.debug(ind)
            _log.debug(len(ds.time))
            if len(ind) > 2:
                time = ds.time[ind[1:]] + np.diff(ds.time[ind]) / 2
                print('TIME', time)
                time = np.hstack((time[0] - (time[1]-time[0]) / 2, time))
                depth = ds.depth[:-1] - np.diff(ds.depth)
                depth = np.hstack((depth, depth[-1] + np.diff(ds.depth)[-1]))
                pc = ax.pcolormesh(time, depth, ds[k][:, ind],
                    rasterized=True, vmin=min, vmax=max, cmap=cmap, shading='auto')
                ax.contour(ds.time[ind], ds.depth, ds.potential_density[:, ind], colors='0.5',
                       levels=np.arange(22, 28, 0.5)+1000, linewidths=0.5, alpha=0.7)
                fig.colorbar(pc, ax=ax, extend='both', shrink=0.6)
            ax.set_title(ds[k].attrs['long_name'] + ' [' +
                         ds[k].attrs['units'] + ']', loc='left', fontsize=9)
            t0 = ds.time[0]
            t1 = ds.time[-1]
            ax.set_xlim([t0, t1])
            ax.set_ylim([depmax, 0])
            if n == 0:
                ax.set_ylabel('DEPTH [m]')
            ax.label_outer()
            ax.set_facecolor('0.8')


        # map:
        figmap = fig.add_subfigure(axs[-1].get_subplotspec())
        axs[-1].set_visible(False)
        latlims = [46.25, 53.6]
        lonlims = [145.1, 123]
        import matplotlib.colors as mcolors
        bounds = [-6500, -6000, -5000, -4500,-4000, -3500, -3000, -2500, -2000, -1500, -1000, -750, -500, -350, -200, -150,  -100, -75, -50, -25, 0]
        ax = figmap.add_subplot([0, 0, 1, 1])
        norm = mcolors.BoundaryNorm(bounds, 256)
        with xr.open_dataset(topofile) as topo:
            # longitude, latitiude, ROSE
            topo['longitude'] = -(topo.longitude - 360)
            topo = topo.sel(latitude=slice(latlims[0], latlims[1]), longitude=slice(*lonlims))

            ax.contourf(topo.longitude, topo.latitude, topo.ROSE, linewidths=1, cmap='Greens_r',
                    levels=np.arange(0, 20_000, 500), vmax=2000, vmin=-200, alpha=0.7, zorder=-10)

            ax.pcolormesh(topo.longitude, topo.latitude, topo.ROSE, linewidths=0, cmap='Blues_r', zorder=-100,
                        rasterized=True, norm=norm, alpha=0.7)

        ax.set_xlim(lonlims)
        ax.set_ylim(latlims)

        ax.plot(-ds.longitude, ds.latitude, '.', markersize=1, color='C1', zorder=100)
        ax.plot(-ds.longitude[-1], ds.latitude[-1], 'go', zorder=100)
        wp = {}
        wp['P26'] = [145, 50]
        wp['P16'] = [134+4/6, 49+17/60]
        wp['P4'] =  [126+40/60, 48+39/60]
        for w in wp:
            ax.plot(wp[w][0], wp[w][1], 'md', zorder=2, markersize=4)
            ax.text(wp[w][0], wp[w][1], '  ' + w, color='m', zorder=1, fontsize='xx-small', fontstyle='italic', rotation=45)
        #cb.ax.tick_params(labelsize='small')

        ax.set_xlabel('Longitude $\mathrm{[^o W]}$')
        ax.set_ylabel('Latitude $\mathrm{[^o N]}$')
        ax.text(140, 49.14, 'L   i   n   e       P', rotation=-7.0, va='bottom', rotation_mode='anchor', fontstyle='italic',  fontsize='x-small')
        ax.text(129.45, 46.55, 'Southern Line', rotation=33, va='bottom', rotation_mode='anchor', fontstyle='italic',  fontsize='x-small')
        ax.text(133, 49.45, 'Calvert Line', rotation=30, va='bottom', rotation_mode='anchor', fontstyle='italic', fontsize='x-small')
        ax.set_aspect(1/np.cos(49 * np.pi / 180))

        with xr.open_dataset('../../CProofTracks.nc') as dst:
            dst = dst.where((dst.longitude>-150) & (dst.longitude< -120)
                            & (dst.latitude>35) & (dst.latitude<60),
                            drop=True)
            print(dst)
            ax.plot(-dst.longitude[::10], dst.latitude[::10], '.', color='0.5', markersize=1, zorder=-90)
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position("right")
        # info at top

        t0 = str(ds.time[0].values)[:16].replace('T', ' ')
        t1 = str(ds.time[-1].values)[:16].replace('T', ' ')
        # active?
        now = np.datetime64(datetime.utcnow())
        if now - ds.time[-1] < np.timedelta64(3, 'D'):
            active = "(ongoing)"
        else:
            active = ""

        print(ds.attrs)
        fig.suptitle(f'C-PROOF {ds.attrs["deployment_name"]}   {ds.attrs["project"]}   {t0} to {t1} {active}')

        fig.savefig((config['figdir'] + '/' +
                     'overview_pcolor_%s.png'%ds.attrs['deployment_name']), dpi=200)



def ts_dens_plots(fname, plottingyaml):
    """
    Temperature anomalies from LineP means and standard deviations.
    """
    _log.info('Anomaly plots!')
    with open(plottingyaml) as fin:
        config = yaml.safe_load(fin)
        if 'starttime' in config.keys():
            starttime = np.datetime64(config['starttime'])
        else:
            starttime = None

    try:
        os.mkdir(config['figdir'])
    except:
        pass

    fig, ax = plt.subplots(layout='constrained', figsize=(4.0, 4.0))
    cmap = plt.cm.RdBu_r

    with (xr.open_dataset(fname, decode_times=True) as ds0,
          xr.open_dataset('../../../meanProfilesLineP.nc', decode_times=True) as dsmean):
        pass

