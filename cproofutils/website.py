import yaml
import glob
import os
import xarray as xr
import numpy as np
from jinja2 import Environment, FileSystemLoader
import geojson
import json
import datetime
import simplekml
import pyglider.utils as pygu
from pathlib import Path

import logging

import argopandas as argo
import pandas as pd

"""
Utilities to make smart directories for the websites.

A driving script may look like

```
import pyglider.website as pyweb
import logging

logging.basicConfig(level=logging.DEBUG)
# uses the templat in `.templates/deploymentsIndex.html` to render
# index.html in each of the subdirectories:

pyweb.geojson_deployments('./')

if 1:
    pyweb.index_deployments('./dfo-walle652/')
    pyweb.index_deployments('./dfo-bb046/')
```

where `dfo-walle652` will contain subdirectories, each with a deployment:

```
./dfo-walle652/dfo-walle652-20210903
./dfo-walle652/dfo-walle652-20220101
...
```

Output are index files in each subdirectory and a geojson of all
the deployments (if that is desired)

```
./deployments.geojson   # from the first script....
./dfo-walle652/index.html
./dfo-walle652/dfo-walle652-20210903/index.html
...

Requires an html template.  The ones used for cproof are in
../example_html_templates of this repo and you need
to supply

`deploymentsIndex.html`
`deploymentsInfoNew.html`
and maybe `deploymentsInfo.html`

These are example templates, and written in jinga templating language.  In
particular, see the strings near the bottom with double braces: eg:
"{{ att['deployment_name'] }}"  These are what get filled by calling
`pyweb.index_deployments`.

Note there are some hardcoded things in this tool, so please read the source
code!

"""

_log = logging.getLogger(__name__)

def ticker(dir):

#jpnote:
# Function to calculate distance traveled by gliders
# For every mission directory where L0-timeseries exists,
# extract distance_over_ground and sum
# Written into html file, and put directly into cproofwesbite/_includes/
# directory, and then added using jekyll {% include %}

    sump = 0                   # set variable to 0
    nprofiles = 0
    totaltime = 0
    t0 = np.datetime64('2200-01-01')
    tf = np.datetime64('1900-01-01')
    totaltime = 0
    t0 = np.datetime64('2200-01-01')
    tf = np.datetime64('1900-01-01')
    gliders = [d for d in glob.glob(dir + '/dfo-*') if os.path.isdir(d)]
    for glider in gliders:
        dirs = [d for d in glob.glob(glider + '/dfo-*') if os.path.isdir(d)]
        for d in dirs:
            _log.info(d)
            nc = sorted(glob.glob(d+'/L0-timeseries/*.nc'), key=os.path.getmtime)
            if len(nc) < 1:
                nc = sorted(glob.glob(d+'/L1-timeseries/*.nc'), key=os.path.getmtime)
                if len(nc) < 1:
                    _log.info('No data!')
                    continue
            with xr.open_dataset(nc[-1]) as ds:
                if ds.time[0].values < t0:
                    t0 = ds.time[0].values
                if ds.time[-1].values > tf:
                    tf = ds.time[-1].values

                dttime = (ds.time[-1] - ds.time[0]) / 1e9 / 24 / 3600
                _log.info('Time %f', dttime)
                totaltime += dttime
                totalkm = 0
                #open most recent of sorted netcdf files (last, [-1])
                # Add last element of dataset to final sum
                totalkm = ds.variables['distance_over_ground'].data[-1]
                # sometimes the distance over ground resets, so we need to find all
                # the indices where it goes back to zero, and add the DOG before:
                list_ = np.where(ds.variables['distance_over_ground'].data==0)[0]

                for i in list_:
                    if i>0:
                        totalkm += ds.variables['distance_over_ground'].data[i-1]
                _log.info(f'Total km: {totalkm}')
                sump += totalkm
                try:
                    newprofiles = len(np.unique(ds.profile_index))
                except:
                    pass
                _log.info(f'Total profiles: {newprofiles}')
                nprofiles += newprofiles

    with open('/Users/cproof/processing/deployments/Totals.html', 'w') as output_file:
        outstr = f"""
Between {str(t0)[:10]} and {str(tf)[:10]}, our gliders have traveled {int(sump):,} km and
made {int(nprofiles):,} CTD, O2, and optics casts, over {int(totaltime):,} days at sea.
"""
        _log.info(outstr)
        output_file.write(outstr)

#######

def update_active_json(glider, mission):
    file_path = './active_deployments_new.json'
    new_data = ( glider + '/'
                + mission + '/figs/overview_pcolor_' + mission +'.png')
    if not os.path.exists(new_data):
        _log.warning('Image %s not found', new_data)
        return
    new_data = 'https://cproof.uvic.ca/gliderdata/deployments/' + new_data
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        # Read existing data
        with open(file_path, 'r') as file:
            try:
                data = json.load(file)
            except json.JSONDecodeError:
                data = []
    else:
        # File does not exist or is empty
        data = []

    url = ('https://cproof.uvic.ca/gliderdata/deployments/' + glider + '/'
                + mission + '/figs/overview_pcolor_' + mission +'.png')
    new_data = {'url': url, 'duration': 600}
    # Append new data
    if isinstance(data, list):
        data.append(new_data)
    elif isinstance(data, dict):
        data.update(new_data)

    # Write updated data back to the file
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)
        _log.info('Wrote to %s', file_path)


def index_deployments(dir, templatedir='./.templates/',
                      listdict={'mission_all.txt':None}):
    """
    Get useful info from deployments under "dir" and add to an
    index.html in "dir"

    The structure is meant to be simple:

    dir/deployment1/deployment.yml
    dir/deployment2/deployment.yml
    dir/deployment3/deployment.yml

    Makes a file `dir/index.html` that has table of the deployments.
    """

    file_loader = FileSystemLoader(templatedir)
    env = Environment(loader=file_loader)
    template = env.get_template('deploymentsIndex.html')

    activef = open('./active_deployments_new.txt', 'a+')

    subdirs = glob.glob(dir + '/*')
    atts = []
    for d in subdirs:
        if os.path.isdir(d):
            if 1:
                _log.info(d)
                nc = glob.glob(d+'/L0-timeseries/*.nc')
                if len(nc) < 1:
                    nc = glob.glob(d+'/L1-timeseries/*.nc')
                    if len(nc) < 1:
                        continue
                with xr.open_dataset(nc[0]) as ds:
                    att = ds.attrs
                    atts.append(att)
            else:
                pass
    if len(atts) > 0:
        output = template.render(atts=atts,
            title=atts[-1]['glider_name'] + atts[-1]['glider_serial'])
        with open(dir + '/index.html', 'w') as fout:
            fout.write(output)

    # now make the individual deployment index pages...
    templateOld = env.get_template('deploymentsInfo.html')
    templateNew = env.get_template('deploymentsInfoNew.html')
    subdirs = glob.glob(dir + '/*')
    atts = []
    for d in subdirs:
        if os.path.isdir(d):
            _log.info(d)
            if 1:
                figs = glob.glob(d + '/figs/*.png')
                for n, fig in enumerate(figs):
                    figs[n] = './figs/' + os.path.split(fig)[1]
                    _log.info(figs[n])
                nc = glob.glob(d+'/L0-timeseries/dfo-*.nc')
                ncdelayed = glob.glob(d + '/L0-timeseries/dfo-*_delayed.nc')
                hasdelayed = len(ncdelayed) > 0
                nccorrected = glob.glob(d + '/L0-timeseries/dfo-*_corrected.nc')
                hascorrected = len(nccorrected) > 0
                template = templateNew

                if len(nc) < 1:
                    # try old style
                    nc = glob.glob(d+'/L1-timeseries/dfo-*.nc')
                    template = templateOld
                    if len(nc) < 1:
                        continue
                with xr.open_dataset(nc[0]) as ds:
                    _log.debug('Opening %s', nc[0])
                    att = ds.attrs
                    #_log.debug(type(ds.keys()))
                    keys = []
                    units = []
                    for key in ds.keys():
                        #_log.debug(ds[key].attrs)
                        try:
                            unit = ds[key].attrs['units']
                        except KeyError:
                            unit = 'no units'
                        keys.append(key + ' [' + unit +']')
                depname = att['deployment_name']
                glider = depname[:-9]
                ovp = glob.glob(d+f'/figs/overview_pcolor_{depname}.png')
                hasoverviewplot = len(ovp) > 0
                if hasoverviewplot:
                    hasoverviewplot = f'./figs/overview_pcolor_{depname}.png'
                else:
                    hasoverviewplot = ''
                output = template.render(deploy_name=depname,
                    title=att['glider_name'] + att['glider_serial'],
                    figs=figs, att=att, keys=keys,
                    hasdelayed=hasdelayed, hascorrected=hascorrected,
                    hasoverviewplot=hasoverviewplot)
                with open(d + '/index.html', 'w') as fout:
                    fout.write(output)
                # output to active_deployments.txt if active...
                if ds['time'].max() > np.datetime64(datetime.datetime.now()) - np.timedelta64(2, 'D'):
                    _log.debug('ACTIVE: %s', depname)
                    activef.write(d[2:]+'\n')
                    update_active_json(glider, depname)
                else:
                    _log.debug('INACTIVE: %s %s', depname, ds['time'][-1].values)

            else:
                pass
    activef.close()

def geojson_deployments(dir, outfile='cproof-deployments.geojson'):
    props = ['deployment_start', 'deployment_end', 'platform_type',
             'glider_model', 'glider_name', 'glider_serial',
             'deployment_name', 'project', 'institution', 'comment']
    subdirs = glob.glob(dir + '/*')
    features = []

    #kml = simplekml.Kml()

    np.random.seed(20190101)
    _log.debug(f'subdirs, {subdirs}')
    colornum = 0;
    filen = 0
    for d in subdirs:
        _log.info(d)
        if os.path.isdir(d):
            subdirs2 = glob.glob(d + '/*')
            for d2 in subdirs2:
                _log.info(d2)
                kml = simplekml.Kml()  #note: put kml in for loop
                if os.path.isdir(d2):
                    try:
                        nc = glob.glob(d2+'/L0-gridfiles/*grid_delayed.nc')
                        if len(nc) <1:
                            nc = glob.glob(d2+'/L0-gridfiles/*grid.nc')
                        if len(nc) < 1:
                            # old style
                            nc = glob.glob(d2+'/L2-gridfiles/*grid.nc')
                        if len(nc) < 1:
                            _log.info(f'Could not find grid file {d2}')
                            continue
                        _log.info('Using %s', nc)
                        with xr.open_dataset(nc[0]) as ds:
                            ds = ds.where(
                                    (ds.longitude>-150) & (ds.longitude< -120)
                                     & (ds.latitude>35) & (ds.latitude<60),
                                     drop=True)
                            ds.attrs['deployment_start'] = str(ds.time[0].values)[:16]
                            ds.attrs['deployment_end'] = str(ds.time[-1].values)[:16]
                            # probably should keep track of other stuff, but lets do this for now
                            dsn = xr.Dataset(coords={'time':ds.time},
                                             data_vars={'longitude':('time', ds.longitude.values),
                                                        'latitude':('time', ds.latitude.values)})
                            # make the navdataset:
                            if filen == 0:
                                dsnav = dsn.copy()
                            else:
                                dsnav = xr.concat((dsnav, dsn), dim='time')
                            filen += 1
                            _log.info(f'opened {nc[0]}')
                            att = ds.attrs
                            if ds.longitude.shape[0] > 2:
                                line = np.vstack((ds.longitude, ds.latitude)).T
                                ls = geojson.LineString(line.tolist())
                                feat = geojson.Feature(geometry=ls)
                                for prop in props:
                                    if prop in ds.attrs.keys():
                                        feat.properties[prop] = ds.attrs[prop]
                                    else:
                                        feat.properties[prop] = ''

                                # get URL....
                                feat.properties['url'] = ('' +
                                    'https://cproof.uvic.ca/gliderdata/deployments/' +
                                    d2[2:])
                                # get color:
                                cols = np.random.randint(0, 200, 3)
                                # cols = pygu.get_html_non_blue(colornum)
                                colornum += 1
                                feat.properties['color'] = '#%02X%02X%02X' % (cols[0], cols[1], cols[2])
                                if ds['time'][-1] > np.datetime64(datetime.datetime.now()) - np.timedelta64(2, 'D'):
                                    feat.properties['active'] = True
                                else:
                                    feat.properties['active'] = False

                                features += [feat]

                                # make the kml:
                                pnt = kml.newpoint(coords=[line[-1]])
                                pnt.style.iconstyle.icon.href = 'https://cproof.uvic.ca/deployments/assets/images/slocum_glider.png'
                                coords = []
                                for thelon, thelat  in zip(ds.longitude.values, ds.latitude.values):
                                    coords += [(thelon, thelat)]
                                pnt.timestamp.when = f'{ds.time.values[-1]}'[:-3]
                                ls = kml.newlinestring(coords=coords,
                                    name=att['deployment_name'],
                                    )
                                ls.timespan.begin = f'{ds.time.values[0]}'[:-3]
                                ls.timespan.end = f'{ds.time.values[-1]}'[:-3]
                                ls.style.linestyle.color = 'ee' + '%02X%02X%02X' %  (cols[2], cols[1], cols[0])
                                ls.style.linestyle.width = 3;
                                kml.save(d2[2:]+'/'+att['deployment_name']+'.kml')

                            else:
                                _log.info(f'No good data {d2}')

                    except:
                        _log.info(f'Could not find grid file {d2}')
    dsnav.to_netcdf('CProofTracks.nc')
    feature_collection = geojson.FeatureCollection(features)
    with open(outfile, 'w') as fout:
        s = geojson.dumps(feature_collection)
        fout.write(s)
######
# julia putko addition: add to different geojson file: Glider + ARGO for the C-PROOF website home page
######
    floatsmeds_bgc_floats = argo.bio_prof[argo.bio_prof['file'].str.contains('meds')]
    my_wmos = [4902549,4902550, 4902551, 4902552,4902553,4902554,4902555,4902583,4902584,4902585,4902586,4902587,4902588,4902589,4902612,4902613,4902616]
   # Adding additional ARGO float numbers

    df = pd.DataFrame()
    _log.info('FLoats')
    coords = []
    for f in argo.float(my_wmos):
        _log.info(f)

        for row in f.prof:
           # if f.prof.latitude > 0:  #make not nan
            lat = f.prof.latitude
            lon = f.prof.longitude

            line = np.vstack((f.prof.longitude, f.prof.latitude)).T
              # line = np.vstack((f.prof.longitude, f.prof.latitude)).T
            mask = np.all(np.isnan(line), axis=1)
            line = line[~mask]
            ls = geojson.LineString(line.tolist())
            feat = geojson.Feature(geometry=ls)

            # get color:
            cols = np.random.randint(0, 200, 3)
            colornum += 1
            feat.properties['color'] = '#%02X%02X%02X' % (cols[0], cols[1], cols[2])
            feat.properties['name'] = 'argo'
            feat.properties['active'] = True
        features += [feat]
####

    #jpnote: Addition of separate outfile that contains argo float data as well as glider mission data
    feature_collection = geojson.FeatureCollection(features)
    with open('cproof-deployments_all.geojson', 'w') as fout:
        s = geojson.dumps(feature_collection)
        fout.write(s)


def link_all_grids(dir):
    dir = Path(dir)
    subdirs = list(dir.glob('*'))

    for d in subdirs:
        _log.info(d)
        if not d.is_dir():
            continue
        subdirs2 = list(d.glob('*'))

        for d2 in subdirs2:
            _log.info(d2)
            if not d2.is_dir():
                continue
            ncr = None
            d2 = d2 / Path('L0-gridfiles')
            nc = list(d2.glob('dfo-*grid_corrected.nc'))
            if len(nc) < 1:
                nc = list(d2.glob('dfo-*grid_delayed.nc'))
            if len(nc) < 1:
                nc = list(d2.glob('dfo-*grid.nc'))
            if len(nc) == 1:
                ncr = nc[0]

            if ncr:
                print('name', str(ncr.name))
                target = Path('/Users/cproof/processing/deployments/allgrids/'+str(ncr.name))
                print('Target', target)
                print('Exists', target.exists())
                print('Exists', target.is_symlink())
                if target.is_symlink():
                    target.unlink()
                _log.debug('liniking %s', ncr.resolve())
                target.symlink_to(ncr.resolve())


def make_wget_txt(dir, outfile='mission_all.txt', comment_str=None,
                  prefix='https://cproof.uvic.ca/gliderdata/deployments/'):

    outf = open(outfile, 'w')
    subdirs = glob.glob(dir + '/*')

    for d in subdirs:
        _log.info(d)
        if not os.path.isdir(d):
            continue
        subdirs2 = glob.glob(d + '/*')
        for d2 in subdirs2:
            _log.info(d2)
            if not os.path.isdir(d2):
                continue
            ncr = None
            nc = glob.glob(d2+'/L0-gridfiles/*grid_delayed.nc')
            if len(nc) == 1:
                ncr = nc[0]
            else:
                nc = glob.glob(d2+'/L0-gridfiles/*grid.nc')
                if len(nc) == 1:
                    ncr = nc[0]
            if ncr:
                if comment_str:
                    with xr.open_dataset(ncr) as ds:
                        print(ds.attrs['comment'])
                        if (comment_str in ds.attrs['comment']):
                            outf.write(prefix+ncr)
                            outf.write('\n')
                            _log.info(f'"{comment_str}" in {ncr}')
                        else:
                            _log.info(f'"{comment_str}" not in {ncr}')
                else:
                    outf.write(prefix+ncr)
                    outf.write('\n')

            else:
                _log.info(f'Nothing at {d2}')

    outf.close()

