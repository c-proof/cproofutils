"""
waypoint_distance
"""

import numpy as np

import logging

mpernm = 1852.0

logger = logging.getLogger(__name__)

def get_simple_distance(shiplon, shiplat, wplon=None, wplat=None,
                           central_lat=None):
    """
    Calculate distance along a line based on waypoint lats and lons.

    shiplon: array-like of floats
        Array of ship longitudes in degrees

    shiplat: array-like of floats
        Array of ship latitudes in degrees

    wplon : array-like of floats (optional)
        Array of decimal degrees longitude defining the line.  If None
        then distance is simply first difference between ship positions.

    wplat : array-like of floats (optional)
        Array of decimal degrees latitude defining the line. If None
        then distance is simply first difference between ship positions.
        Note both wplon and wplat must be provided if one is provided

    central_lat : float (optional)
        central location for distance calculations.

    This uses the simple calculation of
    dx = (lon1 - lon1) * mpernm * 60 * cos(central_lat).  This is inacurate for
    cruises where latitude varies considerably.
    """
    wplon = np.asarray(wplon)
    wplat = np.asarray(wplat)
    shiplon = np.asarray(shiplon)
    shiplat = np.asarray(shiplat)

    lon0 = wplon[0]
    lat0 = wplat[0]
    if central_lat is None:
        central_lat = lat0
    cos_lat = np.cos(central_lat * np.pi / 180)

    shipx = (shiplon - lon0) * 60 * mpernm * cos_lat
    shipy = (shiplat - lat0) * 60 * mpernm

    alongx =  np.zeros_like(shipy)
    acrossx = np.zeros_like(shipy)

    if wplon is None and wplat is None:
        # ship track defines waypoints, and this is just dumb along-track
        # distance....
        alongx[1:] = np.sqrt(np.diff(shipx)**2 + np.diff(shipy)**2)
        return alongx, acrossx

    wpxs = (wplon - lon0) * 60 * mpernm * cos_lat
    wpys = (wplat - lat0) * 60 * mpernm

    along_s = np.zeros((len(wpxs)-1, len(shipx)))
    across_s = np.zeros((len(wpxs)-1, len(shipx)))
    wlens = np.zeros(len(wpxs))
    for ind in range(len(wpxs)-1):
        x0 = wpxs[ind]
        y0 = wpys[ind]
        x1 = wpxs[ind + 1]
        y1 = wpys[ind + 1]
        wlen = np.sqrt((x1 - x0)**2 + (y1 - y0)**2)
        along = ((shipx - x0) * (x1 - x0) + (shipy - y0) * (y1 - y0))
        along = along / wlen
        across = ((shipx - x0) * (y1 - y0) - (shipy - y0) * (x1 - x0))
        across = across / wlen
        # if along < 0 set across to distance from first
        # wp (except for first section)
        if ind > 0:
            bad = along < 0
            across[bad] = np.sqrt((shipx[bad] - x0)**2 + (shipy[bad] - y0)**2)
            along[ind] = 0
        # if along > total length, set to distance from waypoint,
        # except for last section:
        if ind < len(wpxs) - 1:
            bad = along > wlen
            across[bad] = np.sqrt((shipx[bad] -x1)**2 + (shipy[bad] - y1)**2)
            along[bad] = wlen

        across_s[ind, :] = across
        along_s[ind, :] = along

        wlens[ind+1] = wlen + wlens[ind]

    segment = np.argmin(np.abs(across_s), axis=0)
    logger.debug(segment)

    alongx = along_s[segment, np.arange(len(shipx))] + wlens[segment]
    acrossx = across_s[segment, np.arange(len(shipx))]

    logger.debug(alongx)
    logger.debug(acrossx)


    return alongx, acrossx, segment
