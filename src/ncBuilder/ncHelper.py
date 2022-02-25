# -*- coding: utf-8 -*-
import netCDF4 as nc
import numpy as np
import pandas as pd
import datetime as dt

def load_time(nc_variable_time):
    '''
    Returns an array of datetimes.

    Calculates the datetimes given the time variable from a netCDF4 file,
    extracting the delta that the values represent and returns the appropriate
    datetimes that represent the original array.

    Parameters
    ----------
    nc_variable_time : netCDF4.Variable
        A netCDF4 variable that represents the time of the file, it must
        contain the `DELTA since DATE` in the units subkey.
        
    Returns
    -------
    : np.ndarray
        Returns a array of datetimes representing the times given the header
        `hours/minutes/days since 1900-01-01 00:00:00.0` and the float delta
        in the array.
        
    '''
    
    delta, _, *date = nc_variable_time.units.split(' ')
    head_time = pd.to_datetime('T'.join(date)).to_pydatetime()
    return np.array([head_time + dt.timedelta(**{delta: float(i)})\
                     for i in nc_variable_time[:]], dtype=dt.datetime)


def get_lats_lons(nc_something, skip=1, trim=0):
    '''
    Return latitudes and longitudes from the given file.

    Checks if the file is a path or a nc_file and try
    to load its lats and lons.

    Parameters
    ----------
    nc_something : str, nc.Dataset
        File to be opened.
    skip : int
        Stride value to skip when reading the array to reduce dimensionality.
    trim : int
        Value to cut the edges of the array if needed.

    Returns
    -------
    : tuple
        Returns a tuple of np.ndarray which contain the latitudes and the
        longitudes of the given file.
    
    '''
    
    nc_file = _verify_nc(nc_something)
    trim_slice = slice(trim, -trim, skip) if trim > 0 else slice(None, None, skip)
    try:
        # some legacy files contain 'lat' and 'lon' as dimensions
        lats = nc_file.variables['lat'][trim_slice]
        lons = nc_file.variables['lon'][trim_slice]
    except:
        lats = nc_file.variables['latitude'][trim_slice]
        lons = nc_file.variables['longitude'][trim_slice]
    return np.array(lats).astype(np.float32), np.array(lons).astype(np.float32)


def _verify_nc(nc_something):
    '''
    Returns a nc.Dataset.

    Given the object given it'll decide what to do with it if valid.

    Parameters
    ----------
    nc_something : str, nc.Dataset
        File to be opened.

    Returns
    -------
    nc_file : nc.Dataset
        nc.Dataset loaded
    
    '''
    
    if isinstance(nc_something, str):
        nc_file = nc.Dataset(nc_something, 'r')
    elif isinstance(nc_something, nc.Dataset):
        nc_file = nc_something
    else:
        raise TypeError('invalid object, must be str (readable .nc/.nc4 file) or nc.Dataset')
    return nc_file


def get_idx_pos(lat0, lon0, lats, lons):
    '''
    Return the i and j indeces of the nearest point.

    Searchs for the nearest indices given the targets and spaces.

    Parameters
    ----------
    lat0 : float
        Target latitude.
    lon0 : float
        Target longitude.
    lats : np.ndarray, list
        Space of latitudes in degrees.
    lons : np.ndarray, list
        Space of longitudes in degrees.

    Returns
    -------
    minindexX : int
        Index that represents which point of lats is the nearest to lat0.
    minindexY : int
        Index that represents which point of lons is the nearest to lon0.
    
    Notes
    -----
    lats and lons must be in degrees, it's not asserted.
    
    '''
    
    latvals = np.radians(lats[:])
    lonvals = np.radians(lons[:])
    lat0_rad = np.radians(lat0)
    lon0_rad = np.radians(lon0)
    minindexX = np.abs(latvals[:] - lat0_rad).argmin()
    minindexY = np.abs(lonvals[:] - lon0_rad).argmin()
    return minindexX, minindexY
