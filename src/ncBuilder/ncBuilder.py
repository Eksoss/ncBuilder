# -*- coding: utf-8 -*-
import os
import datetime as dt
import netCDF4 as nc
import numpy as np


def _create_nc_dimension(nc_file, shape):
    '''
    Creates the dimensions of the file (x, y, z, w)
        time = UNLIMITED ;
        level = UNLIMITED ;
        latitude = shape[0] ;
        longitude = shape[1] ;

    Parameters
    ----------
    nc_file : netCDF4.Dataset
        Dataset file to create the dimensions.
    shape : tuple
        Shape of the grid, shape=(latitude.size, longitude.size).

    Returns
    -------
    : bool
        Returns True if the dimensions are created correctly, else returns
        False.
        
    '''
    
    try:
        nc_file.createDimension("latitude", shape[0])
        nc_file.createDimension("longitude", shape[1])
        nc_file.createDimension("level")
        nc_file.createDimension("time")
        return True
    except:
        return False


def _verify_timeseries(nc_file, _nc_config, _variable_kw):
    '''
    Modifies the input dict to include a chunksizes kwarg.
    
    Parameters
    ----------
    nc_file : netCDF4.Dataset
        Dataset file that will be added the new dimension and its own variable.
    _nc_config : dict
        Dict containing general configuration for creating a new nc.Variable.
    _varaible_kw : dict
        Dict containing specific configuration for creating a new nc.Variable.
        
    '''
    
    if _variable_kw.pop('timeseries', False):
        _sizes = [nc_file[_dim].size for _dim in _nc_config['dims']]
        _chunksizes = _sizes[:1] + [np.gcd(64, i) for i in _sizes[1:]]
        _variable_kw['chunksizes'] = tuple(_chunksizes)


def create_nc_dimension(nc_file, var, size=None, variable_kw={}):
    '''
    Creates a new dimension and its variable associated.
    
    Returns True if dimension creation was successful else returns False.
    
    Parameters
    ----------
    nc_file : netCDF4.Dataset
        Dataset file that will be added the new dimension and its own variable.
    var : str
        Dimension/Variable name.
    size : int, optional
        Size of the dimension.
    variable_kw : dict, optional
        Dictionary containing additional kwargs used on create_nc_variable.
        
    Returns
    -------
    : bool
        Returns a bool representing if the process was successful of not.
    
    '''
    
    try:
        nc_file.createDimension(var, size)
        comp_lvl = variable_kw.pop('comp_lvl', 6)
        zlib = variable_kw.pop('zlib', True)
        fill_value = variable_kw.pop('fill_value', np.nan)
        _ = variable_kw.pop('dims', None)
        create_nc_variable(nc_file,
                           var,
                           comp_lvl,
                           zlib,
                           fill_value,
                           variable_kw,
                           dims=(var, ))
        return True
    except:
        return False


def create_nc_variable(nc_file, var, comp_lvl=6, zlib=True, fill_value=np.nan, variable_kw={}, **kwargs):
    '''
    Creates the new variable needed into the given dataset.

    Parameters
    ----------
    nc_file : netCDF4.Dataset
        Dataset file that will be added a new variable.
    var : str
        Variable name.
    comp_lvl : int, optional
        Compression level.
    zlib : bool, optional
        Determines if data will be compressed using gzip within the netCDF file.
    fill_value : float, optional
        Fills the missing values with the parsed value.
    variable_kw : dict, optional
        Dictionary containing additional kwargs used on createVariable, as described on its documentation.

    Optional Parameters
    -------------------
    dims : tuple
        Representing which dimensions from the nc_file itself will be
        used for the new variable.
    dtype : dtype
        The dtype of the new variable.
    long_name : str
        Long name of the variable.
    standard_name : str
        The standard name of the variable.
    units : str
        Unit of the variable.
    timeseries : bool
        Verifies if variable is preferably a timeseries, so chunksizes is applied.
        
    Notes
    -----
    The createVariable documentation can be found here: https://unidata.github.io/netcdf4-python/#Dataset.createVariable
        
    '''
    
    _verify_timeseries(nc_file, kwargs, variable_kw)
    nc_file.createVariable(var,
                           kwargs.get('dtype', np.float32),
                           kwargs.get('dims', ('time',
                                               'level',
                                               'latitude',
                                               'longitude')),
                           zlib=zlib,
                           complevel=int(comp_lvl),
                           fill_value=fill_value,
                           **variable_kw)
    nc_file.variables[var].long_name = kwargs.get('long_name', var)
    nc_file.variables[var].standard_name = kwargs.get('standard_name', var)
    nc_file.variables[var].units = kwargs.get('units', var)


def update_nc(nc_file, var, data, dims=[slice(None), ]):
    '''
    Updates the variable given, within the slice parsed.

    Parameters
    ----------
    nc_file : netCDF4.Dataset
        A dataset containing the dimensions and variable already available to
        be inputed.
    var : str
        Key that'll be used to input the data.
    data : np.ndarray
        Array that will be input inside the dataset variable, if it's smaller
        than the full shape it's needed to parse the dims variable to get the
        proper slice to input the data. It must have the same dimensions as
        target variable.
    dims : tuple, list, optional
        It's a tuple/list of slices that represent the position which the data must
        be input inside the dataset variable.

    Returns
    -------
    : bool

    '''
    
    nc_file.variables[var][dims] = data[:]
    nc_file.sync()
    return True


def initialize_nc(dtNow, filePath='./', fileSuffix='.nc'):
    return nc.Dataset(_file_path(dtNow, filePath, fileSuffix), 'w')


def open_file(dtNow, filePath='./', fileSuffix='.nc'):
    return nc.Dataset(_file_path(dtNow, filePath, fileSuffix), 'r+')


def _file_path(dtNow, filePath='./', fileSuffix='.nc'):
    if isinstance(dtNow, (dt.datetime, dt.date)):
        return os.path.join(f'{filePath}', dtNow.strftime(f'%Y%m%d{fileSuffix}'))
    elif isinstance(dtNow, str):
        return os.path.join(f'{filePath}', f'{dtNow}{fileSuffix}')
    else:
        raise 'Invalid dtNow type, must be str, dt.datetime or dt.date'
        

def create_nc(nc_file, lat, lon, comp_lvl=6, **kwargs):
    '''
    Create the initial variables and dimensions.

    Prepares the dataset file creating its dimensions, creating parsed variables
    and setting up the basic metadata.

    Parameters
    ----------
    nc_file : netCDF4.Dataset
        The netCDF4.Dataset that will be prepared with dimensions and (optional)
        variables.
    lat : np.ndarray
        Array containing increasing values for latitudes.
    lon : np.ndarray
        Array containing increasing valeus for longitudes.
    comp_lvl : int
        Compression level for the netCDF4 variables.

    Optional Parameters
    -------------------
    time : np.ndarray, list
        A array-like containing the datetimes of the time dimension,
        else it's used the dtNow value for 1-sized dimension.
    level : np.ndarray, list
        A array-like containing the float values of the level dimension,
        else it's used only a 1000.
    vars : dict
        A dictionary containing each variable that will be used
        (those can be added later too), parsing a dictionary with the
        dimensions that will be used such as time, level, latitude or
        longitude. Moreover you can give the dtype, long_name,
        standard_name and units of the variable.
    
    Examples
    --------
    Needs dict of variables to be used:
    {'vars':
        {'temp':
            {'dims': ('time', 'level', 'latitude', 'longitude'),
             'dtype': np.float64,
             'long_name': 'temperature',
             'standard_name': 'temperature in C',
             'units': 'C',
             'variable_kw': {'least_significant_digit': 2,
                             'timeseries': True,
                            },
            },
        },
        
        {'rain':
            {'dims': ('time', 'latitude', 'longitude'),
             'dtype': np.float32,
             'long_name': 'precipitation',
             'standard_name': 'precipitation in mm (kg/m^2)',
             'units': 'mm'
            },
        },
    }

    Notes
    -----
    The time dimension will be composed of floats, with the default delta is
    in hours, but it won't be a problem as a float can represent different
    fractions of time.
    
    '''
    
    is_dimensions = _create_nc_dimension(nc_file, [len(lat), len(lon)])

    if is_dimensions:
        time_arr = kwargs.get('time',
                              np.array([dt.datetime.now().replace(hour=0,
                                                                  minute=0,
                                                                  second=0)]))
        level_arr = kwargs.get('level', np.array([1000., ]))
        header_time = time_arr[0]
        float_time_arr = np.array([(time - header_time).total_seconds() / 3600.
                                   for time in time_arr], dtype=np.float64)

        time_units = header_time.strftime('hours since %Y-%m-%d %H:%M:%S')
        
        # The dimensions and its variables are hard coded for consistency
        create_nc_variable(nc_file,
                           'time',
                           comp_lvl=4,
                           zlib=False,
                           fill_value=None,
                           dtype='f8',
                           dims=('time', ),
                           long_name='Time',
                           standard_name='times',
                           units=time_units)
        
        create_nc_variable(nc_file,
                           'level',
                           comp_lvl=4,
                           zlib=False,
                           fill_value=None,
                           dtype='f8',
                           dims=('level', ),
                           long_name='Level',
                           standard_name='air_pressure',
                           units='hPa')
        
        create_nc_variable(nc_file,
                           'latitude',
                           comp_lvl=4,
                           zlib=False,
                           fill_value=None,
                           dtype='f8',
                           dims=('latitude', ),
                           long_name='Latitude',
                           standard_name='latitude',
                           units='degrees_north')
        
        create_nc_variable(nc_file,
                           'longitude',
                           comp_lvl=4,
                           zlib=False,
                           fill_value=None,
                           dtype='f8',
                           dims=('longitude', ),
                           long_name='Longitude',
                           standard_name='longitude',
                           units='degrees_east')

        # dimension variables are updated with their respective arrays
        update_nc(nc_file, 'time', float_time_arr)
        update_nc(nc_file, 'level', level_arr)
        update_nc(nc_file, 'latitude', lat)
        update_nc(nc_file, 'longitude', lon)
        
        '''
          Those are the original way the dimension variables were created,
        it'll be kept here for conference until 0.0.9
          The project will probably change the way the dimensions are created
        by parsing a 'dimensions' key word for its custom creation
        
        '''
        
        # nc_file.createVariable('time', 'f8', ('time', ))
        # nc_file.variables['time'].long_name = 'Time'
        # nc_file.variables['time'].standard_name = 'times'
        # nc_file.variables['time'].units = time_units
        # nc_file.variables['time'][:] = float_time_arr
        
        # nc_file.createVariable('level', 'f8', ('level', ))
        # nc_file.variables['level'].long_name = 'Level'
        # nc_file.variables['level'].standard_name = 'air_pressure'
        # nc_file.variables['level'].units = 'hPa'
        # nc_file.variables['level'][:] = level_arr

        # nc_file.createVariable('latitude', 'f8', ('latitude', ))
        # nc_file.variables['latitude'].long_name = 'Latitude'
        # nc_file.variables['latitude'].standard_name = 'latitude'
        # nc_file.variables['latitude'].units = 'degrees_north'
        # nc_file.variables['latitude'][:] = lat

        # nc_file.createVariable('longitude', 'f8', ('longitude', ))
        # nc_file.variables['longitude'].long_name = 'Longitude'
        # nc_file.variables['longitude'].standard_name = 'longitude'
        # nc_file.variables['longitude'].units = 'degrees_east'
        # nc_file.variables['longitude'][:] = lon

        for var, var_kw in kwargs.get('vars', {}).items():
            create_nc_variable(nc_file, var, comp_lvl, **var_kw)
        
    nc_file.Conventions = "CF-1.4"
    nc_file.Metadata_Conventions = "Unidata Dataset Discovery v1.0"

    nc_file.HISTORY = "Created by Chico and Zang's mathmagical utility at " + \
        dt.datetime.now().strftime('%Y-%m-%dT%H:%M')
    nc_file.sync()

    return lat, lon
