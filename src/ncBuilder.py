# -*- coding: utf-8 -*-

import os
import datetime as dt
import netCDF4 as nc
import numpy as np


def create_NC_dimension(nc_file, shape, size=0):
    """
    Cria as dimensoes do arquivo (x, y, z)
        Time = UNLIMITED ; bottom_top = 40 ; lat = 192 ; lon = 192 ;
    """
    try:
        nc_file.createDimension("latitude", shape[0])
        nc_file.createDimension("longitude", shape[1])
        nc_file.createDimension("level")
        nc_file.createDimension("time")
        return True
    except:
        return False


def create_NC_variable(nc_file, var, comp_lvl, **kwargs):
    nc_file.createVariable(var, kwargs.get('dtype', np.float32), kwargs.get('dims', ('time', 'level', 'latitude', 'longitude')),
                                   zlib=True, complevel=int(comp_lvl), fill_value=np.nan)
    nc_file.variables[var].long_name = kwargs.get('long_name', var)
    nc_file.variables[var].standard_name = kwargs.get('standard_name', var)
    nc_file.variables[var].units = kwargs.get('units', var)


def update_NC(nc_file, varName, data, dims=[slice(None, None), ]):
    """
    _ix, _iy = int
    _id =  str
    """
    nc_file.variables[varName][dims] = data[:]
    nc_file.sync()
    return True


def initiNC(dtNow, filePath='NC', fileSuffix='_probdens.nc'):
    return nc.Dataset(os.path.join(f'{filePath}', dtNow.strftime(f'%Y%m%d{fileSuffix}')), 'w')


def openFile(dtNow, filePath='NC', fileSuffix='_probdens.nc'):
    return nc.Dataset(os.path.join(f'{filePath}', dtNow.strftime(f'%Y%m%d{fileSuffix}')), 'r+')


def filePath(dtNow, filePath='NC', fileSuffix='_probdens.nc'):
    return os.path.join(f'{filePath}', dtNow.strftime(f'%Y%m%d{fileSuffix}'))
    

def createNC(nc_file, dtNow, lat, lon, comp_lvl=6, **kwargs):
    '''
    Needs dict of variables to be used:
    {'vars':
        {'temp':
            {'dims': ('time', 'level', 'latitude', 'longitude'),
             'dtype': np.float64,
             'long_name': 'temperature',
             'standard_name': 'temperature in C',
             'units': 'C'
            }
        }
    }
    '''
    """
    Cria as variaveis vazias para serem preenchidas
        Create the initial variables and dimensions
    """

    is_dimensions = create_NC_dimension(nc_file, [len(lat), len(lon)])

    if is_dimensions:
        time_arr = kwargs.get('time', np.array([dtNow.replace(hour=0,  minute=0, second=0)])) # dt.datetime
        level_arr = kwargs.get('level', np.array([1000., ]))
        header_time = time_arr[0]
        float_time_arr = np.array([(time - header_time).total_seconds() / 3600. for time in time_arr], dtype=np.float64) # hours since
        
        nc_file.createVariable('time', 'f8', ('time', ))
        nc_file.variables['time'].long_name = 'Time'
        nc_file.variables['time'].standard_name = 'times'
        nc_file.variables['time'].units = header_time.strftime('hours since %Y-%m-%d %H:%M:%S')
        nc_file.variables['time'][:] = float_time_arr
        
        nc_file.createVariable('level', 'f8', ('level', ))
        nc_file.variables['level'].long_name = 'Level'
        nc_file.variables['level'].standard_name = 'air_pressure'
        nc_file.variables['level'].units = 'hPa'
        nc_file.variables['level'][:] = level_arr

        nc_file.createVariable('latitude', 'f8', ('latitude', ))
        nc_file.variables['latitude'].long_name = 'Latitude'
        nc_file.variables['latitude'].standard_name = 'latitude'
        nc_file.variables['latitude'].units = 'degrees_north'
        nc_file.variables['latitude'][:] = lat

        nc_file.createVariable('longitude', 'f8', ('longitude', ))
        nc_file.variables['longitude'].long_name = 'Longitude'
        nc_file.variables['longitude'].standard_name = 'longitude'
        nc_file.variables['longitude'].units = 'degrees_east'
        nc_file.variables['longitude'][:] = lon

        for key, key_dict in kwargs.get('vars', {}).items():
            create_NC_variable(nc_file, key, comp_lvl, **key_dict)
        
    nc_file.Conventions = "CF-1.4"
    nc_file.Metadata_Conventions = "Unidata Dataset Discovery v1.0"

    nc_file.HISTORY = "Created by Chico and Zang's mathmagical utility at " + \
        dt.datetime.now().strftime('%Y-%m-%dT%H:%M')
    nc_file.sync()

    return(lat, lon)
