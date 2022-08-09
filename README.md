# ncBuilder
netCDF4 tool to help creating new files from scratch, and loading variables from existing ones.



## example
```
from ncBuilder import (
    ncBuilder,
    ncHelper,
)
import datetime as dt
import numpy as np
import pandas as pd

dtNow = dt.datetime.now()
times = pd.date_range('2021-01-01', periods=10).to_pydatetime()
lat = np.arange(-30., 10. + 0.05, 0.05)
lon = np.arange(-60., -30 + 0.05, 0.05)
temp = np.random.rand(times.size, lat.size, lon.size)

nc_file = ncBuilder.initialize_nc(dtNow, '/path/to/file', '_suffix.nc')
# or
nc_file = nc.Dataset('your_nc_file.nc', 'w')
# nc_file = ncBuilder.initialize_nc('file_name_here', '/path/to/file', '_suffix.nc')

# creating dimensions
ncBuilder.create_nc(nc_file,
                    lat,
                    lon,
                    time=times)
# or create variables while creating dimensions
ncBuilder.create_nc(nc_file,
                    lat,
                    lon,
                    time=times, 
                    vars={'temp': {'dims': ('time', 'latitude', 'longitude'),
                                   'units': 'ºC',
                                   'standard_name': 'temperature in ºC',
                                   },
                          })

# creating a single variable
ncBuilder.create_nc_variable(nc_file,
                             'temp',
                             comp_lvl=6,
                             dims=('time',
                                   'latitude',
                                   'longitude'),
                             dtype=np.float64,
                             units='°C',
                             standard_name='temperature in °C')

# updating a variable
ncBuilder.update_nc(nc_file, 'temp', temp)

# updating partially a variable
ncBuilder.update_nc(nc_file, 'temp', temp[:, 5:, :-5], dims=(slice(None), 
                                                             slice(5, None),
                                                             slice(None, -5)))

# loading time
times = ncHelper.load_time(nc_file.variables['time'])

# getting lats and lons
lats, lons = ncHelper.get_lats_lons(nc_file)

# getting nearest point
idx, jdx = ncHelper.get_idx_pos(lat0, lon0, lats, lons)

# closing the nc_file
nc_file.close()
```

## version 0.0.7a
- Packages are now installed with the lib, they were not configured properly.

## version 0.0.7
- Added create_nc_dimension, so new custom dimensions can be easily added.
- Minor modifications to the inner creation of dimensions to accomodate the new function.

## version 0.0.6
- Added zlib, fill_value and variable_kw when creating a new variable, so it's more customizable
- Changed the way dimension variables are created for more consistency of the process

## version 0.0.5
- Added ncHelper, which contains load_time, get_lats_lons and get_idx_pos

## version 0.0.3
- Modified function names to use snake_style pattern
- Comments added to the main functions
- Example updated to follow the new patterns
