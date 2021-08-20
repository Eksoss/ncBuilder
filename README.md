# ncBuilder
netCDF4 tool to help creating new files from scratch

## version 0.0.3
- Modified function names to use snake_style pattern
- Comments added to the main functions
- Example updated to follow the new patterns

## example
```
from ncBuilder import ncBuilder
import datetime as dt
import numpy as np
import pandas as pd

dtNow = dt.datetime.now()
times = pd.date_range('2021-01-01', periods=10).to_pydatetime()
lat = np.arange(-30., 10. + 0.05, 0.05)
lon = np.arange(-60., -30 + 0.05, 0.05)
temp = np.random.rand(times.size, lat.size, lon.size)

nc_file = ncBuilder.initialize_nc(dtNow, '/path/to/file', '_suffix.nc')
# nc_file = ncBuilder.initialize_nc('file_name_here', '/path/to/file', '_suffix.nc')

# creating dimensions
ncBuilder.create_nc(nc_file, lat, lon, time=times)

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
ncBuilder.update_nc(nc_file, 'temp', temp[:, 5:, :-5], dims=(slice(None), slice(5, None), slice(None, -5), ))

# closing the nc_file
nc_file.close()
```
