# ncBuilder
netCDF4 tool to help creating new files from scratch

## example
```
from ncBuilder import ncBuilder
import datetime as dt
import numpy as np
import pandas as pd

dtNow = dt.datetime.now()
nc_file = ncBuilder.initiNC(dtNow, '/path/to/file', '_suffix.nc')

times = pd.date_range('2021-01-01', periods=10)
lat = np.arange(-30., 10. + 0.05, 0.05)
lon = np.arange(-60., -30 + 0.05, 0.05)
temp = np.random.rand(times.size, lat.size, lon.size)

ncBuilder.createNC(nc_file, dtNow, lat, lon, time=times)
ncBuilder.create_NC_variable(nc_file,
                             'temp',
                             comp_lvl=6,
                             dims=('time',
                                   'latitude',
                                   'longitude'),
                             dtype=np.float64,
                             units='C',
                             standard_name='temperature in °C')
ncBuilder.update_NC(nc_file, 'temp', temp)
nc_file.close()
```
