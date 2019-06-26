import numpy as np
import matplotlib.pyplot as plt
#from pyresample.geometry import  areadef
import pyresample as pyr
from pyresample import kd_tree, geometry


lons = np.linspace(-55, -21, 100)
lats = np.linspace(71, 68, 100)
lonslons, latslats = np.meshgrid(lons, lats)
z = lonslons *latslats *np.sin((lonslons  +55)/34 *2 * np.pi)
#grid_def = GridDefinition(lons=lons, lats=lats)

area_id = 'ease'
description = 'ease grid north'
proj_id = 'WGS 84'
proj_string = 'EPSG:6931'
width = 1000
height = 1000
area_extent = (-2286563, -2336644, -544973, -880249)
area_def = geometry.AreaDefinition(area_id, description, \
                                    proj_id, proj_string, width, 
                                    height, area_extent)\
                                    
x = np.linspace(-2286563, -544973, 100)
y = np.linspace(-2336644, -880249, 100)
xx, yy = np.meshgrid(x,y)
z1 = xx * yy * np.sin((xx + 544973)/1741590 * 2 * np.pi)



#%%

#>>> area_def = geometry.AreaDefinition('areaD', 'Europe (3km, HRV, VTC)', 'areaD',
#...                                {'a': '6378144.0', 'b': '6356759.0',
#...                                 'lat_0': '50.00', 'lat_ts': '50.00',
#...                                 'lon_0': '8.00', 'proj': 'stere'},
#...                                800, 800,
#...                                [-1370912.72, -909968.64,
#...                                 1029087.28, 1490031.36])
#>>> channel1 = np.fromfunction(lambda y, x: y*x, (50, 10))
#>>> channel2 = np.fromfunction(lambda y, x: y*x, (50, 10)) * 2
#>>> channel3 = np.fromfunction(lambda y, x: y*x, (50, 10)) * 3
#>>> data = np.dstack((channel1, channel2, channel3))
#>>> lons = np.fromfunction(lambda y, x: 3 + x, (50, 10))
#>>> lats = np.fromfunction(lambda y, x: 75 - y, (50, 10))
swath_def = geometry.SwathDefinition(lons=lonslons, lats=latslats)


result = kd_tree.resample_nearest(swath_def, z, area_def, \
                                  radius_of_influence=34000, fill_value=np.nan)
