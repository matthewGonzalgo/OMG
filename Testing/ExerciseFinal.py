import numpy as np
import matplotlib.pyplot as plt
import pyresample as pyr
from pyresample import kd_tree, geometry

# longitude and latitude arrays of original swath
lons = np.linspace(-55, -21, 100)
lats = np.linspace(71, 68, 100)
lonslons, latslats = np.meshgrid(lons, lats)
# Z value is arbitrary in this case, just made it look cool
z = lonslons *latslats *np.sin((lonslons  +55)/34 *2 * np.pi) 
swath_def = geometry.SwathDefinition(lons=lonslons, lats=latslats)

# Area of original swath
'''
area_id = 'WGS 84'
description = 'lat lon system'
proj_id = 'WGS 84'
proj_string = 'EPSG:4326'
width = 1000
height = 1000
area_extent = (-55, 71, -21, 68)
area_def_original = geometry.AreaDefinition(area_id, description, \
                                    proj_id, proj_string, width, 
                                    height, area_extent)\
'''

# Area of new common grid in different coordinate system
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

#xe = xedge                                     
xe = np.linspace(-2286563, -544973, 1001)
ye = np.linspace(-2336644, -880249, 1001)
xx, yy = np.meshgrid(xe,ye)


# result = z values of the new coordinate system
result = kd_tree.resample_nearest(swath_def, z, area_def, \
                                  radius_of_influence=34000, fill_value=np.nan)

