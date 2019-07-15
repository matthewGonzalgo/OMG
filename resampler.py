import numpy as np
import os
# import glob
import matplotlib.pyplot as plt
from pyresample import kd_tree, geometry 
import xarray as xr
import utm
import time
import argparse # For command line interface
import math

'''
Take data and annotation files as arguments
Store data array from data file
Parse annotation file for number of longitude and latitude lines, starting longitude and latitude for area extent, and spacing
Use number of longitude and latitude lines to reshape the the data array
Find minimum data value and make a new data array that fills values lower than minimum with np.nan
Use parsed starting longitude and latitude corners to define area extent of original area
'''

#
# Parses data from anonotation file
# This function assumes every annotation contains information in the same format
#
def get_variables(annotation):
    variables = []
    with open(annotation, "r") as f:
        lines = f.readlines()
        idx = 0
        for line in lines:
            if line.startswith("GRD Latitude Lines"): # Start of data should start at index 62 but checking just in case
                print(idx) 
                startofdata = idx
                break
            idx += 1

        # Assuming all data is on successive lines and there are 14 data points
        for r in range(14):
            ln = lines[startofdata + r].split()
            
            # The value will be at the 5th index for the first 6 lines and at the 6th index for the last 8
            if r < 6:
                variables.append(ln[5]) 
            else:
                variables.append(ln[6])

    # Convert the list of "strings" into floats
    variables = list(map(float, variables)) 
    return variables  

#
# Plots both the original data input and the resampled data
#
def plot_data(original, xspace, yspace, resampled):
    start_time = time.time()
    plt.figure(figsize=(10,10));plt.imshow(original, vmin=0, vmax=500, cmap='jet');plt.colorbar()
    print("Original plotting time: --- %s seconds ---" % (time.time() - start_time))

    start_time = time.time()
    plt.figure(figsize=(10, 10));plt.pcolormesh(xspace, yspace, resampled);plt.colorbar();plt.grid()
    print("Resampled plotting time: --- %s seconds ---" % (time.time() - start_time))

    plt.show()
    print("Plotting complete")



def resample(grid, variables, resolution, test):
    # Assuming data and annotation file have the same name, just different type and contain only one period

    # lat/lon_lines required to be int not float
    lat_lines = int(variables[0])
    lon_lines = int(variables[1])
    print(lat_lines)
    print(lon_lines)
    lat_start = variables[2]
    lon_start = variables[3]
    lat_space = variables[4]
    lon_space = variables[5]
    # ul_lat = variables[6]
    # ul_lon = variables[7]
    ur_lat = variables[8]
    ur_lon = variables[9]
    ll_lat = variables[10]
    ll_lon = variables[11]
    # lr_lat = variables[12]
    # lr_lon = variables[13]
    
    (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, n_xcells, n_ycells) = \
        create_utms(resolution, ll_lat, ll_lon, ur_lat, ur_lon)
    
    lats = np.linspace(lat_start - 0.5 * lat_space, (lat_start - 0.5 * lat_space) + (lat_space * (lat_lines + 1)), lat_lines + 1)
    lons = np.linspace(lon_start - 0.5 * lon_space, (lon_start - 0.5 * lon_space) + (lon_space * (lon_lines + 1)), lon_lines + 1)

    # Option to downside for testing purposes
    grid = grid[::test, ::test]
    lats = lats[::test]
    lons = lons[::test]
    if test is not 1:
        lon_lines = math.ceil(lon_lines / test)
        lat_lines = math.ceil(lat_lines / test)
        '''ur_lat = math.ceil(ur_lat / test)
        ur_lon = math.ceil(ur_lon / test)
        ll_lat = math.ceil(ll_lat / test)
        ll_lon = math.ceil(ll_lon / test)'''

    #
    # Original Area definition:
    #
    area_id = 'WGS84'
    description = 'lat-lon'
    proj_id = annotation[annotation.find("gr"):annotation.find(".")]
    proj_string = 'EPSG:4326'
    width = lon_lines
    height = lat_lines
    area_extent = (ll_lon, ll_lat, ur_lon, ur_lat)
    area_original = geometry.AreaDefinition(area_id, description, proj_id, proj_string, width, height, area_extent)

    print(area_original.shape)

    # (x_low, x_high, dx, y_low, y_high, dy) = get_utm_range(ll_lat, ll_lon, ur_lat, ur_lon)
    # (x_lower, x_higher, y_lower, y_higher) = scale_dimensions(.2, x_low, x_high, y_low, y_high, dx, dy)
    tag = str(utm.from_latlon(ll_lat, ll_lon)[2])
    print("Tag: " + tag)
    #
    # New Area definition we have defined:
    #
    area_id_new = area_id
    description_new = 'UTM ' + tag + "N"
    proj_id_new = proj_id + '_new'
    proj_string_new = 'EPSG:326' + tag
    width_new = n_xcells
    height_new = n_ycells

    print(width_new)
    print(height_new)

    area_extent_new = (x_utm_corners[0], y_utm_corners[-1], x_utm_centers[-1], y_utm_corners[0])
    area_new = geometry.AreaDefinition(area_id_new, description_new, proj_id_new, proj_string_new, width_new, height_new, area_extent_new)

    
    print("New area shape:")
    print(area_new.shape)
    # print("get_lonlats:")
    # print(area_new.get_lonlats())
    print("get_lonlats shape:")
    print(area_new.get_lonlats()[0].shape, area_new.get_lonlats()[1].shape)

    print("area_original shape:")
    print(area_original.shape)
    print("Grid shape:")
    print(grid.shape)

    #wf = lambda r: 1 - r/200
    start_time = time.time()
    #result = kd_tree.resample_custom(area_original, grid, area_new, radius_of_influence=resolution / 2, fill_value=np.nan, weight_funcs=wf)
    result = kd_tree.resample_nearest(area_original, grid, area_new, radius_of_influence=resolution * 10, fill_value=np.nan)

    print("Result calculation time: --- %s seconds ---" % (time.time() - start_time))
    print("Result shape:")
    print(result.shape)

    '''
    xe = np.linspace(x_low, x_high, width_new) 
    ye = np.linspace(y_high, y_low, height_new)
    '''
    xx, yy = np.meshgrid(x_utm_corners, y_utm_corners)
    

    
    print("Resample complete")

    
    plot_data(grid, xx, yy, result)

    
    return (result, area_new, x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers)


#
# Gets the max and min values in x and y directions for utm 
# as well as the change in both x and y.
# Returns the information as a tuple
#
def get_utm_range(ll_lat, ll_lon, ur_lat, ur_lon):
    x_low = utm.from_latlon(ll_lat, ll_lon)[0]
    x_high = utm.from_latlon(ur_lat, ur_lon)[0]
    dx = x_high - x_low

    y_low = utm.from_latlon(ll_lat, ll_lon)[1]
    y_high = utm.from_latlon(ur_lat, ur_lon)[1]
    dy = y_high - y_low

    return (x_low, x_high, dx, y_low, y_high, dy)

# 
# Increases an area in all directions by a set factor given the area's corners
# Returns the new corners as a tuple
#
def scale_dimensions(factor, x_low, x_high, y_low, y_high, dx, dy):
    x_lower = x_low - factor * dx
    x_higher = x_high + factor * dx
    y_lower = y_low - factor * dy
    y_higher = y_high + factor * dy

    return (x_lower, x_higher, y_lower, y_higher)


# 
# Creates utm coordinate arrays to use for saving with the xarray Dataset
#
def create_utms(resolution, ll_lat, ll_lon, ur_lat, ur_lon):
 
    (x_low, x_high, dx, y_low, y_high, dy) = get_utm_range(ll_lat, ll_lon, ur_lat, ur_lon)
    (x_lower, x_higher, y_lower, y_higher) = scale_dimensions(.2, x_low, x_high, y_low, y_high, dx, dy)
    
    # Simplifying dimensions 
    real_x_low = round(x_lower, -3)
    real_y_low = round(y_lower, -3)


    Dx = x_higher - real_x_low
    Dy = y_higher - real_y_low


    n_xcells = (math.ceil(Dx / resolution))
    n_ycells = (math.ceil(Dy / resolution))

    real_x_high = real_x_low + n_xcells * resolution
    real_y_high = real_y_low + n_ycells * resolution
    
    print(dx)
    print(dy)
    print(Dx)
    print(Dy)
    print(n_xcells)
    print(n_ycells)

    n_xcorners = n_xcells + 1
    n_ycorners = n_ycells + 1

    x_utm_corners = np.linspace(real_x_low, real_x_high, n_xcorners)
    y_utm_corners = np.linspace(real_y_low, real_y_high, n_ycorners)

    x_utm_centers = np.linspace(real_x_low + .5 * resolution, \
    real_x_high - .5 * resolution, n_xcells)
    y_utm_centers = np.linspace(real_y_low + .5 * resolution, \
    real_y_high - .5 * resolution, n_ycells)
    # x_utms = np.linspace(real_x_low, real_x_high, 76)
    # y_utms = np.linspace(real_y_low, real_y_high, 258)

    print(x_utm_corners.shape, y_utm_corners.shape)
    print("utm arrays created")
    return (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, n_xcells, n_ycells)

#
# Saves resampled swath data in netCDF format using xarray
# Requires the resampled data, utm arrays, and the new resampled area to save the swath as an xarray Dataset
# The resolution and full path of the data file are used for 
#
def save_resample(resampled, utm_x, utm_y, new_area, resolution, name, path):
    new_grid_lon, new_grid_lat = new_area.get_lonlats()
    swath = xr.Dataset( \
    {'elevation': (['y','x'], resampled)}, \
    coords = {'y': utm_y, \
    'x': utm_x, \
    'longitude': (('y','x'), new_grid_lon), \
    'latitude' : (('y','x'), new_grid_lat)})

    

    swath.elevation.attrs['description'] = name
    swath.elevation.attrs['units'] = 'meters'

    swath.attrs['xmin'] = np.min(utm_x)
    swath.attrs['ymax'] = np.max(utm_y)
    swath.attrs['dx_spacing'] = resolution
    swath.attrs['dy_spacing'] = resolution
    swath.attrs['grid_mapping'] = 'crs'    
    swath.attrs['no_data'] = -9999.0
    swath.attrs['srid'] = "urn:ogc:def:crs:EPSG::32624"
    swath.attrs['proj4text'] = new_area.proj4_string
    swath.attrs['Projection'] = new_area.description
    swath.attrs['proj_id'] = new_area.proj_id
    swath.attrs['Insitution'] = 'JPL'
    swath.attrs['author'] = 'Matthew G.'
    swath.attrs['nx'] = len(utm_x)-1
    swath.attrs['ny'] = len(utm_y)-1
    swath.attrs['_FillValue'] = np.NaN
    swath.attrs['_CoordinateTransformType'] = "Projection"
    swath.attrs['_CoordinateAxisTypes'] = "GeoX GeoY"

    swath['x'].attrs['units']='meters'
    swath['x'].attrs['long_name'] = 'X'
    swath['x'].attrs['coverage_content_type'] = 'coordinate'
    swath['x'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['x'].attrs['axis'] = 'X'
    swath['x'].attrs['valid_range'] = (np.min(utm_x), np.max(utm_x))

    swath['y'].attrs['units']='meters'
    swath['y'].attrs['long_name'] = 'Y'
    swath['y'].attrs['coverage_content_type'] = 'coordinate'
    swath['y'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['y'].attrs['axis'] = 'Y'
    swath['y'].attrs['valid_range'] = (np.min(utm_y), np.max(utm_y))

    swath.attrs['cdm_data_type']  = 'Grid'
    swath.attrs['geospatial_lat_units'] = "degree_north"
    swath.attrs['geospatial_lon_units'] = "degree_east"
    swath.attrs['geospatial_x_units'] = "meters"
    swath.attrs['geospatial_y_units'] = "meters"
    swath.attrs['geospatial_bounds_crs'] = "EPSG:32624"
    swath.attrs['geospatial_x_resolution'] = str(resolution) + " meters"
    swath.attrs['geospatial_y_resolution'] = str(resolution) + " meters"
    
    print(swath)
    swath.to_netcdf(path + name + ".nc")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #
    # Format for adding an argument:
    # parser.add_argument("-x", "--fullname", action="store", help="comment", default=default value, dest="variable_name",  type=datatype, required=T/F)
    #
    parser.add_argument("-d", "--data", action="store", help="Complete path to data file", dest="data", type=str, required=True)
    parser.add_argument("-t", "--test", action="store", help="Varible to scale down data and decrease runtime for testing", default=1, dest="test", type=int, required=False)
    parser.add_argument("-r", "--res", action="store", help="Resolution", default=100, dest="res", type=int, required=False)
    parser.add_argument("-s", "--save", action="store", help="Path to save", default="/Users/gonzalgo/Desktop/", dest="save",  type=str, required=False)

    args = parser.parse_args()
    name = args.data[args.data.find("greenl"):args.data.find(".")]
    print("Name: " + name)
    annotation = name + '.ann'
    print("Annotation:" + annotation)
    vars = get_variables(annotation)
    print(vars)

    lat_lines = int(vars[0])
    lon_lines = int(vars[1])
    print(lat_lines)
    print(lon_lines)

    g = np.fromfile(args.data, dtype = '<f4')
    grid = np.reshape(g, (lat_lines, lon_lines))
    print(grid.min()) # Should be -10000.0
    grid_nan = np.where(grid > grid.min(), grid, np.nan)


    (plot, new_area, x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers) = resample(grid_nan, vars, args.res, args.test)
    name = args.data[args.data.find("greenl"):args.data.find(".")]
    save_resample(plot, x_utm_centers, y_utm_centers, new_area, args.res, name, args.save)

    print("Done")