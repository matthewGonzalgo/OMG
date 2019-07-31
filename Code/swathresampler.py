import numpy as np
import os
import glob
import matplotlib.pyplot as plt
from pyresample import kd_tree, geometry 
import xarray as xr
import utm
import time
import argparse # For command line interface
import math
from datetime import datetime
from createdirectories import create_directory

'''
This is a resampler for an entire collection of swath data.
The main argument is the directory containing all the years data for a given swath
Assumes the directory and file names are consistent with the 'greenl' naming convention
Store data array from data file
Parse all annotation files in a given directory for number of longitude and latitude lines, starting longitude and latitude for area extent, and spacing
Keep track of the largest longitude and latitude bounds based on the corners
Use number of longitude and latitude lines to reshape the the data array
Create utm arrays from the latlons
Find minimum data value and make a new data array that fills values lower than minimum with np.nan
Use parsed starting longitude and latitude corners to define area extent of original area
'''

#
# Parses data from anonotation file
# This function assumes every annotation contains information in the same format
# Creates an array of variables and returns a dictionary with proper labelling
#
def get_variables(annotation):
    variables = []
    with open(annotation, "r") as f:
        lines = f.readlines()
        idx = 0
        for line in lines:
            if line.startswith("GRD Latitude Lines"):
                print("Ususally 62: " + str(idx)) 
                startofdata = idx
                break
            idx += 1

        print("Reading lines to get variables:")
        # Assuming all data is on successive lines and there are 14 data points
        for r in range(14):
            ln = lines[startofdata + r].split()
            # The value will be at the 5th index for the first 6 lines and at the 6th index for the last 8
            # There is a weird instance where 'x10' mysteriously appears
            print(ln)
            if r < 6:
                i = 5
            else:
                i = 6
            print(i)
            if '\x10' in ln[i]:
                print('\x10 found in ' + ln[i] + '!')
                ln[i] = ln[i].replace("\x10", '')
            variables.append(ln[i])

    # Convert the list of "strings" into floats
    print("Raw input variables:")
    print(variables)
    print('\n')

    variables = list(map(float, variables)) 

    vardict =	{
    # lat/lon_lines required to be int not float, may not even be necessary though

    'lat_lines' : int(variables[0]),
    'lon_lines' : int(variables[1]),
    'lat_start' : variables[2],
    'lon_start' : variables[3],
    'lat_space' : variables[4],
    'lon_space' : variables[5],
    'ul_lat' : variables[6],
    'ul_lon' : variables[7],
    'ur_lat' : variables[8],
    'ur_lon' : variables[9],
    'll_lat' : variables[10],
    'll_lon' : variables[11],
    'lr_lat' : variables[12],
    'lr_lon' : variables[13],
    }

    print("Variable dictionary created:")
    print(vardict)
    print('\n')
    
    return vardict  


#
# Plots both the original data input and the resampled data
# May need to change saving locations on other machines
#
def save_and_show(original, xspace, yspace, lats, lons, resampled, resolution, name, save_path):
    print("Saving original and resampled plots")
    start_time = time.time()
    plt.figure(figsize=(10,10));plt.imshow(original, vmin=0, vmax=500, cmap='jet');plt.colorbar()
    print("Original plotting time: --- %s seconds ---" % (time.time() - start_time))
    plt.savefig(save_path + name + "original.png")
    print('Saved original')

    start_time = time.time()
    plt.figure(figsize=(10, 10));plt.pcolormesh(lons, lats, original);plt.colorbar();plt.grid()
    print("Original latlon plotting time: --- %s seconds ---" % (time.time() - start_time))

    plt.savefig(save_path + name + "original_latlon.png")
    print('Saved original latlon')

    start_time = time.time()
    plt.figure(figsize=(10, 10));plt.pcolormesh(xspace, yspace, resampled);plt.colorbar();plt.grid()
    print("Resampled plotting time: --- %s seconds ---" % (time.time() - start_time))
    plt.savefig(save_path + name + "resampled.png")
    print('Saved resampled')

    print("Save complete\n")
    
    # Show plots for the purpose of checking
    # plt.show()
    


def make_latlon_edges(lat_start, lon_start, lat_space, lon_space, lat_lines, lon_lines):
    # Area def uses edges
    lat_edges = np.linspace(lat_start - 0.5 * lat_space, (lat_start - 0.5 * lat_space) + (lat_space * (lat_lines + 1)), lat_lines + 1)
    lon_edges = np.linspace(lon_start - 0.5 * lon_space, (lon_start - 0.5 * lon_space) + (lon_space * (lon_lines + 1)), lon_lines + 1)

    return lat_edges, lon_edges

def make_latlon_centers(lat_start, lon_start, lat_space, lon_space, lat_lines, lon_lines):
    # Area def uses edges
    lat_centers = np.linspace(lat_start, lat_start + (lat_space * (lat_lines)), lat_lines)
    lon_centers = np.linspace(lon_start, lon_start + (lon_space * (lon_lines )), lon_lines)

    return lat_centers, lon_centers

def resample(grid, annotation, vardict, lat_centers, lon_centers, area_new, test, resolution):
    print("Grid shape:")
    print(grid.shape)


    num_grid_rows = grid.shape[0]
    num_grid_columns = grid.shape[1]

    #
    # Original Area definition:
    #
    #area_id = 'WGS84'
    #description = 'lat-lon'
    #proj_id = annotation[annotation.find("gr"):annotation.find(".")]
    #proj_string = 'EPSG:4326'
    #width = num_grid_columns
    #height = num_grid_rows 


    # Hack to deal with two cases: west/east on first column 
    #if vardict['ll_lon'] > vardict['ur_lon']:
    #area_extent = (lon_edges_l[0], lat_edges_b[0], lon_edges_r[-1], lat_edges_t[-1])
    #else:   
    #    area_extent = (vardict['ll_lon'], vardict['ur_lat'], vardict['ur_lon'], vardict['ll_lat'])
        
    #area_original = geometry.AreaDefinition(area_id, description, proj_id, proj_string, width, height, area_extent)

    
    #print("area_original shape (must match grid shape):")
    #print(area_original.shape)

    lon_centers_xx, lat_centers_xx = np.meshgrid(lon_centers, lat_centers);
    area_original = geometry.SwathDefinition(lons=lon_centers_xx, lats=lat_centers_xx)
    
    # print("get_lonlats shape:")
    # print(area_new.get_lonlats()[0].shape, area_new.get_lonlats()[1].shape)

    print(datetime.now())
    print("--- Resampling ---")
    wf = lambda r: 1
    start_time = time.time()
    # Multiplying radius of influence by test to account for skipping data when downsizing
    result = kd_tree.resample_custom(area_original, grid, area_new, radius_of_influence=test * math.sqrt(2 * (resolution / 2)**2), fill_value=np.nan, weight_funcs=wf)
    #result = kd_tree.resample_nearest(area_original, grid, area_new, radius_of_influence=test * math.sqrt(2 * (resolution / 2)**2), fill_value=np.nan)

    print("Result calculation time: --- %s seconds ---" % (time.time() - start_time))
    print("Resulting shape:")
    print(result.shape)

    '''
    xe = np.linspace(x_low, x_high, width_new) 
    ye = np.linspace(y_high, y_low, height_new)
    '''
    
    print("Resample complete\n")
    return result

def downsize(grid, lat_edges, lon_edges, lat_centers, lon_centers, test):
    lat_centers = lat_centers[::test]
    lon_centers = lon_centers[::test]
    lat_edges_b = lat_edges[::test]
    lon_edges_l = lon_edges[::test]
    lat_edges_t = lat_edges[1::test]
    lon_edges_r = lon_edges[1::test]
 
    grid = grid[::test, ::test]
    return grid, lat_centers, lon_centers, lat_edges_b, lon_edges_l, lat_edges_t, lon_edges_r

#
# Gets the max and min values in x and y directions for utm 
# as well as the change in both x and y.
# Returns the information as a tuple
#
def get_utm_range(max_lat, max_lon, min_lat, min_lon):
    min_x = utm.from_latlon(min_lat, min_lon)[0]
    max_x = utm.from_latlon(max_lat, max_lon)[0]
    dx = abs(max_x - min_x) # Shouldn't have to be absolute value

    min_y = utm.from_latlon(min_lat, min_lon)[1]
    max_y = utm.from_latlon(max_lat, max_lon)[1]
    dy = abs(max_y - min_y) # Shouldn't have to be absolute value

    return (min_x, max_x, dx, min_y, max_y, dy)


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
def create_utms(resolution, max_lat, max_lon, min_lat, min_lon):
    print("--- Creating utm area ---")

    print('(resolution, max_lat, max_lon, min_lat, min_lon)')
    print(resolution, max_lat, max_lon, min_lat, min_lon)

    (min_x, max_x, dx, min_y, max_y, dy) = get_utm_range(max_lat, max_lon, min_lat, min_lon)
    print('(min_x, max_x, dx, min_y, max_y, dy)')
    print(min_x, max_x, dx, min_y, max_y, dy)

    '''
    x_low = min(ll_x, ur_x)
    x_high = max(ll_x, ur_x)
    y_low = min(ll_y, ur_y)
    y_high = max(ll_y, ur_y)
    '''
    (x_lower, x_higher, y_lower, y_higher) = scale_dimensions(.1, min_x, max_x, min_y, max_y, dx, dy)
    print('--- Scaled dimensions ---')
    print('(x_lower, x_higher, y_lower, y_higher)')
    print(x_lower, x_higher, y_lower, y_higher)

    # Simplifying dimensions 
    real_x_low = round(x_lower, -3)
    real_y_low = round(y_lower, -3)
    print('Simplified dimensions')
    print('Low:')
    print(real_x_low, real_y_low)

    Dx = abs(x_higher - real_x_low)
    Dy = abs(y_higher - real_y_low)


    n_xcells = (math.ceil(Dx / resolution))
    n_ycells = (math.ceil(Dy / resolution))

    real_x_high = real_x_low + n_xcells * resolution
    real_y_high = real_y_low + n_ycells * resolution
    print('High:')
    print(real_x_high, real_y_high)

    print('Real Dx/Dy:')
    print(Dx)
    print(Dy)

    print('Number of x/y cells:')
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

    print('x/y utm corner shapes')
    print(x_utm_corners.shape, y_utm_corners.shape)
    print("utm arrays created\n")
    return (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, n_xcells, n_ycells)


#
# Saves resampled swath data in netCDF format using xarray
# Requires the resampled data, utm arrays, and the new resampled area to save the swath as an xarray Dataset
# The resolution and full path of the data file are used for 
#
def save_resample(resampled, utm_x, utm_y, new_area, resolution, name, save_path):
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
    swath.attrs['no_data'] = np.NaN
    swath.attrs['srid'] = "urn:ogc:def:crs:" + new_area.proj4_string
    swath.attrs['proj4text'] = new_area.proj4_string
    swath.attrs['Projection'] = new_area.description
    swath.attrs['proj_id'] = new_area.proj_id
    swath.attrs['Insitution'] = 'JPL'
    swath.attrs['Mission'] = 'Oceans Melting Greenland'
    swath.attrs['Mission website'] = 'https://omg.jpl.nasa.gov/portal/'
    swath.attrs['DOI'] = '10.5067/OMGEV-ICEGA'
    swath.attrs['Citation'] = 'OMG Mission. 2016. Glacier elevation data from the GLISTIN-A campaigns. Ver. 0.1. OMG SDS, CA, USA. Dataset accessed [YYYY-MM-DD] at http://dx.doi.org/10.5067/OMGEV-ICEGA.'
    swath.attrs['Mission Citation'] = '10.5670/oceanog.2016.100'
    swath.attrs['author'] = 'Matthew Gonzalgo & Forrest Graham'
    swath.attrs['Date created'] = str(datetime.now())
    swath.attrs['Processessing code repository'] = 'https://github.com/matthewGonzalgo/OMG'
    swath.attrs['Original data source URL'] = 'https://uavsar.jpl.nasa.gov/cgi-bin/data.pl'
    swath.attrs['File naming convention document'] = 'https://uavsar.jpl.nasa.gov/science/documents/topsar-format.html'
    swath.attrs['Note on swath ID number'] = 'Started numbers from Cape Farewell and moved counterclockwise around Greenland'
    swath.attrs['nx'] = len(utm_x)-1
    swath.attrs['ny'] = len(utm_y)-1
    swath.attrs['_FillValue'] = np.NaN
    swath.attrs['_CoordinateTransformType'] = "Projection"
    swath.attrs['_CoordinateAxisTypes'] = "GeoX GeoY"

    swath['x'].attrs['units']='meters'
    swath['x'].attrs['long_name'] = 'X coordinate of grid cell center'
    swath['x'].attrs['coverage_content_type'] = 'coordinate'
    swath['x'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['x'].attrs['axis'] = 'X'
    swath['x'].attrs['valid_range'] = (np.min(utm_x), np.max(utm_x))

    swath['y'].attrs['units']='meters'
    swath['y'].attrs['long_name'] = 'Y coordinate of grid cell center'
    swath['y'].attrs['coverage_content_type'] = 'coordinate'
    swath['y'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['y'].attrs['axis'] = 'Y'
    swath['y'].attrs['valid_range'] = (np.min(utm_y), np.max(utm_y))

    swath.attrs['cdm_data_type']  = 'Grid'
    swath.attrs['geospatial_lat_units'] = "degree_north"
    swath.attrs['geospatial_lon_units'] = "degree_east"
    swath.attrs['geospatial_x_units'] = "meters"
    swath.attrs['geospatial_y_units'] = "meters"
    swath.attrs['geospatial_bounds_crs'] = new_area.proj4_string
    swath.attrs['geospatial_x_resolution'] = str(resolution) + " meters"
    swath.attrs['geospatial_y_resolution'] = str(resolution) + " meters"
    
    print(swath)
    
    saved = save_path + name + ".nc"
    swath.to_netcdf(saved)
    print("Swath saved to " + saved + "\n")

def max_area(directory):
    # Values initially set to max and min to set up initial comparison 
    max_area_dict =	{
        'max_lat' : -float('inf'),
        'max_lon' : -float('inf'),
        'min_lat' : float('inf'),
        'min_lon' : float('inf')
        }
    os.chdir(directory)
    vardict_dict = {}
    for annotation in glob.glob("*.ann"):
        vardict = get_variables(annotation)
        vardict_dict.update({annotation : vardict}) # Save vardicts for later going through each grid
        varmax_lat = max(vardict['ur_lat'], vardict['ll_lat'], vardict['lr_lat'], vardict['ul_lat'])
        varmax_lon = max(vardict['ur_lon'], vardict['ll_lon'], vardict['lr_lon'], vardict['ul_lon'])
        varmin_lat = min(vardict['ur_lat'], vardict['ll_lat'], vardict['lr_lat'], vardict['ul_lat'])
        varmin_lon = min(vardict['ur_lon'], vardict['ll_lon'], vardict['lr_lon'], vardict['ul_lon'])

        if varmax_lat > max_area_dict['max_lat']:
            max_area_dict['max_lat'] = varmax_lat

        if varmax_lon > max_area_dict['max_lon']:
            max_area_dict['max_lon'] = varmax_lon

        if varmin_lat < max_area_dict['min_lat']:
            max_area_dict['min_lat'] = varmin_lat

        if varmin_lon < max_area_dict['min_lon']:
            max_area_dict['min_lon'] = varmin_lon

    print("Max corners from all years of swath data:")
    print(max_area_dict)
    print('\n')
    return max_area_dict, vardict_dict


if __name__ == '__main__':
    start_time = time.time()

    parser = argparse.ArgumentParser()
    #
    # Format for adding an argument:
    # parser.add_argument("-x", "--fullname", action="store", help="comment", default=default value, dest="variable_name",  type=datatype, required=T/F)
    #
    parser.add_argument("-i", "--input", action="store", help="Complete path to input directory", dest="dir", type=str, required=True)
    parser.add_argument("-t", "--test", action="store", help="Varible to scale down data and decrease runtime for testing", default=1, dest="test", type=int, required=False)
    parser.add_argument("-r", "--res", action="store", help="Resolution", default=10, dest="res", type=int, required=False)
    parser.add_argument("-x", "--suffix", action="store", help="filename suffix", default="", dest="suf",  type=str, required=False)
    parser.add_argument("-o", "--output", action="store", help="Complete path to output directory", dest="save_dir",  type=str, required=True)
    args = parser.parse_args()
    
    # Indicate test runs in their save name
    if args.test != 1:
        args.suf = "testing" + str(args.test) + args.suf
    #
    # Exception handling
    #
    assert os.path.isdir(args.dir), "Directory does not exist"
    # Checking for proper naming convention
    assert not args.dir.endswith('/'), "Do not end directory with slash"
    dir_name = args.dir[args.dir.rfind('/') + 1:]
    print("dir_name: " + dir_name)
    assert dir_name.find('greenl_') != -1 or dir_name.find('grland_') != -1, "Bad directory name"
    print('Good directory')

    assert args.test >= 1, "test must be greater than or equal to 1"

    assert args.res > 0, "resolution must be greater than 0"


    max_area_dict, vardict_dict = max_area(args.dir)
    
    (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, n_xcells, n_ycells) = create_utms(args.res, \
                                                                                        max_area_dict['max_lat'], \
                                                                                        max_area_dict['max_lon'], \
                                                                                        max_area_dict['min_lat'], \
                                                                                        max_area_dict['min_lon'])


    tag = str(utm.from_latlon(max_area_dict['max_lat'], max_area_dict['min_lon'])[2])
    print("Tag: " + tag)

    # Gets name of directory
    swath_id = args.dir[args.dir.rfind('/') + 1:4]
    swath_directory_name = args.dir[args.dir.rfind('/') + 1:]
    print("Swath directory name: " + swath_directory_name)
    #
    # New Area definition we have defined for the swath (all years):
    #
    area_id_new = 'WGS84 / UTMzone ' + tag + 'N'
    description_new = 'UTM ' + tag + 'N ' + swath_directory_name
    proj_id_new = area_id_new
    proj_string_new = '+proj=utm +zone=' + str(tag) + ' +ellps=WGS84 +datum=WGS84 +units=m +no_defs' 
    width_new = n_xcells
    height_new = n_ycells

    print(width_new)
    print(height_new)

    print('(x_utm_corners[0], y_utm_corners[0], x_utm_corners[-1], y_utm_corners[-1])')
    print(x_utm_corners[0], y_utm_corners[0], x_utm_corners[-1], y_utm_corners[-1])

    area_extent_new = (x_utm_corners[0], y_utm_corners[0], x_utm_corners[-1], y_utm_corners[-1]) # may need to revisit
    area_new = geometry.AreaDefinition(area_id_new, description_new, proj_id_new, proj_string_new, width_new, height_new, area_extent_new)

    print("New area shape:")
    print(area_new.shape)
   

    # For data plotting
    xx, yy = np.meshgrid(x_utm_corners, y_utm_corners)

    for grd in glob.glob('*.grd'):
        print("Swath data grd file name: " + grd)
        # 
        # Some swaths have different numbers in different years because they had different headings,
        # so can't assume the directory must match the file names. Instead assume the files are
        # located in the proper directory.
        #
 
        g = np.fromfile(grd, dtype = '<f4')

        ann = grd[:grd.find(".")] + '.ann'
        vdict = vardict_dict.get(ann)       

        grid = np.reshape(g, (vdict['lat_lines'], vdict['lon_lines']))
        print('Min value, usually -10000.0: ' + str(grid.min())) 
        grid_nan = np.where(grid > grid.min(), grid, np.nan)


        print("Annotation:" + ann)
        n = grd[grd.rfind('/') + 1:grd.find(".")]
        print("Name: " + n)
        name = str(args.res) + "m_" + n + args.suf
        print("Modified name: " + name)

        
    
        lat_edges, lon_edges = make_latlon_edges(vdict['lat_start'], vdict['lon_start'], vdict['lat_space'], vdict['lon_space'], vdict['lat_lines'], vdict['lon_lines'])
        print('lat_edges[0], lon_edges[0], lat_edges[-1], lon_edges[-1]')
        print(lat_edges[0], lon_edges[0], lat_edges[-1], lon_edges[-1])
        
        lat_centers, lon_centers = make_latlon_centers(vdict['lat_start'], vdict['lon_start'], vdict['lat_space'], vdict['lon_space'], vdict['lat_lines'], vdict['lon_lines'])
        #grid_nan, lat_centers, lon_centers = downsize(grid_nan, lat_centers, lat_centers, args.test)
        grid_nan, lat_centers, lon_centers, lat_edges_b, lon_edges_l, lat_edges_t, lon_edges_r = \
            downsize(grid_nan, lat_edges, lon_edges, lat_centers, lon_centers, args.test)
        
        plot = resample(grid_nan, ann, vdict, lat_centers, lon_centers, area_new, args.test, args.res)


        # Creating and saving this new directory within the current directory of data
        fig_dir = args.save_dir + args.dir[args.dir.rfind("/"):] + "_plots"
        print("fig_dir: " + fig_dir)
        create_directory(fig_dir)

        # For the purpose of saving,
        fig_dir = fig_dir + "/"
        print("fig_dir/: " + fig_dir)
        
        save_and_show(grid_nan, xx, yy, lat_centers, lon_centers, plot, args.res, name, fig_dir)

        # Also creating and saving this new directory within the current directory
        netCDF_dir = args.save_dir + args.dir[args.dir.rfind("/"):] + "_netCDF"
        print("netCDF_dir/: " + netCDF_dir)
        create_directory(netCDF_dir)
        
        # For the purpose of saving,
        netCDF_dir = netCDF_dir + "/"
        save_resample(plot, x_utm_centers, y_utm_centers, area_new, args.res, name, netCDF_dir)

# Below is old
# (plot, new_area, x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers) = resample(grid_nan, vars, args.res, args.test, name)

# save_resample(plot, x_utm_centers, y_utm_centers, new_area, args.res, name, args.save)
    print("Total runtime of multiresampler: --- %s seconds ---" % (time.time() - start_time))
    print("Done!\n")
