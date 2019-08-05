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
import copy as copy

'''
This is a resampler for an entire collection of swath data.
The main argument is the directory containing all the years data for a given swath
Assumes the directory and file names are consistent with the 'greenl' 'grland' naming convention
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

    fig_max_points_on_one_side = 2000


    ## plot original data   

    ### calculate how many indices to skip
    nx_orig = original.shape[0]
    ny_orig = original.shape[1]

    skip_x_orig = np.ceil(nx_orig / fig_max_points_on_one_side)
    skip_y_orig = np.ceil(ny_orig / fig_max_points_on_one_side)
    print('original  skip x y: ', skip_x_orig, skip_y_orig)
    skip_orig = int(np.maximum(skip_x_orig, skip_y_orig))
    print('original skipping: ', skip_orig)

    ### imshow

    plt.figure(figsize=(10,10));
    
    if skip_orig > 1:
        plt.imshow(original[::skip_orig, ::skip_orig],vmin=0, vmax=500, cmap='jet'); \
        plt.title('Original data [skipping every ' + str(skip_orig) + ' cells]')
    else:
        plt.imshow(original, vmin=0, vmax=500, cmap='jet')
        plt.title('Original data')
        
    plt.colorbar()
    print("Original plotting time: --- %s seconds ---" % (time.time() - start_time))
    plt.savefig(save_path + name + "original_array.png")
    print('Saved original')


    ### pcolormesh

    start_time = time.time()
    plt.figure(figsize=(10, 10));
    if skip_orig > 1:
        plt.pcolormesh(lons[::skip_orig], \
                       lats[::skip_orig], \
                   original[::skip_orig, ::skip_orig], vmin=0, vmax=500, cmap='jet');
        plt.title('Original data [skipping every ' + str(skip_orig) + ' cells]')
    else:
        plt.pcolormesh(lons, lats, original, vmin=0, vmax=500, cmap='jet');
        plt.title('Original data')
    
    plt.colorbar();plt.grid()
    print("Original latlon plotting time: --- %s seconds ---" % (time.time() - start_time))

    plt.savefig(save_path + name + "original_mapped.png")
    print('Saved original latlon')

    ## Plot Resampled Data
    ### calculate how many indices to skip
    
    nx_resampled = resampled.shape[0]
    ny_resampled = resampled.shape[1]

    skip_x_resampled = np.ceil(nx_resampled / fig_max_points_on_one_side)
    skip_y_resampled = np.ceil(ny_resampled / fig_max_points_on_one_side)
    print('resampled  skip x y: ', skip_x_resampled, skip_y_resampled)
    skip_resampled = int(np.maximum(skip_x_resampled, skip_y_resampled))
    print('resampled skipping: ', skip_resampled)

    start_time = time.time()

    ### imshow

    plt.figure(figsize=(10, 10));

    if skip_resampled > 1:
        plt.imshow(resampled[::skip_resampled, ::skip_resampled], \
                vmin=0, vmax=500, cmap='jet', origin='lower');
        plt.title('resampled data [skipping every ' + str(skip_resampled) + ' cells]')
    else:
        plt.imshow(resampled, \
            vmin=0, vmax=500, cmap='jet', origin='lower');
        plt.title('resampled data')
        
    plt.colorbar();plt.grid()
    print("Resampled plotting time: --- %s seconds ---" % (time.time() - start_time))
    plt.savefig(save_path + name + "resampled_array.png")
    print('Saved resampled')

    ### pcolormesh
    
    plt.figure(figsize=(10, 10));

    if skip_resampled > 1:
        plt.pcolormesh(xspace[::skip_resampled, ::skip_resampled],\
                       yspace[::skip_resampled, ::skip_resampled],\
                    resampled[::skip_resampled, ::skip_resampled],vmin=0, vmax=500, cmap='jet');
        plt.title('resampled data [skipping every ' + str(skip_resampled) + ' cells]')
    else:
        plt.pcolormesh(xspace, yspace, resampled, vmin=0, vmax=500, cmap='jet');
        plt.title('resampled data')
        
    plt.colorbar();plt.grid()
    print("Resampled plotting time: --- %s seconds ---" % (time.time() - start_time))
    plt.savefig(save_path + name + "resampled_mapped.png")
    print('Saved resampled')

    print("Save complete\n")
    
    

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

    area_id = 'WGS84'
    description = 'lat-lon'
    proj_string = 'proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
    proj_id = 'EPSG:4326'
    width  = np.shape(lon_centers)[0]
    height = np.shape(lat_centers)[0]

    print(lon_centers[0], lat_centers[-1], lon_centers[0], lat_centers[0])

    # define the area_extent as [min_lon, min_lat, max_lon, max_lat]
    
    ## if the last lat value is smaller than the first lat value, then we have to switch
    ## which lat value to use when defining area_extent
    if lat_centers[-1] < lat_centers[0]:
        area_extent = (lon_centers[0], lat_centers[-1], lon_centers[-1], lat_centers[0])
    
    ## if the last lat value is smaller than the first lat value, then we have to switch
    ## which lat value to use when defining area_extent
    else:
        area_extent = (lon_centers[0], lat_centers[0], lon_centers[-1], lat_centers[-1])

    area_original = geometry.AreaDefinition(area_id, description, \
        proj_id, proj_string, width, height, area_extent)

    print(datetime.now())
    print("--- Resampling ---")
    wf = lambda r: 1
    start_time = time.time()
    
    # Multiplying radius of influence by test to account for skipping data when downsizing
    result = kd_tree.resample_custom(area_original, grid, area_new, \
        radius_of_influence=test * math.sqrt(2 * (resolution / 2)**2), \
        fill_value=np.nan, weight_funcs=wf)

    ## if the last lat value is smaller than the first lat value, then we have to flip 
    ## the result in the up/down direction (rows)
    
    if lat_centers[-1] < lat_centers[0]:
        result = np.flipud(result)

    
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
# Function for finding the nth occurence of a string.
# This is used to easily get the year collected from the name in save_resample.
# Could potentially be used elsewhere in the script.
# 
def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+len(needle))
        n -= 1
    return start

#
# Saves resampled swath data in netCDF format using xarray
# Requires the resampled data, utm arrays, and the new resampled area to save the swath as an xarray Dataset
# The resolution and full path of the data file are used for 
#
def save_resample(resampled, utm_x, utm_y, new_area, resolution, name, save_path):
    new_grid_lon, new_grid_lat = new_area.get_lonlats()

    print ('NGLAT 0/-1 ')
    print (new_grid_lat[0,0], new_grid_lat[-1,0])
    print ('NGLON 0/-1 ')
    print (new_grid_lon[0,0], new_grid_lon[-1,0])

    # another manifestation of the orientation problem
    if new_grid_lat[0,0] > new_grid_lat[-1,0]:
        new_grid_lat = np.flipud(new_grid_lat)

    #da = xr.DataArray(g, dims= ('lat','lon'), coords={'lat': lat, 'lon': lon})
    swath = xr.DataArray(resampled, \
                dims = ('y','x'), \
                coords = {'y': utm_y, \
                          'x': utm_x, \
                          'lon': (('y','x'), new_grid_lon), \
                          'lat': (('y','x'), new_grid_lat)})

    swath.name = 'elevation'
    swath = swath.to_dataset()


    print ('SAVING .... ')
    print ('new area : ')
    print (new_area)

    # 
    # Based on the naming convention we defined, the date that the 
    # swath was collected comes after the 5th underscore and ends
    # before the 6th underscore.
    # The year as well as the specific date will be saved as a swath attributes.
    #
    date_collected = name[find_nth(name, '_', 5) + 1:find_nth(name, '_', 6)]
    print('data collected : ')
    print(date_collected)
    year_collected = date_collected[0:2]
    day_collected = date_collected[2:4]
    month_collected = date_collected[4:6]

    swath.elevation.attrs['description'] = name
    swath.elevation.attrs['units'] = 'meters'

    swath.attrs['xmin'] = int(np.min(utm_x))
    swath.attrs['ymax'] = int(np.max(utm_y))
    swath.attrs['spacing'] = int(resolution)
    #swath.attrs['grid_mapping'] = 'crs'    
    swath.attrs['no_data'] = np.NaN
    #swath.attrs['srid'] = "urn:ogc:def:crs:" + new_area.proj4_string
    swath.attrs['proj4text'] = new_area.proj4_string
    swath.attrs['proj4string'] = new_area.proj4_string
    swath.attrs['Projection'] = new_area.description
    swath.attrs['proj4'] = new_area.proj_id
    swath.attrs['Year collected'] = 2000 + int(year_collected)
    swath.attrs['Date collected (day/month/year)'] = day_collected + '/' + month_collected + '/' + year_collected
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

    swath['x'].attrs['units']='meter'
    swath['x'].attrs['long_name'] = "Cartesian x-coordinate";
    swath['x'].attrs['coverage_content_type'] = 'coordinate'
    swath['x'].attrs['standard_name'] = 'projection_x_coordinate'
    swath['x'].attrs['axis'] = 'X'
    swath['x'].attrs['valid_range'] = (np.min(utm_x), np.max(utm_x))

    swath['y'].attrs['units']='meter'
    swath['y'].attrs['long_name'] =  "Cartesian y-coordinate";
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



def remove_outliers(data, block_size, z_max):
    nr = data.shape[0]
    nc = data.shape[1]    
    
    num_blocks_c = int(np.ceil(nc / block_size))
    num_blocks_r = int(np.ceil(nr / block_size))
    
    print(nr, nc, num_blocks_r, num_blocks_c)
    
    data_new = np.zeros(data.shape)*np.nan
    
    for bc in range(num_blocks_c):
        startc = (bc)*block_size
        endc   = (bc+1  )*block_size
        endc   = np.minimum(endc, nc)
        
        for br in range(num_blocks_r):
            startr = (br)*block_size
            endr   = (br+1  )*block_size
            endr = np.minimum(endr, nr)
            
            #print (bc, br, startc, endc, startr, endr)

            tmp_data = data[startr:endr, startc:endc]
            num_not_nan = np.count_nonzero(~np.isnan(tmp_data.ravel()))
            
            # require at least 20 not nan points here.
            if num_not_nan > 20:
                tmp_mean = np.nanmean(tmp_data)
                tmp_std  = np.nanstd(tmp_data)

                if np.abs(tmp_mean) > 0 and tmp_std > 0:
                    tmp_z    = (tmp_data - tmp_mean) / tmp_std
                
                    tmp_data_clean = np.where(np.abs(tmp_z) < z_max, tmp_data, np.nan)
                    data_new[startr:endr, startc:endc] = tmp_data_clean
                else:
                    data_new[startr:endr, startc:endc] = tmp_data
            else:
                data_new[startr:endr, startc:endc] = tmp_data

    return data_new


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
    
    (x_utm_corners, y_utm_corners, x_utm_centers, y_utm_centers, \
        n_xcells, n_ycells) = create_utms(args.res, \
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
    proj_id_new = "+init=epsg:326" + tag;
    proj_string_new = '+proj=utm +zone=' + str(tag) + ' +ellps=WGS84 +datum=WGS84 +units=m +no_defs' 
    width_new = n_xcells
    height_new = n_ycells

    area_extent_new = (x_utm_corners[0], y_utm_corners[0], \
        x_utm_corners[-1], y_utm_corners[-1]) # may need to revisit

    area_new = geometry.AreaDefinition(area_id_new, description_new, \
        proj_id_new, proj_string_new, width_new, height_new, area_extent_new)

    print("New area:")
    print(area_new)
    print("New area shape:")
    print(area_new.shape)
   

    print(width_new)
    print(height_new)

    print('(x_utm_corners[0], y_utm_corners[0], x_utm_corners[-1], y_utm_corners[-1])')
    print(x_utm_corners[0], y_utm_corners[0], x_utm_corners[-1], y_utm_corners[-1])


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

        # first remove points that are known bad data.  Here, anything less than 1000m
        # is definitely bad.

        bad_data_flag = -1000

        print('Min value, usually -10000.0: ' + str(grid.min())) 
        grid_nan = np.where(grid > bad_data_flag, grid, np.nan)

        # now scan through the data are remove outliers.
        # pull out blocks of block_size x block_size, and remove points that are 
        # more than z_max standard deviations from the mean within that block
        block_size = 300
        z_max = 5
        print ('... removing outliers ', str(block_size), str(z_max))
        grid_remove_outliers = remove_outliers(grid_nan, block_size, z_max)

        grid = grid_remove_outliers
        
        print("Annotation:" + ann)
        n = grd[grd.rfind('/') + 1:grd.find(".")]
        print("Name: " + n)
        name = str(args.res) + "m_" + n + args.suf
        print("Modified name: " + name)

        
    
        lat_edges, lon_edges = make_latlon_edges(vdict['lat_start'], \
            vdict['lon_start'], vdict['lat_space'], vdict['lon_space'], \
            vdict['lat_lines'], vdict['lon_lines'])

        print('lat_edges[0], lon_edges[0], lat_edges[-1], lon_edges[-1]')
        print(lat_edges[0], lon_edges[0], lat_edges[-1], lon_edges[-1])
        
        lat_centers, lon_centers = make_latlon_centers(vdict['lat_start'], \
            vdict['lon_start'], vdict['lat_space'], vdict['lon_space'], \
            vdict['lat_lines'], vdict['lon_lines'])

        grid, lat_centers, lon_centers, lat_edges_b, lon_edges_l, lat_edges_t, lon_edges_r = \
            downsize(grid, lat_edges, lon_edges, lat_centers, lon_centers, args.test)
        



        print('lat_centers[0], lat_centers[-1], lon_centers[0], lon_centers[-1]')
        print(lat_centers[0], lat_centers[-1], lon_centers[0], lon_centers[-1])

        plot = resample(grid, ann, vdict, lat_centers, lon_centers,\
                         area_new, args.test, args.res)


        # Creating and saving this new directory within the current directory of data
        create_directory(args.save_dir)

        fig_dir = args.save_dir + args.dir[args.dir.rfind("/"):] + "_plots"
        print("fig_dir: " + fig_dir)
        create_directory(fig_dir)

        # For the purpose of saving,
        fig_dir = fig_dir + "/"
        print("fig_dir/: " + fig_dir)
        
        save_and_show(grid, xx, yy, lat_centers, lon_centers, plot, args.res, name, fig_dir)

        # Also creating and saving this new directory within the current directory
        netCDF_dir = args.save_dir + args.dir[args.dir.rfind("/"):] + "_netCDF"
        print("netCDF_dir/: " + netCDF_dir)
        create_directory(netCDF_dir)
        
        # For the purpose of saving,
        netCDF_dir = netCDF_dir + "/"
        save_resample(plot, x_utm_centers, y_utm_centers, area_new, args.res, name, netCDF_dir)

# save_resample(plot, x_utm_centers, y_utm_centers, new_area, args.res, name, args.save)
    print("Total runtime of multiresampler: --- %s seconds ---" % (time.time() - start_time))
    print("Done!\n")
