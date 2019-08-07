import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from glob import glob
import argparse
from time import time
from createdirectories import create_directory
import os 

""" 
This script creates 'difference netCDF' files between two consecutive years that
a swath was collected. It requires both the directory containing the netCDF 
files for all years of one swath and the directory to save the differences to.

A subdirectory will be created within the outer outoput directory specificed 
which will contain both the plots as png files in addition to the
netCDF files.

'Newer year' and 'Older year' attributes are created for the difference netCDFs.
"""

#
# Gets the netCDF files for all the differences between consecutive years of a
# given swath and returns them as a list.
#
def get_differences(CDF_dir):
    start_time = time()
    print('--- Getting differences ---')
    G = glob(CDF_dir + '/*nc')
    print('::: Number of netCDFs :::')
    print(len(G))

    datasets = []
    for g in G:
        print('Opening dataset' + g)
        datasets.append(xr.open_dataset(g))
    #print(datasets)
    print('\n')

    years_cdf_dict = {}
    # Maps the year collected to the appropriate dataset
    for d in datasets:
        print('...Checking dataset type... :')
        print(type(d))

        print('...Collecting info...')
        print(d.x)
        print('Date collected: ' + d.attrs['Date collected (day-month-year)'])
        years_cdf_dict.update({d.attrs['Year collected'] : d})

    print(years_cdf_dict)
    print('\n')

    # Years in descending order
    years = sorted(years_cdf_dict.keys(), reverse=True)
    print(years)
    differences = []
    for i in range(len(years) - 1):
        newer_year = years[i]
        older_year = years[i + 1]
        print('--- newer, older ---')
        print(newer_year, older_year)

        newer_cdf = years_cdf_dict.get(newer_year)
        older_cdf = years_cdf_dict.get(older_year)
        #print(newer_cdf, older_cdf)

        difference_cdf = newer_cdf - older_cdf
        
        # Assumes CDF_dir is in the format swathresampler creates
        name = CDF_dir[CDF_dir.rfind('/') + 1:]
        print('Name: ' + name)
        id = name[:name.find('_')]
        print('ID: ' + id)
        
        # Makes the CDF description the arbitrary swath ID
        difference_cdf.elevation.attrs['description'] = id
        difference_cdf.elevation.attrs['units'] = 'meters'

        difference_cdf.attrs['Newer year'] = newer_year
        difference_cdf.attrs['Older year'] = older_year
        print('... new cdf ...')
        print(difference_cdf) 
        differences.append(difference_cdf)

    print("\nDifference calculation time: --- %s seconds ---" % (time() - start_time))
    print('--- Differences collected: ---')
    print(len(differences))
    print('Should be one less than ::: Number of netCDFs :::')
    print(len(G))
    print('\n')

    return differences

#
# Saves both a plot of the difference and the netCDF file
def save_differences(differences, save_dir):
    start_time = time()
    print('--- Saving differences ---')

    #
    # A subdirectory will be created to contain the .png and .nc files.
    # The description for every diff in differences should be the same since
    # they all come from the same swath.
    #
    swathid_dir = save_dir + '/' + differences[0].elevation.attrs['description']
    create_directory(swathid_dir)

    for diff in differences:
        # Just making the plot look nicer by matching the dimensions of the 
        # figure appropriately.
        print(diff.dims)
        if diff.dims['y'] > diff.dims['x']:
            plt.figure(figsize=(10, 20))
        else:
            plt.figure(figsize=(20, 10)) 

        diff.elevation.plot(vmin=-20, vmax=20, cmap='bwr_r')

        year_diff_str = str(diff.attrs['Newer year']) + '-' + str(diff.attrs['Older year'])
        id = diff.elevation.attrs['description']
        plt.title('Swath ID ' + id + ' ' + year_diff_str + ' elevation difference')


        saveas = swathid_dir + '/' + id + '_' + year_diff_str
        
        pngsave = saveas + '.png'
        plt.savefig(saveas + '.png')
        print('Saved plot to ' + pngsave)
        
        cdfsave = saveas + '.nc'
        diff.to_netcdf(cdfsave)
        print('Saved netCDF to ' + cdfsave)

        print('Saved ' + year_diff_str + '\n')
        #plt.close()

    print("Saving runtime: --- %s seconds ---" % (time() - start_time))
    print("--- Saving complete ---\n")



if __name__ == '__main__':
    start_time = time()
    parser = argparse.ArgumentParser()
    #
    # Format for adding an argument:
    # parser.add_argument("-x", "--fullname", action="store", help="comment", default=default value, dest="variable_name",  type=datatype, required=T/F)
    #
    parser.add_argument("-i", "--input", action="store", help="Complete path to input CDF_dir containing netCDF files for all years of a swath", dest="input", type=str, required=True)
    parser.add_argument("-o", "--output", action="store", help="Complete path to output CDF_dir to contain netCDF files for every difference between consecutive years", dest="output", type=str, required=True)
    args = parser.parse_args()

    assert args.input[-1] != '/', "Don't end the directory with '/'"
    assert args.output[-1] != '/', "Don't end the directory with '/'"

    difs = get_differences(args.input)
    save_differences(difs, args.output)
    print("Total runtime of swathdifferencer: --- %s seconds ---" % (time() - start_time))
    print("Done!")