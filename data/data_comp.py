# data_comp.py
""" Compares the data in two identical datasets.
The following are specified as command-line arguments: (1) the path to the NetCDF 
dataset, (2) the name of the variable to be plotted, and (3) the indices to be plotted 
for each of the variable's dimensions. (To specify which two dimensions to plot against, 
enter '-1' instead of an index.) """

from sys import argv # argv is the list of space-separated arguments from command prompt.
import numpy as np
from scipy.io import netcdf # Scipy's NetCDF library
from collections import OrderedDict

def nc_read(file):
    """ Reads a simple format of NetCDF datasets with n dimensions (each having a 
    coordinate variable) and m (non-coordinate) variables. The dimension names and 
    lengths, the coordinate variables, the variable names and the variables are returned. 
    Note that variables are loaded into memory and returned as arrays, so this is not 
    adequate for very large datasets. """
    np.set_printoptions(linewidth=1000)
    
    with netcdf.netcdf_file(file, 'r') as f: # Open the NetCDF file
        dataset_vars = OrderedDict()
        
        for el in f.variables:
            dataset_vars[el] = f.variables[el][:].copy() # Copy the variable data into memory
        
    return dataset_vars

if __name__ == "__main__":
    np.set_printoptions(threshold=2000)
    # Unpack the argv list
    file1 = argv[1] # Element 1 is dataset name (a string)
    file2 = argv[2] # Element 2 is dataset name (a string)
    varname = argv[3]
    
    # Dictionaries containing all the dataset variables
    dataset1_vars = nc_read(file1)
    dataset2_vars = nc_read(file2)
    
    # Select the variable specified in the command-line arguments
    var1 = dataset1_vars[varname]
    var2 = dataset2_vars[varname]
    
    # Calculate relative and absolute differences
    reldiff = np.abs((var1-var2)/var1)
    absdiff = np.abs((var1-var2))
    
    # Choose relative and absolute tolerances
    reltol = 0.05
    abstol = 1.e-6
    
    big_reldiff = (reldiff>reltol)
    big_absdiff = (absdiff>abstol)
    big_both = np.logical_and(big_reldiff, big_absdiff)
    
    print("Number of entries: {}".format(absdiff.size))
    print("With big relative difference: {}".format(reldiff[big_reldiff].size))
    print("With big absolute difference as well: {}".format(reldiff[big_both].size))
    print("reltol = {}".format(reltol))
    
    print("\nreldiff\t\tvar1\t\tvar2\t\tabsdiff")
    
    # Build an array for comparing the different values
    compare = np.transpose( np.array([reldiff.flatten(), var1.flatten(), var2.flatten(), absdiff.flatten()]) )
    compare = compare[compare[:,0].argsort()] # Sort according to first column
    #print(compare[-1580:, :])
    
    for el in compare:
        if el[0]>reltol and el[3]>abstol:
            print(el)
    
    
    #for idx, val in enumerate(reldiff.flatten()):
    #    if (val>reltol) and (absdiff.flatten()[idx]>abstol):
    #        print( "{}\t{}\t{}\t{}".format(reldiff.flatten()[idx], var1.flatten()[idx], 
    #                                      var2.flatten()[idx], absdiff.flatten()[idx]) )
    
    #print("\nAbsolute difference: \n{}".format(absdiff))
    #print("\nRelative difference: \n{}".format(reldiff))
