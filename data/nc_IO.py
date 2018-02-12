# nc_IO.py
""" Routines for interacting with NetCDF datasets. """

import numpy as np
from scipy.io import netcdf

def nc_read(Filename):
    """ Reads a simple format of NetCDF datasets with n dimensions (each having a 
    coordinate variable) and m (non-coordinate) variables. The dimension names and 
    lengths, the coordinate variables, the variable names and the variables are returned. 
    Note that variables are loaded into memory and returned as arrays, so this is not 
    adequate for very large datasets. """
    np.set_printoptions(linewidth=1000)
    
    with netcdf.netcdf_file(Filename, 'r') as File: # Open the NetCDF file
        # Print some relevant info
        #print("\nAttributes: ", File.GlobalAtt)
        #print("\nDimensions: ", File.dimensions)
        #print("\nVariables: ", File.variables)
        
        # Print out the global attributes, assuming they are named 'GlobalAttributes'
        try:
            print("\nGlobal attributes:", File.GlobalAttributes)
        except AttributeError:
            print("\n\tEXCEPTION: No attribute named 'GlobalAttributes' found.\n")
        
        # Get the dimension names and lengths. Note that the dimensions are held in 
        # ordered dictionaries, so their orders are not lost
        dim_names = []
        dim_lengths = []
        for el in File.dimensions:
            dim_names.append(el)
            dim_lengths.append(File.dimensions[el])
        
        print("dim_names:", dim_names)
        print("dim_lengths:", dim_lengths)
        
        # Note that by default, the opened datasets are not loaded into memory.
        # Copy (i.e. load to memory) the dataset variables as numpy arrays held in lists
        coord_vars = []
        var_names = []
        vars = []
        for el in File.variables:
            if (el in dim_names):
                coord_vars.append(File.variables[el][:].copy())
            else:
                var_names.append(el)
                vars.append(File.variables[el][:].copy())
        
        print("coord_vars:", coord_vars)
        print("var_names:", var_names)
        print("vars:", vars)
        
    return dim_names, dim_lengths, coord_vars, var_names, vars
