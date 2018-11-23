# data_anal.py
""" The following are specified as command-line arguments: (1) the path to the NetCDF 
dataset, (2) the name of the variable to be plotted, and (3) the indices to be plotted 
for each of the variable's dimensions. (To specify which two dimensions to plot against, 
enter '-1' instead of an index.) """

from sys import argv # argv is the list of space-separated arguments from command prompt.
import numpy as np
import matplotlib.pyplot as plt

import plot_routines # Routine for plotting data
import nc_IO # Routine for reading data from a NetCDF file

if __name__ == "__main__":
    
    # Unpack the argv list
    filename = argv[1] # Element 1 is dataset name (a string)
    varname  = argv[2] # Element 2 is variable name (a string)
    dims = argv[3:] # Subsequent elements are slices to plot. Build a list.
    for idx, val in enumerate(dims): # Make the list elements integers
        dims[idx] = int(val)
    
    # Check that exactly two dims were chosen for plotting
    if ( np.sum(np.array(dims)<0) != 2 ):
        print("\tWARNING: NUMBER OF PLOTTING DIMENSIONS SHOULD BE 2")
    
    
    var, var_dims, coord_vars = nc_IO.nc_read(filename, varname) # Get data
    # "var" is the variable
    # "var_dims" is a tuple with the dims on which the var is defined
    # "coord_vars" is a dictionary containing coordinate variables related to these dims (if applicable)
    
    # Check that the right number of dimensions were provided on the command line
    if (len(var_dims) != len(dims)):
        print("\tWARNING: INCORRECT NUMBER OF DIMENSIONS WAS PROVIDED")
    
    
    # Build a tuple to index the variable
    tup = dims # Starts off as a list
    plotting_dims = [] # This list will contain the names of the plotting dimensions
    for idx, val in enumerate(tup):
        if (val<0): # Want whole slice if value is negative
            tup[idx] = slice(None,None,None) # Equivalent to ":"
            plotting_dims.append(var_dims[idx])
    tup = tuple(tup) # Redefine as tuple
    
    
    print("Plotting dimensions: {}".format(plotting_dims))
    
    
    # We plot every var in the dataset against the first two dimensions
    Labels = [plotting_dims[0], plotting_dims[1], varname] # Axis labels
    
    # Define the extents of the axes if they have coordinate variables
    if (plotting_dims[0] in coord_vars):
        X = coord_vars[plotting_dims[0]]
        dx = (X[1] - X[0])/2.
        horiz_extent = (X[0]-dx, X[-1]+dx)
    else:
        horiz_extent = () # Leads to default behaviour
    
    if (plotting_dims[1] in coord_vars):
        Y = coord_vars[plotting_dims[1]]
        dy = (Y[1] - Y[0])/2.
        vert_extent = (Y[0]-dy, Y[-1]+dy)
    else:
        vert_extent = () # Leads to default behaviour
    
    fig, ax = plt.subplots()
    plot_routines.ColorPlot(fig, ax, var[tup], Labels, horiz_extent, vert_extent)
    plt.show()
