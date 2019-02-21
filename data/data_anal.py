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


def parse_args(args):
    """ Parses through the arguments to get the filename, varname, and dims. """
    # Unpack the argv list
    filename = argv[1] # Element 1 is dataset name (a string)
    varname  = argv[2] # Element 2 is variable name (a string)
    dims = argv[3:] # Subsequent elements are slices to plot. Build a list.
    for idx, val in enumerate(dims): # Make the list elements integers
        dims[idx] = int(val)
    
    return filename, varname, dims

def main():
    filename, varname, dims = parse_args(argv)
    
    var, var_dims, coord_vars = nc_IO.nc_read_var(filename, varname) # Get data
    # "var" is the variable; "var_dims" is a tuple with the dims on which the var is defined;
    # "coord_vars" is a dictionary containing coordinate variables related to these dims (if applicable)
    
    # Check that the right number of dimensions were provided on the command line
    assert (len(var_dims) == len(dims)), "Incorrect number of dimensions was provided."
    
    
    # Build a tuple to index the variable
    tup = dims # Starts off as a list
    plotting_dims = [] # This list will contain the names of the plotting dimensions
    for idx, val in enumerate(tup):
        if (val<0): # Want whole slice if value is negative
            tup[idx] = slice(None,None,None) # Equivalent to ":"
            plotting_dims.append(var_dims[idx])
    tup = tuple(tup) # Redefine as tuple
    
    
    print("Plotting dimensions: {}".format(plotting_dims))
    # Check that either one or two dims were chosen for plotting
    assert (len(plotting_dims)==1 or len(plotting_dims)==2), "Number of plotting dimensions should be 1 or 2."
    
    if ( len(plotting_dims) == 2 ): # We plot the var against the two plotting dimensions
        
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
        
        plt.close()
        fig, ax = plt.subplots()
        plot_routines.ColorPlot(fig, ax, var[tup], Labels, horiz_extent, vert_extent)
        plt.show()
    
    elif ( len(plotting_dims) == 1 ): # We plot the var against the single plotting dimension
        
        Labels = [plotting_dims[0], varname] # Axis labels
        
        # The horizontal axis for the plot
        if (plotting_dims[0] in coord_vars): # Get the coordinate variable if possible
            horiz_axis = coord_vars[plotting_dims[0]]
        else: # Otherwise just plot vs. index
            horiz_axis = range(len(var[tup]))
        
        plt.close()
        fig, ax = plt.subplots()
        plot_routines.Plot(ax, horiz_axis, var[tup], Labels)
        plt.show()
    
    return

if __name__ == "__main__":
    main()
