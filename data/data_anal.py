# data_anal.py
""" This module acts as a driver that calls the function which reads the data from a file 
    (which the user must specify as a command-line argument) according to specific 
    formatting and the function that plots the data. """

from sys import argv # argv is the list of space-separated arguments from command prompt.
import numpy as np
from numpy import sqrt, pi

import plot_routines # Routine for plotting data
import nc_IO # Routine for reading data from a NetCDF file


# Use the nc_IO module to import the dataset
Filename = argv[1]
dim_names, dim_lengths, coord_vars, var_names, vars = nc_IO.nc_read(Filename)


AxLabels = [dim_names[0], dim_names[1]]
plot_routines.ColorPlot1(coord_vars[0], coord_vars[1], vars[0], AxLabels=AxLabels)
