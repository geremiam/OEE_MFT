# plot_routines.py
""" Various routines for plotting specific data formats. """

import matplotlib.pyplot as plt
import numpy as np

def ColorPlot1(row, col, z, Labels=[], clim=[]):
    """ Makes a colorplot of the 2D array z, with the row variable as the vertical axis 
    and the column variable as the horizontal axis. 'row' and 'col' are the coordinate 
    variables for the array 'z'. The entries of 'Labels' are the labels for the rows, 
    the columns, and z, respectively. """
    
    fig, axes = plt.subplots()
    
    extent = [col[0], col[-1], 
              row[0], row[-1]] # Determines the ranges of the axes
    
    if (Labels != []):
        axes.set_ylabel(Labels[0])
        axes.set_xlabel(Labels[1])
    
    cs = axes.imshow(z, extent=extent, interpolation="nearest", origin="lower", 
                     aspect='auto')
    
    if (clim != []): # Sets colorbar limits
        cs.set_clim(clim)
    cbar = fig.colorbar(cs, ax=axes)
    if (Labels != []):
        cbar.set_label(Labels[2])
    
    plt.show()

def ColorPlot2(row, col, z, Labels=[], clim=[]):
    """ In development... This is similar to the previous routine, but uses pcolormesh 
    instead of imshow. 
    https://ocefpaf.github.io/python4oceanographers/blog/2015/01/05/pcolor/ """
    
    fig, axes = plt.subplots()
    
    if (Labels != []):
        axes.set_xlabel(Labels[1])
        axes.set_ylabel(Labels[0])
    
    axes.pcolormesh(z)
    
    plt.show()
