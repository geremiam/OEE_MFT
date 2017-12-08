// init_routines.cc
/* This module contains routines for initializing arrays to specific values. */

#include <iostream> // For input/output to command line
#include "init_routines.h" // Include module's header for consistency check


double InitArray(const double xMin, const double xMax, const int NumPts, double* Array, 
                 const bool endpoint)
{
    // Assigns an equally spaced, NumPts-point grid, from xMin to xMax, to Array.
    double dx=0;
    
    if (endpoint)
        dx = (xMax-xMin)/(NumPts-1);
    else
        dx = (xMax-xMin)/(NumPts);
    
    for (int i=0; i<NumPts; ++i) Array[i] = xMin + i*dx;
    
    return dx;
}

float InitArray(const float xMin, const float xMax, const int NumPts, float* Array, 
                const bool endpoint)
{
    // Assigns an equally spaced, NumPts-point grid, from xMin to xMax, to Array.
    float dx=0;
    
    if (endpoint)
        dx = (xMax-xMin)/(NumPts-1);
    else
        dx = (xMax-xMin)/(NumPts);
    
    for (int i=0; i<NumPts; ++i) Array[i] = xMin + i*dx;
    
    return dx;
}