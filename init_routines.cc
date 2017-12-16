// init_routines.cc
/* This module contains routines for initializing arrays to specific values. */

#include <iostream> // For input/output to command line
#include <complex> // For complex numbers
#include "init_routines.h" // Include module's header for consistency check


double LinInitArray(const double xMin, const double xMax, const int NumPts, double* Array, 
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

float LinInitArray(const float xMin, const float xMax, const int NumPts, float* Array, 
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

void ValInitArray(const int NumPts, double*const Array, const double Value)
{
    // Initializes all values in the array to a single value (zero by default).
    for (int i=0; i<NumPts; ++i)
        Array[i] = Value;
}

void ValInitArray(const int NumPts, float*const Array, const float Value)
{
    // Initializes all values in the array to a single value (zero by default).
    for (int i=0; i<NumPts; ++i)
        Array[i] = Value;
}

void ValInitArray(const int NumPts, std::complex<double>*const Array, 
                  const std::complex<double> Value)
{
    // Initializes all values in the array to a single value (zero by default).
    for (int i=0; i<NumPts; ++i)
        Array[i] = Value;
}

void ValInitArray(const int NumPts, std::complex<float>*const Array, 
                  const std::complex<float> Value)
{
    // Initializes all values in the array to a single value (zero by default).
    for (int i=0; i<NumPts; ++i)
        Array[i] = Value;
}