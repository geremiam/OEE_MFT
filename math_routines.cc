// math_routines.cc
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/

#include <iostream> // Input/output to command line
#include <cmath> // Numerous math functions
#include "math_routines.h" // Include header file for consistency check

bool nF0(double const energy)
{
    /* Zero-temperature Fermi function. Overloaded for floats. */
    bool ans = true;
    if (energy > 0)
        ans = false;
    
    return ans;    
}
bool nF0(float const energy)
{
    /* Zero-temperature Fermi function. Overloaded for doubles. */
    bool ans = true;
    if (energy > 0)
        ans = false;
    
    return ans;    
}
double nF(const double beta, const double energy)
{
    /* Fermi function. Overloaded for floats. */
    return 1./( 1. + std::exp(beta*energy) );
}
float nF(const float beta, const float energy)
{
    /* Fermi function. Overloaded for doubles. */
    return 1./( 1. + std::exp(beta*energy) );
}