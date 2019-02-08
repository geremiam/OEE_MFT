// misc_routines.cc
/*
Description
*/
#include <iostream> // Input/output to command line
#include <cmath>
#include "misc_routines.h" // Include header file for consistency check

void copy_array(const int len, const double*const source, double*const target)
{
    // Copied the elements of "source" into "target".
    for (int i=0; i<len; ++i)
        target[i] = source[i];
}

bool check_bound_array(const int len, const double maxval, const double*const array)
{
    // Checks that each element of "array" is not greater (in absolute value) than maxval.
    if (maxval<=0.)
      std::cout << "WARNING: maxval is not greater than zero. 'false' is always returned." << std::endl;
    
    bool IsSmall = true; // Start it off as true
    
    for (int i=0; i<len; ++i) // If a single element is greater than the boud, it becomes false.
      IsSmall = IsSmall && ( std::abs(array[i]) < maxval );
    
    return IsSmall;
}
