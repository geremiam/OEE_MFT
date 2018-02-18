// math_routines.cc
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/

#include <iostream> // Input/output to command line
#include <cmath> // Numerous math functions
#include <algorithm> // Useful algorithms. Here, we use std::sort()
#include "math_routines.h" // Include header file for consistency check

double MaxArrayValue(const double* const Array, const int ArrayLen)
{
    // Routine for finding the maximum value in an Array of length 'ArrayLen'.
    const double* const max_p = std::max_element(Array,Array+ArrayLen);
    return *max_p;
}

double nF0(double const energy)
{
    /* Zero-temperature Fermi function. Overloaded for floats. */
    double ans = 1.;
    if (energy > 0.)
        ans = 0.;
    
    return ans;    
}
float  nF0(float const energy)
{
    /* Zero-temperature Fermi function. Overloaded for doubles. */
    float ans = 1.f;
    if (energy > 0.f)
        ans = 0.f;
    
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
    return 1.f/( 1.f + std::exp(beta*energy) );
}

double FermiEnerg_cpy(const int num_states, const int filled_states, 
                      double const *const energies, const bool print_bounds)
{
    /* Routine that finds the Fermi energy given the set of energy eigenvalues (of which 
    there are num_states) and the number of filled sates. The values in 'energies' are 
    copied in order to not be overwritten. Remember to properly round filled_states to an 
    integer when using. */
    
    // Allocate an array to hold a copy of energies
    double *const energies_copy = new double [num_states];
    // Copy the values
    for (int i=0; i<num_states; ++i)
        energies_copy[i] = energies[i];
    
    /* We call the sort function from the standard library. This function reorders the 
    values in the array energies_copy in ascending order. */
    std::sort(energies_copy, energies_copy+num_states);
    
    if (print_bounds)
    {
    std::cout << "Max energy, Min energy: " << energies_copy[0] << " " 
              << energies_copy[num_states-1] << "\t";
    }
        
    double FermiEnerg=-666; // To hold the Fermi energy
    
    if (filled_states<0)
        std::cout << "ERROR: in function FermiEnerg, filled_states cannot be negative."
                  << std::endl;
    else if (filled_states==0)
    {
        std::cout << "WARNING: in function FermiEnerg, filled_states is zero."
                  << std::endl;
        FermiEnerg = energies_copy[0];
    }
    else if (filled_states>num_states)
        std::cout << "ERROR: in function FermiEnerg, filled_states is too large."
                  << std::endl;
    else
        FermiEnerg = energies_copy[filled_states-1];
    
    
    // Deallocate memory for energies_copy
    delete [] energies_copy;
    
    return FermiEnerg;
}

double FermiEnerg(const int num_states, const int filled_states, 
                  double*const energies, const bool print_diagnostics)
{
    /* Routine that finds the Fermi energy given the set of energy eigenvalues (of which 
    there are num_states) and the number of filled sates. NOTE THAT VALUES IN 'energies' 
    ARE REORDERED. Remember to properly round filled_states to an integer when using. */
    
    /* We call the sort function from the standard library. This function reorders the 
    values in the array energies in ascending order. */
    std::sort(energies, energies+num_states);
        
    double FermiEnerg=-666; // To hold the Fermi energy
    
    if (filled_states<0)
        std::cout << "ERROR: in function FermiEnerg, filled_states cannot be negative."
                  << std::endl;
    else if (filled_states==0)
    {
        std::cout << "WARNING: in function FermiEnerg, filled_states is zero."
                  << std::endl;
        FermiEnerg = energies[0];
    }
    else if (filled_states>num_states)
        std::cout << "ERROR: in function FermiEnerg, filled_states is too large."
                  << std::endl;
    else
        FermiEnerg = energies[filled_states-1]; // Highest occupied state
    
    
    if (print_diagnostics)
    {
        std::cout << "Energies: ";
        
        std::cout << "min=" << energies[0] << " ";
        
        std::cout << "prev=";
        if ( (filled_states>1) && (filled_states<=num_states) )
            std::cout << energies[filled_states-2] << " ";
        else
            std::cout << "-- ";
        
        std::cout << "FE=" << FermiEnerg << " ";
        
        std::cout << "next=";
        if ( (filled_states>0) && (filled_states<num_states-1) )
            std::cout << energies[filled_states] << " ";
        else
            std::cout << "-- ";
        
        std::cout << "max=" << energies[num_states-1] << std::endl;
    }
    
    return FermiEnerg;
}
