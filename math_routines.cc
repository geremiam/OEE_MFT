// math_routines.cc
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/

#include <iostream> // Input/output to command line
#include <cmath> // Numerous math functions
#include <complex> // Includes complex numbers
#include <algorithm> // Useful algorithms. Here, we use std::sort()
#include "math_routines.h" // Include header file for consistency check

double log_1p_exp(const double x)
{
    // Calculates the function log(1+exp(x)).
    // The function log1p takes care of x << -1.
    // If x >> 1, then answer is approximately x; use the conditional operator for this.
    return (x>log_max_double) ? x : log1p(exp(x));
}

void MinMaxArray(const double* const Array, const int ArrayLen, double& min, double& max)
{
    // Routine for finding the maximum value in an Array of length 'ArrayLen'.
    const double* const min_p = std::min_element(Array,Array+ArrayLen);
    const double* const max_p = std::max_element(Array,Array+ArrayLen);
    min = *min_p;
    max = *max_p;
}

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
double nF(const double T, const double energy)
{
    // Fermi function for doubles. Use long doubles to avoid overflow.
    const long double arg = (long double)(energy)/(long double)(T);
    const long double ans = 1.l/( 1.l + std::exp(arg) );
    return (double)(ans);
}
float nF(const float T, const float energy)
{
    // Fermi function for doubles. Use long doubles to avoid overflow.
    const long double arg = (long double)(energy)/(long double)(T);
    const long double ans = 1.l/( 1.l + std::exp(arg) );
    return (float)(ans);
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

std::complex<double> TraceMF(const int num_bands, 
                             const std::complex<double>*const*const evecs, 
                             const std::complex<double>*const mat, const double*const occs)
{
    /* Calculates traces of the form that comes up in mean-field computations. evecs is a 
    2D array whose columns are the evecs, mat is the matrix that defines the MF in the 
    form of a 1D array in row-major storage, and occs is a 1D array giving the occupations 
    of the bands (with the same choice of order as evecs) (given by the Fermi function).*/
    // Not great to implement mat. mul. by hand, be we will do it for simplicity
    // In row-major storage, "mat[b][c]" is mat[b*num_bands + c].
    std::complex<double> accumulator = {0.,0.};
    
    for (int a=0; a<num_bands; ++a)
      for (int b=0; b<num_bands; ++b)
        for (int c=0; c<num_bands; ++c)
          accumulator += conj(evecs[b][a]) * mat[b*num_bands+c] * evecs[c][a] * occs[a];
    
    return accumulator;
}

void Occupations(const int arrlen, const double mu, const double*const energies, 
                 double*const occs, const bool zerotemp, const double T)
{
    /* Calculates the occupations based on energies, mu, etc. */
    if (zerotemp)
      for (int i=0; i<arrlen; ++i)
        occs[i] = nF0(energies[i]-mu);
    else
      for (int i=0; i<arrlen; ++i)
        occs[i] = nF(T, energies[i]-mu);
}
