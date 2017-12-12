// driver.cc
#include <iostream> // For terminal input and output
#include <complex> // For complex numbers
#include "math_routines.h" // Various math functions, such as the Fermi distribution
#include "alloc_dealloc.h" // Allocation of multidimensional arrays
#include "init_routines.h" // Array initialization routines





int main()
{
    
    /* We allocate and initialize the momentum arrays. The parameters (number of grid 
    points and bounds) are obtained from the ham.cc module. */
    double* kx_grid = new double [kx_pts];
    InitArray(kx_bounds[0], kx_bounds[1], kx_pts, kx_grid);
    
    double* ky_grid = new double [ky_pts];
    InitArray(ky_bounds[0], kx_bounds[1], kx_pts, kx_grid);
    
    /* We allocate a 3D array to hold the energy eigenvalues (not including mu). The 
    first two indices correspond to kx and ky momentum and the third to the other quantum 
    numbers (as specified in the ham module). */
    double*** evals = Alloc3D_d(kx_pts, ky_pts, bands_num);
    ZeroInitArray(kx_pts*ky_pts*bands_num, evals);
    
    
    
    
    
    Dealloc3D(evals, kx_pts);
    delete [] ky_grid;
    delete [] kx_grid;
    
    return 0;
}