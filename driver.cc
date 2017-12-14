// driver.cc
#include <iostream> // For terminal input and output
#include <complex> // For complex numbers
#include "math_routines.h" // Various math functions, such as the Fermi distribution
#include "alloc_dealloc.h" // Allocation of multidimensional arrays
#include "init_routines.h" // Array initialization routines
#include "diag_routines.h" // Routines for finding evals and evecs

#include "ham1.h" // Choose correct Hamiltonian


int main()
{
    
    params_t par;
    
    /* We allocate and initialize the momentum arrays. The parameters (number of grid 
    points and bounds) are obtained from the ham.cc module. */
    double*const kx_grid = new double [par.kx_pts];
    LinInitArray(par.kx_bounds[0], par.kx_bounds[1], par.kx_pts, kx_grid);
    
    double*const ky_grid = new double [par.ky_pts];
    LinInitArray(par.ky_bounds[0], par.ky_bounds[1], par.ky_pts, ky_grid);
    
    /* We allocate a 3D array to hold the energy eigenvalues (not including mu). The 
    first two indices correspond to kx and ky momentum and the third to the other quantum 
    numbers (as specified in the ham module). */
    double*const*const*const evals = Alloc3D_d(par.kx_pts, par.ky_pts, par.bands_num);
    ZeroInitArray(par.kx_pts*par.ky_pts*par.bands_num, evals[0][0]);
    
    
    
    
    
    Dealloc3D(evals, par.kx_pts);
    delete [] ky_grid;
    delete [] kx_grid;
    
    return 0;
}