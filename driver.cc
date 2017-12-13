// driver.cc
#include <iostream> // For terminal input and output
#include <complex> // For complex numbers
#include "math_routines.h" // Various math functions, such as the Fermi distribution
#include "alloc_dealloc.h" // Allocation of multidimensional arrays
#include "init_routines.h" // Array initialization routines

#include "ham1.h" // Choose correct Hamiltonian



int main()
{
    
    params_t params;
    Set_params(params);
    
    /* We allocate and initialize the momentum arrays. The parameters (number of grid 
    points and bounds) are obtained from the ham.cc module. */
    double*const kx_grid = new double [params.kx_pts];
    LinInitArray(params.kx_bounds[0], params.kx_bounds[1], params.kx_pts, kx_grid);
    
    double*const ky_grid = new double [params.ky_pts];
    LinInitArray(params.ky_bounds[0], params.ky_bounds[1], params.ky_pts, ky_grid);
    
    /* We allocate a 3D array to hold the energy eigenvalues (not including mu). The 
    first two indices correspond to kx and ky momentum and the third to the other quantum 
    numbers (as specified in the ham module). */
    double*const*const*const evals = Alloc3D_d(params.kx_pts, params.ky_pts, params.bands_num);
    ZeroInitArray(params.kx_pts*params.ky_pts*params.bands_num, evals[0][0]);
    
    
    
    
    
    Dealloc3D(evals, params.kx_pts);
    delete [] ky_grid;
    delete [] kx_grid;
    
    return 0;
}

double Evaluate_mu(const int kx_pts, const int ky_pts, const int bands_num, 
                   const double electrons_num, double*** evals)
{
    /* Returns mu for specified number of electrons. */
    double mu=0;
    return mu;
}