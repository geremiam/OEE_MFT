// ham1.cc
/* Defines the Hamiltonian (in row-major full storage layout) and relevant parameters. 
*/
#include <iostream>
#include <complex>
#include "init_routines.h"
#include "alloc_dealloc.h"
#include "kspace.h"
//#include "ham1.h" // Include module's header file for consistency check

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

/* Make sure the rectangular zone used here is equivalent to the first Brillouin zone. */
const double a = 1.;
const int kx_pts = 100;
const double kx_bounds [2] = {-pi/a, pi/a};
const int ky_pts = 100;
const double ky_bounds [2] = {-pi/a, pi/a};

const int bands_num = 2; // The number of bands, i.e. the order of the matrix for each k
const int ham_array_rows = bands_num; // Same as matrix order for full storage
const int ham_array_cols = bands_num; // Same as matrix order for full storage

/* We define the parameter space for the Hamiltonian, which will be scanned for the 
parameter study. */
const double rho = 1.; // Average electron density (FIXED)
const int U_pts = 10;
const double U_bounds [2] = {0., 1.};

double Dispersion(const double kx, const double ky)
{
    /* Non-interacting dispersion relation, excluding the chemical potential. */
    return 0.5 * (kx*kx + ky*ky);
}

void HamEval(const double kx, const double ky, const double M, 
             std::complex<double>** ham_array)
{
    /* Given the parameters kx, ky, and M, calculate the k-space Hamiltonian and assign 
    it to ham_array in full storage layout.*/
    
    const double epsilon = Dispersion(kx, ky);
    
    ham_array[0][0] = epsilon + U*rho/2. - U*M;
    ham_array[0][1] = 0.;
    ham_array[1][0] = 0.;
    ham_array[1][1] = epsilon + U*rho/2. + U*M;
}

double Evaluate_M(const int kx_pts, const int ky_pts, const int bands_num, 
                 const double mu, double*** evals)
{
    double accumulator = 0;
    for (int i=0; i<kx_pts; ++i)
        for (int j=0; j<ky_pts; ++j)
            accumulator += 0;
    
    // I think it's incorrect to divide by total number of states.
    return accumulator / (kx_pts*ky_pts*bands_num);
}


int main(int argc, char* argv[])
{
    if (int i=0; i<argc; ++i)
        std::cout << argv[i] << std::endl;
    
    // Declare (and construct) and instance of kspace_t
    kspace_t kspace(kx_pts, kx_bounds, ky_pts, ky_bounds, bands_num);
    
    // Declare an array to hold the Hamiltonian
    std::complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, ham_array[0]); // Initialize to zero
    
    // Define the parameter arrays
    double*const U_grid = new double [U_pts]; // Allocate memory
    const bool endpoint = true; // Include endpoint (not an important choice)
    LinInitArray(U_bounds[0], U_bounds[1], U_pts, U_grid, endpoint); // Initialize
    
    // Define the order parameter array, which spans the parameter space
    double*const M_grid = new double [U_pts]; // Allocate memory
    const double M_startval = 0.5; // Choose a starting value
    ValInitArray(U_pts, M_grid, M_startval); // Initialize to the starting value
    
    double M =      M_startval;
    double Mprime = M_startval;
    
    do
    {
        M = Mprime
        Mprime = 
    } while (std::abs(Mprime-M) < tol);
    
    
    // Deallocate memory
    delete [] M_grid;
    delete [] U_grid;
    Dealloc2D(ham_array);
    return 0;
}