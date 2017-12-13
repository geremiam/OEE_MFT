// ham1.cc
/* Defines the Hamiltonian (in row-major full storage layout) and relevant parameters. 
*/
#include <iostream>
#include <complex>
#include "ham1.h" // Include module's header file for consistency check
#include "math_routines.h"

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

/* Make sure the rectangular zone used here is equivalent to the first Brillouin zone. */
const double a = 1.;
const int kx_pts = 100;
const double kx_bounds [2] = {-pi/a, pi/a};
const int ky_pts = 100;
const double ky_bounds [2] = {-pi/a, pi/a};

/* TEMPORARY: Parameters for the Hamiltonian. We will have to figure out how to change 
these in order to do a parameter study. */
const double rho = 1.; // Average electron density
const double U = 0.1;


const int bands_num = 2; // The number of bands, i.e. the order of the matrix for each k
const int ham_array_rows = bands_num; // Same as matrix order for full storage
const int ham_array_cols = bands_num; // Same as matrix order for full storage

void Set_params(params_t& params)
{
    params.kx_pts = kx_pts;
    params.kx_bounds = kx_bounds;
    params.ky_pts = ky_pts;
    params.ky_bounds = ky_bounds;
    
    params.bands_num = bands_num;
    params.ham_array_rows = ham_array_rows;
    params.ham_array_cols = ham_array_cols;
}

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
    
    return accumulator / (kx_pts*ky_pts*bands_num);
}