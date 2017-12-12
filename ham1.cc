// ham1.cc
/* Defines the Hamiltonian (in row-major full storage layout) and relevant parameters. 
*/
#include <iostream>
#include "ham1.h" // Include module's header file for consistency check
#include "math_routines.h"

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

/* Make sure the rectangular zone used here is equivalent to the first Brillouin zone. */
const int kx_pts = 100;
const int kx_bounds [2] = {-pi/a, pi/a};
const int ky_pts = 100;
const int ky_bounds [2] = {-pi/a, pi/a};

/* TEMPORARY: Parameters for the Hamiltonian. We will have to figure out how to change 
these in order to do a parameter study. */
const double rho = 1.; // Average electron density
const double U = 
const double a = 1.;

const int bands_num = 2; // The number of bands, i.e. the order of the matrix for each k
const int ham_array_rows = bands_num; // Same as matrix order for full storage
const int ham_array_cols = bands_num; // Same as matrix order for full storage

double Dispersion(const double kx, const double ky)
{
    /* Non-interacting dispersion relation, excluding the chemical potential. */
    return 0.5 * (kx*kx + ky*ky);
}

void HamEval(const double kx, const double ky, const double M, 
             std::complex<double>** ham_array)
{
    /* Given the parameters kx, ky, and M, calculate the k-space Hamiltonian and assign 
    it to ham_array in packed storage layout.*/
    
    epsilon = Dispersion(kx, ky);
    
    ham_array[0][0] = epsilon + U*rho/2. - U*M;
    ham_array[0][1] = 0.;
    ham_array[1][0] = 0.;
    ham_array[1][1] = epsilon + U*rho/2. + U*M;
}

double EvaluateM(const int kx_pts, const int ky_pts, const int bands_num, 
                 const double mu, double*** evals)
{
    double accumulator = 0;
    for (int i=0; i<kx_pts; ++i)
        for (int j=0; j<ky_pts; ++j)
            accumulator += 
    
    return accumulator / (kx_pts*ky_pts*bands_num);
}