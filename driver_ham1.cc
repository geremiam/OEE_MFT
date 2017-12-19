// ham1.cc
/* Defines the Hamiltonian (in row-major full storage layout) and relevant parameters. 
*/
#include <iostream>
#include <complex> // For complex numbers
#include <cmath> // For many math functions
#include "init_routines.h" // Array initialization
#include "alloc_dealloc.h" // Multidim array allocation
#include "math_routines.h" // Various custom math functions
#include "diag_routines.h" // Routines for finding evals and evecs
#include "kspace.h" // Defines a class for holding a band structure
#include "nc_IO.h" // Class for creating simple NetCDF datasets

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

/* Make sure the rectangular zone used here is equivalent to the first Brillouin zone. */
const double a = 1.;
const int kx_pts = 200;
const double kx_bounds [2] = {-pi/a, pi/a};
const int ky_pts = 200;
const double ky_bounds [2] = {-pi/a, pi/a};

const int bands_num = 2; // The number of bands, i.e. the order of the matrix for each k
const int ham_array_rows = bands_num; // Same as matrix order for full storage
const int ham_array_cols = bands_num; // Same as matrix order for full storage

/* We define the parameter space for the Hamiltonian, which will be scanned for the 
parameter study. */
const double t = 1.;
const double rho = 1.; // Average electron density (FIXED)
const int U_pts = 11;
const double U_bounds [2] = {0., 4.};

double Dispersion(const double kx, const double ky)
{
    /* Non-interacting dispersion relation, excluding the chemical potential. */
    return -2.*t*(std::cos(a*kx)+std::cos(a*ky));
}

void Evaluate_ham(const double kx, const double ky, const double U, const double rho, 
                  const double M, std::complex<double>*const*const ham_array)
{
    /* Given the parameters kx, ky, and M, calculate the k-space Hamiltonian and assign 
    it to ham_array in full storage layout. */
    const double epsilon = Dispersion(kx, ky);
    
    ham_array[0][0] = epsilon + U*rho/2.;
    ham_array[0][1] = - U*M;
    ham_array[1][0] = - U*M;
    ham_array[1][1] = epsilon + U*rho/2.;
}

// Evaluates the contribution to the OP from a single k (see notes)
double Evaluate_M_term(const double mu, const double*const evals, 
                       const std::complex<double>*const*const evecs)
{
    // Not good to implement matrix mult. by hand... we will do it for simplicity
    const double sigma1 [bands_num][bands_num] = {{0., 1.}, 
                                                  {1., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int a=0; a<bands_num; ++a)
        for (int b=0; b<bands_num; ++b)
            for (int c=0; c<bands_num; ++c)
                accumulator += std::conj(evecs[b][a]) * sigma1[b][c] * evecs[c][a] * nF0(evals[a] - mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(2.*kx_pts*ky_pts));
    if (imag_part>1.e-15)
        std::cout << "WARNING: M has nonzero imaginary part: " << imag_part << std::endl;
    
    return std::real(accumulator/(2.*kx_pts*ky_pts));
}



// ######################################################################################
int main(int argc, char* argv[])
{
    // Declare (and construct) and instance of kspace_t
    kspace_t kspace(kx_pts, kx_bounds, ky_pts, ky_bounds, bands_num);
    
    // Declare an array to hold the Hamiltonian
    std::complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); // Initialize to zero
    
    // Declare an array to hold its evecs
    std::complex<double>*const*const evecs = Alloc2D_z(bands_num, bands_num);
    ValInitArray(bands_num*bands_num, &(evecs[0][0])); // Initialize to zero
    
    // Define the parameter arrays
    double*const U_grid = new double [U_pts]; // Allocate memory
    const bool endpoint = true; // Include endpoint (not an important choice)
    LinInitArray(U_bounds[0], U_bounds[1], U_pts, U_grid, endpoint); // Initialize
    
    // TEMPORARY: for testing, take a single U value
    //const double U = 0.6;
    
    // Define the order parameter array, which spans the parameter space
    double*const M_grid = new double [U_pts]; // Allocate memory
    const double M_startval = 0.1; // Choose a starting value
    ValInitArray(U_pts, M_grid, M_startval); // Initialize to the starting value
    
    
    
    const double tol = 1.e-15;
    std::cout << "tol = " << std::scientific << tol << std::endl;
    
    for (int h=0; h<U_pts; ++h) // Loop over values of U
    {
    
        double M =      M_startval;
        double Mprime = M_startval;
        
        std::cout << "U = " << U_grid[h] << std::endl; // Print current value of U
        
        do // Iterate until self-consistency is achieved
        {
            M = Mprime;
            std::cout << "M = " << M << "\t";
            
            // Given the parameters, diagonalize the Hamiltonian at each grid point
            for (int i=0; i<kx_pts; ++i)
                for (int j=0; j<ky_pts; ++j)
                    {
                        Evaluate_ham(kspace.kx_grid[i], kspace.ky_grid[j], U_grid[h], rho, M, ham_array);
                        simple_zheev(bands_num, &(ham_array[0][0]), &(kspace.energies[i][j][0]));
                    }
            
            // Use all the energies to compute the chemical potential
            const int num_states = kx_pts*ky_pts*bands_num;
            const int filled_states = kx_pts*ky_pts*rho; // Be careful about lattice basis
            double mu = FermiEnerg(num_states, filled_states, &(kspace.energies[0][0][0]));
            std::cout << "mu = " << mu << "\t";
            
            // Use all the occupation numbers and the eigenvectors to find the order parameter
            // It is probably best to diagonalize a second time to avoid storing the evecs
            double accumulator = 0;
            
            for (int i=0; i<kx_pts; ++i)
                for (int j=0; j<ky_pts; ++j)
                    {
                        Evaluate_ham(kspace.kx_grid[i], kspace.ky_grid[j], U_grid[h], rho, M, ham_array);
                        simple_zheev(bands_num, &(ham_array[0][0]), &(kspace.energies[i][j][0]), 
                                     true, &(evecs[0][0]));
                        accumulator += Evaluate_M_term(mu, 
                                                         &(kspace.energies[i][j][0]), evecs);
                    }
            Mprime = accumulator;
            // Print out final M value
            std::cout << "Mprime = " << Mprime << "\t";
            std::cout << "abs(Mprime-M) = " << std::abs(Mprime-M) << std::endl;
        } while (std::abs(Mprime-M) > tol);
        
        M_grid[h] = Mprime;
        std::cout << std::endl;
    }
    
    
    
    // Print out M array
    std::cout << std::endl;
    std::cout << "M_grid = " << std::endl;
    for (int h=0; h<U_pts; ++h)
        std::cout << M_grid[h] << " ";
    std::cout << std::endl;
    
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    const std::string GlobalAttr = "The global attributes";
    const std::string path="data/ham1/";
    
    const int dims_num_ = 1;
    const std::string dim_names [dims_num_] = {"U"};
    const int dim_lengths [dims_num_] = {U_pts};
    
    const int vars_num_ = 1;
    const std::string var_names [vars_num_] = {"M"};
    
    newDS_t newDS(GlobalAttr, dims_num_, dim_names, dim_lengths,
                  vars_num_, var_names, path);
    
    const double*const coord_vars [dims_num_] = {U_grid};
    newDS.WriteCoordVars(coord_vars);
    
    const double*const vars [vars_num_] = {M_grid};
    newDS.WriteVars(vars);
    
    
    
    // Deallocate memory
    delete [] M_grid;
    delete [] U_grid;
    Dealloc2D(ham_array);
    Dealloc2D(evecs);
    return 0;
}