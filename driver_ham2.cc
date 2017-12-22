// driver_ham2.cc
/* Mean-field theory of antiferromagnetism in the Haldane-Hubbard model. This driver 
defines model-specific parameters and functions, and performs the iteration until self-
consistency is achieved. The data is saved as a NetCDF dataset. */
#include <iostream>
#include <complex> // For complex numbers
#include <cmath> // For many math functions
#include <string>
#include "init_routines.h" // Array initialization
#include "alloc_dealloc.h" // Multidim array allocation
#include "math_routines.h" // Various custom math functions
#include "diag_routines.h" // Routines for finding evals and evecs
#include "kspace.h" // Defines a class for holding a band structure
#include "nc_IO.h" // Class for creating simple NetCDF datasets
using std::to_string;
//using std::cos; using std::sin; using std::conj; // Not sure if these are necessary
using std::polar;

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

/* Make sure the rectangular zone used here is equivalent to the first Brillouin zone. */
const double a = 1.; // We take a to be the NN distance
const int kx_pts = 173;
const double kx_bounds [2] = {-(2.*pi)/(3.*sqrt(3.)*a), (4.*pi)/(3.*sqrt(3.)*a)};
const int ky_pts = 200;
const double ky_bounds [2] = {-(2.*pi)/(3.*a), (2.*pi)/(3.*a)};

const int bands_num = 4; // The number of bands, i.e. the order of the matrix for each k
const int ham_array_rows = bands_num; // Same as matrix order for full storage
const int ham_array_cols = bands_num; // Same as matrix order for full storage

/* We define the parameter space for the Hamiltonian, which will be scanned for the 
parameter study. */
const double t1 = 1.; // t1 is set to 1, and so we measure all other energies in units of t1
const double t2 = 0.; // NNN hopping amplitude
const double phi = 0.; // Flux phase in the Haldane model
const double eps = 0.; // Potential difference between A and B sublattices
const double rho = 1.; // Average electron density (FIXED)
const int U_pts = 24; const double U_bounds [2] = {0., 7.}; // Hubbard interaction
const double M_startval = 0.1; // Choose a starting value

// Class that defines the parameter space for this Hamiltonian
class pspace_t
{
private:
    // Private copy constructor (prohibits copy creation)
    pspace_t(const pspace_t&);
    // Private assignment operator (prohibits assignment)
    const pspace_t& operator=(const pspace_t&);
    
    // Hubbard interaction
    const int U_pts;
    const double*const U_bounds; // Two-component array
public:
    double*const U_grid; // U coordinate variable
    
    double*const M_grid; // M variable
    
    // Constructor declaration
    pspace_t(const int U_pts_, const double*const U_bounds_)
        :U_pts(U_pts_), U_bounds(U_bounds_), U_grid(new double [U_pts_]),
        M_grid(new double [U_pts_])
    {
        /* The initialization list initialized the parameters and allocates memory for 
        the arrays. */
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(U_bounds[0], U_bounds[1], U_pts, U_grid, endpoint);
        
        // The other parameters (the order parameters) span the coordinate space
        ValInitArray(U_pts, M_grid); // Initialize to zero (for now)
        
        std::cout << "pspace_t instance created." << std::endl;
    }
    // Destructor declaration
    ~pspace_t() {
        delete [] U_grid;
        delete [] M_grid;
        std::cout << "pspace_t instance deleted." << std::endl;
    }
    
    void SaveData(const std::string GlobalAttr_, const std::string path_)
    {
        /* Method for saving the data of this class. This method uses the class from the 
        module nc_IO that creates a simple NetCDF class and allows writing of variables.*/
        /* We define parameters required to create the dataset. Don't forget to adjust 
        these depending on the parameter space defined above. */
        const int dims_num = 1;
        const std::string dim_names [dims_num] = {"U"};
        const int dim_lengths [dims_num] = {U_pts};
        const int vars_num = 1; // Variables other than coord variables
        const std::string var_names [vars_num] = {"M"};
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(GlobalAttr_, dims_num, dim_names, dim_lengths,
                      vars_num, var_names, path_);
        
        const double*const coord_vars [dims_num] = {U_grid};
        newDS.WriteCoordVars(coord_vars); // Write the coordinate variables
        
        const double*const vars [vars_num] = {M_grid};
        newDS.WriteVars(vars); // Write the variables
    }
};

void Dispersion(const double kx, const double ky, const double t2_, 
                std::complex<double>*const h)
{
    /* Assigns to h the values of the momentum-dependant 2*2 Haldane Hamiltonian. */
    h[0] = +eps - 2.*t2_*( 2.*cos(sqrt(3.)/2.*a*kx-phi)*cos(3./2.*a*ky) + cos(sqrt(3.)*a*kx+phi) );
    h[2] = -t1*( polar(1.0,a*ky) + polar(1.0,-a*ky/2.)*2.*cos(sqrt(3.)/2.*a*kx) );
    h[3] = conj(h[2]);
    h[4] = -eps - 2.*t2_*( 2.*cos(sqrt(3.)/2.*a*kx+phi)*cos(3./2.*a*ky) + cos(sqrt(3.)*a*kx-phi) );
}

void Evaluate_ham(const double kx, const double ky, const double t2_, 
                  const double rho_, const double U_, 
                  const double M_, std::complex<double>*const*const H)
{
    /* Given the parameters kx, ky, and M, calculate the 4*4 k-space Hamiltonian and 
    assign it to ham_array in full storage layout. */
    // We calculate the kinetic energy part of the Hamiltonian
    std::complex<double> h [4] = {0.,0.};
    Dispersion(kx, ky, t2_, h);
    
    H[0][0] = h[0]+U_*rho_/2.; H[0][1] = h[1];        H[0][2] = -U_*M_; H[0][3] = 0.;
    H[1][0] = h[2]; H[1][1] = h[3]+U_*rho_/2.;        H[1][2] = 0.;     H[1][3] = +U_*M_;
    
    H[2][0] = -U_*M_; H[2][1] = 0.;            H[2][2] = h[0]+U_*rho_/2.; H[2][3] = h[1];
    H[3][0] = 0.;     H[3][1] = +U_*M_;        H[3][2] = h[2]; H[3][3] = h[3]+U_*rho_/2.;
}

double Evaluate_M_term(const double mu, const double*const evals, 
                       const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to the OP from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double A [bands_num][bands_num] = {{0., 0., +1., 0.}, 
                                             {0., 0., 0., -1.},
                                             {+1., 0., 0., 0.},
                                             {0., -1., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int a=0; a<bands_num; ++a)
        for (int b=0; b<bands_num; ++b)
            for (int c=0; c<bands_num; ++c)
                accumulator += conj(evecs[b][a])*A[b][c]*evecs[c][a]*nF0(evals[a]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(4.*kx_pts*ky_pts));
    if (imag_part>1.e-15)
        std::cout << "WARNING: M has nonzero imaginary part: " << imag_part << std::endl;
    
    return std::real(accumulator/(4.*kx_pts*ky_pts));
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
    
    pspace_t pspace(U_pts, U_bounds); // Declare object of type pspace (parameter space)
    ValInitArray(U_pts, pspace.M_grid, M_startval); // Initialize to the starting value
    
    // Choose a tolerance for the equality of M and Mprime and print it.
    const double tol = 1.e-6;
    std::cout << "tol = " << std::scientific << tol << std::endl;
    
    for (int h=0; h<U_pts; ++h) // Loop over values of U
    {
    
        double M =      M_startval;
        double Mprime = M_startval;
        
        std::cout << "U = " << pspace.U_grid[h] << std::endl; // Print current value of U
        
        do // Iterate until self-consistency is achieved
        {
            M = Mprime;
            std::cout << "M = " << M << "\t";
            
            // Given the parameters, diagonalize the Hamiltonian at each grid point
            for (int i=0; i<kx_pts; ++i)
            for (int j=0; j<ky_pts; ++j)
            {
                Evaluate_ham(kspace.kx_grid[i], kspace.ky_grid[j], t2, rho, 
                             pspace.U_grid[h], M, ham_array);
                simple_zheev(bands_num, &(ham_array[0][0]), &(kspace.energies[i][j][0]));
            }
            
            // Use all the energies to compute the chemical potential
            const int num_states = kx_pts*ky_pts*bands_num;
            const int filled_states = 2*kx_pts*ky_pts*rho;//Be careful about lattice basis
            double mu = FermiEnerg(num_states, filled_states, &(kspace.energies[0][0][0]));
            std::cout << "mu = " << mu << "\t";
            
            // Use all the occupation numbers and the eigenvectors to find the order parameter
            // It is probably best to diagonalize a second time to avoid storing the evecs
            double accumulator = 0;
            
            for (int i=0; i<kx_pts; ++i)
            for (int j=0; j<ky_pts; ++j)
            {
                Evaluate_ham(kspace.kx_grid[i], kspace.ky_grid[j], t2, rho, 
                             pspace.U_grid[h], M, ham_array);
                simple_zheev(bands_num, &(ham_array[0][0]), &(kspace.energies[i][j][0]), 
                             true, &(evecs[0][0]));
                accumulator += Evaluate_M_term(mu, &(kspace.energies[i][j][0]), evecs);
            }
            Mprime = accumulator;
            // Print out final M value
            std::cout << "Mprime = " << Mprime << "\t";
            std::cout << "abs(Mprime-M) = " << std::abs(Mprime-M) << std::endl;
        } while (std::abs(Mprime-M) > tol);
        
        // We save the converged M value to the array pspace.M_grid.
        pspace.M_grid[h] = Mprime;
        std::cout << std::endl;
    }
    
    
    // Print out M array
    std::cout << std::endl;
    std::cout << "pspace.M_grid = " << std::endl;
    for (int h=0; h<U_pts; ++h)
        std::cout << pspace.M_grid[h] << " ";
    std::cout << std::endl;
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    const std::string GlobalAttr = "Haldane Hubbard model (ham2)"
        ": NN distance a = "+to_string(a)+"; kx_pts = "+to_string(kx_pts)+
        "; kx_bounds = "+to_string(kx_bounds[0])+", "+to_string(kx_bounds[1])+
        "; ky_pts = " + to_string(ky_pts)+
        "; ky_bounds = "+to_string(ky_bounds[0])+", "+to_string(ky_bounds[1])+
        "; bands_num = "+to_string(bands_num)+"; t1 = "+to_string(t1)+
        "; phi = "+to_string(phi)+"; eps = "+to_string(eps)+"; rho = "+to_string(rho)+
        "; M_startval = "+to_string(M_startval); // Define a string of metadata
        
    const std::string path="data/ham2/"; //Choose the path for saving (include final '/')
    
    pspace.SaveData(GlobalAttr, path); // Call saving method
    
    
    // Deallocate memory
    Dealloc2D(ham_array);
    Dealloc2D(evecs);
    return 0;
}