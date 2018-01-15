// driver_ham3.cc
/* Mean-field theory of antiferromagnetism in the Haldane bilayer, where layer II also 
has a Hubbard interaction. This driver defines model-specific parameters and functions, 
and performs the iteration until self-consistency is achieved. The data is saved as a 
NetCDF dataset. */
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

/* Size of the Hamiltonian, or equivalently number of bands (constant) */
const int bands_num = 8; // The number of bands, i.e. the order of the matrix for each k
const int ham_array_rows = bands_num; // Same as matrix order for full storage
const int ham_array_cols = bands_num; // Same as matrix order for full storage

/* Constants entering the definition of the Hamiltonian */
class pars_t {
  private:
  public:
    // Default values of the parameters
    double t1 = 1.; // NN hopping
    double t2 = 0.25; // NNN hopping
    double eps = 0.1; // Potential difference between A and B sublattices
    double phi = pi/2.; // Flux phase in the Haldane model
    pars_t(const double alpha_=1.) { // The energy parameters get scaled by alpha_.
        t1 *= alpha_;
        t2 *= alpha_;
        eps *= alpha_;
    }
};// Class for holding the parameters proper to each layer
const double tperp_0 = 0.3; // Base value of tperp (gets scaled)
const double L = 0.; // bias voltage between layers I and II
const double rho = 1.; // Average (global) electron density, between 0 and 2
/* Range and resolution of the parameter study */
const int alpha_pts = 6; const double alpha_bounds [2] = {1., 2.}; // scaling factor
const int U_pts = 6; const double U_bounds [2] = {0., 10.};//Hubbard interaction strength
/* Starting values for the order parameters */
const double M_startval = 0.1; // Choose a starting value
const double rhoI_startval = 1.; // Choose starting value

// Class that defines the parameter space for this Hamiltonian
class pspace_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspace_t(const pspace_t&);
    // Private assignment operator (prohibits assignment)
    const pspace_t& operator=(const pspace_t&);
    
    // scaling factor
    const int alpha_pts;
    const double*const alpha_bounds; // Two-component array
    // Hubbard interaction
    const int U_pts;
    const double*const U_bounds; // Two-component array
  public:
    double*const alpha_grid; // alpha coordinate variable
    double*const U_grid; // U coordinate variable
    
    double*const*const M_grid; // M variable
    double*const*const rhoI_grid; // rhoI variable
    
    // Constructor declaration
    pspace_t(const int alpha_pts_, const double*const alpha_bounds_, 
             const int U_pts_,  const double*const U_bounds_)
        :alpha_pts(alpha_pts_), alpha_bounds(alpha_bounds_), alpha_grid(new double [alpha_pts_]),
         U_pts(U_pts_),   U_bounds(U_bounds_),   U_grid(new double [U_pts_]),
         M_grid(Alloc2D_d(alpha_pts_, U_pts_)), 
         rhoI_grid(Alloc2D_d(alpha_pts_, U_pts_)) // Important -- Note order: alpha, U
    {
        /* The initialization list initializes the parameters and allocates memory for 
        the arrays. */
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(alpha_bounds[0], alpha_bounds[1], alpha_pts, alpha_grid, endpoint);
        LinInitArray(U_bounds[0],  U_bounds[1],  U_pts,  U_grid,  endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(alpha_pts*U_pts, &(M_grid[0][0]), -99.); // Init to unlikely value
        ValInitArray(alpha_pts*U_pts, &(rhoI_grid[0][0]), -99.); // Likewise
        
        
        std::cout << "Instance of pspace_t created." << std::endl;
    }
    // Destructor declaration
    ~pspace_t() {
        delete [] alpha_grid;
        delete [] U_grid;
        Dealloc2D(M_grid);
        Dealloc2D(rhoI_grid);
        std::cout << "Instance of pspace_t deleted." << std::endl;
    }
    
    void SaveData(const std::string GlobalAttr_, const std::string path_) {
        /* Method for saving the data of this class. This method uses the class from the 
        module nc_IO that creates a simple NetCDF class and allows writing of variables.*/
        /* We define parameters required to create the dataset. Don't forget to adjust 
        these depending on the parameter space defined above. 
        Important: Note that the order of the variables must be kept consistent. */
        const int dims_num = 2;
        const std::string dim_names [dims_num] = {"alpha", "U"};
        const int dim_lengths [dims_num] = {alpha_pts, U_pts};
        const int vars_num = 2; // Variables other than coord variables
        const std::string var_names [vars_num] = {"M", "rhoI"}; // Use same order below
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(GlobalAttr_, dims_num, dim_names, dim_lengths,
                      vars_num, var_names, path_);
        
        const double*const coord_vars [dims_num] = {alpha_grid, U_grid};
        newDS.WriteCoordVars(coord_vars); // Write the coordinate variables
        
        const double*const vars [vars_num] = {&(M_grid[0][0]), &(rhoI_grid[0][0])};
        newDS.WriteVars(vars); // Write the variables
    }
};

void Assign_h(const double kx, const double ky, const pars_t pars,
                std::complex<double>*const h)
{
    double t1 = pars.t1; double t2 = pars.t2; // Define local variables for convenience
    double eps = pars.eps; double phi = pars.phi;
    /* Assigns to h (a 1D array in row-major layout) the values of the momentum-dependant 
    2*2 Haldane Hamiltonian. */
    h[0] = +eps - 2.*t2*( 2.*cos(sqrt(3.)/2.*a*kx-phi)*cos(3./2.*a*ky) + cos(sqrt(3.)*a*kx+phi) );
    h[1] = -t1*( polar(1.,a*ky) + polar(1.,-a*ky/2.)*2.*cos(sqrt(3.)/2.*a*kx) );
    h[2] = conj(h[1]);
    h[3] = -eps - 2.*t2*( 2.*cos(sqrt(3.)/2.*a*kx+phi)*cos(3./2.*a*ky) + cos(sqrt(3.)*a*kx-phi) );
}

void Assign_ham(const double kx, const double ky, const pars_t parsI, const pars_t parsII, 
                const double tperp_, const double L_,  const double U_, 
                const double rhoI_, const double M_, std::complex<double>*const*const H)
{
    /* Given the parameters kx, ky, and M_, calculate the 8*8 k-space Hamiltonian and 
    assign it to H in full storage layout. The diag routine only uses the lower triangle, 
    so we only assign that part. */
    
    // We calculate the kinetic energy parts of the Hamiltonian
    std::complex<double> hI [4] = {0.,0.};
    Assign_h(kx, ky, parsI, hI);
    std::complex<double> hII [4] = {0.,0.};
    Assign_h(kx, ky, parsII, hII);
    
    H[0][0]=hI[0]+U_*rhoI_/2.+L_/2.; 
    H[1][0]=hI[2]; H[1][1]=hI[3]+U_*rhoI_/2.+L_/2.; 
    
    H[2][0]=-U_*M_; H[2][1]=0.;    H[2][2]=hI[0]+U_*rhoI_/2.+L_/2.; 
    H[3][0]=0.; H[3][1]=+U_*M_;    H[3][2]=hI[2]; H[3][3]=hI[3]+U_*rhoI_/2.+L_/2.; 
    
    H[4][0]=tperp_; H[4][1]=0.;    H[4][2]=0.; H[4][3]=0.;    H[4][4]=hII[0]-L_/2.; 
    H[5][0]=0.; H[5][1]=tperp_;    H[5][2]=0.; H[5][3]=0.;    H[5][4]=hII[2]; H[5][5]=hII[3]-L_/2.; 
    
    H[6][0]=0.; H[6][1]=0.;    H[6][2]=tperp_; H[6][3]=0.;    H[6][4]=0.; H[6][5]=0.;    H[6][6]=hII[0]-L_/2.; 
    H[7][0]=0.; H[7][1]=0.;    H[7][2]=0.; H[7][3]=tperp_;    H[7][4]=0.; H[7][5]=0.;    H[7][6]=hII[2]; H[7][7]=hII[3]-L_/2.;
    
}

double Compute_M_term(const double mu, const double*const evals, 
                      const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to M from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double C [bands_num][bands_num] = {{0., 0., +1., 0.,   0., 0., 0., 0.}, 
                                             {0., 0., 0., -1.,   0., 0., 0., 0.},
                                             {+1., 0., 0., 0.,   0., 0., 0., 0.},
                                             {0., -1., 0., 0.,   0., 0., 0., 0.},
                                             
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int b=0; b<bands_num; ++b)
        for (int c=0; c<bands_num; ++c)
            for (int d=0; d<bands_num; ++d)
                accumulator += conj(evecs[c][b])*C[c][d]*evecs[d][b]*nF0(evals[b]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(double)(4*kx_pts*ky_pts));
    if (imag_part>1.e-15)
        std::cerr << "WARNING: M has nonzero imaginary part: " << imag_part << std::endl;
    
    return std::real(accumulator/(double)(4*kx_pts*ky_pts));
}

double Compute_rhoI_term(const double mu, const double*const evals, 
                         const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to rhoI from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double B [bands_num][bands_num] = {{1., 0., 0., 0.,    0., 0., 0., 0.}, 
                                             {0., 1., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 1., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 1.,    0., 0., 0., 0.},
                                             
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int b=0; b<bands_num; ++b)
        for (int c=0; c<bands_num; ++c)
            for (int d=0; d<bands_num; ++d)
                accumulator += conj(evecs[c][b])*B[c][d]*evecs[d][b]*nF0(evals[b]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(double)(2*kx_pts*ky_pts));
    if (imag_part>1.e-15)
        std::cerr << "WARNING: rhoI has nonzero imaginary part: " << imag_part << std::endl;
    
    return std::real(accumulator/(double)(2*kx_pts*ky_pts));
}

// ######################################################################################
int main(int argc, char* argv[])
{
    const bool with_output = true; // Show output for diagnostics
    
    // Choose a tolerance for the equality of the mean fields and print it.
    const double tol = 1.e-6;
    std::cout << "\ntol = " << std::scientific << tol << std::endl << std::endl;
    
    // Declare object of type pspace (parameter space)
    pspace_t pspace(alpha_pts, alpha_bounds, U_pts, U_bounds);
    
    // Loop over values of the parameter space
    /* PARALLELIZATION:, note that different threads do not write to the same parts of 
    pspace. For some reason, kx_bounds and ky_bounds need to be declared as shared even 
    though they are const. In any event, they are only read and cannot be written to. */
    #pragma omp parallel default(none) shared(pspace,kx_bounds,ky_bounds,std::cout)
    {
    
    /* Declare (and construct) an instance of kspace_t (local to each thread). */
    kspace_t kspace(kx_pts, kx_bounds, ky_pts, ky_bounds, bands_num);
    /* Declare an array to hold the Hamiltonian (local to each thread). */
    std::complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); //Initialize to zero
    /* Declare an array to hold its evecs (local to each thread). */
    std::complex<double>*const*const evecs = Alloc2D_z(bands_num, bands_num);
    ValInitArray(bands_num*bands_num, &(evecs[0][0])); // Initialize to zero
    
    
    #pragma omp for schedule(dynamic,1)
    for (int g=0; g<alpha_pts; ++g)
      for (int h=0; h<U_pts;  ++h)
      {
        // Set the OPs to their startvals
        double M =      M_startval;
        double Mprime = M_startval;
        double rhoI =      rhoI_startval;
        double rhoIprime = rhoI_startval;
        
        // Define pars_t objects
        pars_t parsI; // parsI is not scaled
        pars_t parsII(pspace.alpha_grid[g]); // parsII is scaled by alpha
        
        double tperp = tperp_0 * pspace.alpha_grid[g]; // tperp is also scaled by alpha
        
        if (with_output)
            std::cout << "alpha = " << pspace.alpha_grid[g] << ", "
                      << "U = "  << pspace.U_grid[h] << std::endl; //Print current params
        
        do // Iterate until self-consistency is achieved
        {
            M = Mprime;
            rhoI = rhoIprime; // Update mean-field values
            if (with_output) std::cout << "M=" << M << ", rhoI=" << rhoI << "\t";
            
            // Given the parameters, diagonalize the Hamiltonian at each grid point
            for (int i=0; i<kx_pts; ++i)
              for (int j=0; j<ky_pts; ++j)
              {
                Assign_ham(kspace.kx_grid[i], kspace.ky_grid[j], parsI, parsII, tperp, L, pspace.U_grid[h], rhoI, M, ham_array);
                simple_zheev(bands_num, &(ham_array[0][0]), &(kspace.energies[i][j][0]));
              }
            
            // Use all the energies to compute the chemical potential
            // Be careful about lattice basis
            const int num_states = kx_pts*ky_pts*bands_num;
            const int filled_states = int( rho * double(4*kx_pts*ky_pts) );
            double mu = FermiEnerg(num_states, filled_states, &(kspace.energies[0][0][0]), true);
            if (with_output) std::cout << "mu=" << mu << "\t";
            
            // Use all the occupation numbers and the evecs to find the order parameter
            // Probably best to diagonalize a second time to avoid storing the evecs
            double M_accumulator = 0.;
            double rhoI_accumulator = 0.;
            
            for (int i=0; i<kx_pts; ++i)
              for (int j=0; j<ky_pts; ++j)
              {
                Assign_ham(kspace.kx_grid[i], kspace.ky_grid[j], parsI, parsII, tperp, L, pspace.U_grid[h], rhoI, M, ham_array);
                simple_zheev(bands_num, &(ham_array[0][0]), &(kspace.energies[i][j][0]), true, &(evecs[0][0]));
                M_accumulator +=       Compute_M_term(mu, &(kspace.energies[i][j][0]), evecs);
                rhoI_accumulator += Compute_rhoI_term(mu, &(kspace.energies[i][j][0]), evecs);
              }
            Mprime = M_accumulator;
            rhoIprime = rhoI_accumulator;
            // Print out final OP values
            if (with_output)
                std::cout << "Mprime=" << Mprime << ", rhoIprime=" << rhoIprime 
                          << "\tabs(Mprime-M)=" << std::abs(Mprime-M) 
                          << ", abs(rhoIprime-rhoI)=" << std::abs(rhoIprime-rhoI) 
                          << std::endl;
        } while ( (std::abs(Mprime-M)>tol) || (std::abs(rhoIprime-rhoI)>tol) );
        
        // We save the converged OP values to the pspace arrays.
        pspace.M_grid[g][h] = Mprime;
        pspace.rhoI_grid[g][h] = rhoIprime;
        
        if (with_output) std::cout << std::endl;
      }
    
    // Deallocate the memory
    Dealloc2D(ham_array);
    Dealloc2D(evecs);
    }
    
    
    // Print out M array
    if (with_output)
    {
        std::cout << std::endl << "pspace.M_grid = " << std::endl;
        for (int g=0; g<alpha_pts; ++g)
        {
            for (int h=0; h<U_pts;  ++h)
                std::cout << pspace.M_grid[g][h] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    // Print out rhoI array
    if (with_output)
    {
        std::cout << std::endl << "pspace.rhoI_grid = " << std::endl;
        for (int g=0; g<alpha_pts; ++g)
        {
            for (int h=0; h<U_pts;  ++h)
                std::cout << pspace.rhoI_grid[g][h] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    // Define a string of metadata
    const std::string GlobalAttr = "Haldane Hubbard bilayer (ham3)"
        ": NN distance a = "+to_string(a)+"; kx_pts = "+to_string(kx_pts)+
        "; kx_bounds = "+to_string(kx_bounds[0])+", "+to_string(kx_bounds[1])+
        "; ky_pts = " + to_string(ky_pts)+
        "; ky_bounds = "+to_string(ky_bounds[0])+", "+to_string(ky_bounds[1])+
        "; bands_num = "+to_string(bands_num)/*+"; t1 = "+to_string(t1)+
        "; phi = "+to_string(phi)+"; eps = "+to_string(eps)+"; rho = "+to_string(rho)+
        "; M_startval = "+to_string(M_startval)*/+"; tol = "+to_string(tol);
        
    const std::string path="data/ham3/"; //Choose the path for saving (include final '/')
    
    pspace.SaveData(GlobalAttr, path); // Call saving method
    
    return 0;
}
