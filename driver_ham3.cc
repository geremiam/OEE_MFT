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
#include "IO.h" // For PrintMatrix()
#include "ham3.h" // Source code for ham3


/* Range and resolution of the parameter study */
const int alpha_pts = 6; const double alpha_bounds [2] = {1., 3.}; // scaling factor
const int U_pts = 6; const double U_bounds [2] = {0., 10.};//Hubbard interaction strength


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
        
        
        std::cout << "pspace_t instance created.\n";
    }
    // Destructor declaration
    ~pspace_t() {
        delete [] alpha_grid;
        delete [] U_grid;
        Dealloc2D(M_grid);
        Dealloc2D(rhoI_grid);
        std::cout << "pspace_t instance deleted.\n";
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

bool IterativeSearch(double& Mprime, double& rhoIprime, ham3_t& ham3, kspace_t& kspace, 
                    std::complex<double>*const*const evecs, const bool with_output=false)
{
    /* Performs the iterative self-consistent search using the parameters from ham3 and 
    the arrays kspace and evecs. The initial values of Mprime and rhoIprime are used as 
    the starting values for the search; the end values are also output to Mprime and 
    rhoIprime. */
    
    // For clarity, we define references for the MF parameter arguments.
    double& mag_out  = Mprime;
    double& rhoI_out = rhoIprime;
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        ham3.mag = mag_out;
        ham3.rhoI = rhoI_out; // Update mean-field values
        if (with_output) std::cout << "mag=" << ham3.mag << ", rhoI=" << ham3.rhoI << "\t";
        
        // Given the parameters, diagonalize the Hamiltonian at each grid point
        for (int i=0; i<ham3.kx_pts; ++i)
          for (int j=0; j<ham3.ky_pts; ++j)
            {
            ham3.Assign_ham(kspace.kx_grid[i], kspace.ky_grid[j]);
            simple_zheev(ham3.bands_num, &(ham3.ham_array[0][0]), 
                                                        &(kspace.energies[i][j][0]));
            }
        
        // Use all energies to compute chemical potential (elements get reordered)
        double mu = FermiEnerg(ham3.num_states, ham3.filled_states, 
                                                        &(kspace.energies[0][0][0]));
        
        // Use all the occupation numbers and the evecs to find the order parameter
        // Probably best to diagonalize a second time to avoid storing the evecs
        double mag_accumulator = 0.;
        double rhoI_accumulator = 0.;
        
        for (int i=0; i<ham3.kx_pts; ++i)
          for (int j=0; j<ham3.ky_pts; ++j)
            {
            ham3.Assign_ham(kspace.kx_grid[i], kspace.ky_grid[j]);
            simple_zheev(ham3.bands_num, &(ham3.ham_array[0][0]), 
                                  &(kspace.energies[i][j][0]), true, &(evecs[0][0]));
            mag_accumulator  += ham3.ComputeTerm_mag(mu,  &(kspace.energies[i][j][0]), evecs);
            rhoI_accumulator += ham3.ComputeTerm_rhoI(mu, &(kspace.energies[i][j][0]), evecs);
          }
        mag_out  = mag_accumulator;
        rhoI_out = rhoI_accumulator;
        // Print out final OP values
        if (with_output)
            std::cout << "mag_out=" << mag_out << ", rhoI_out=" << rhoI_out 
                      << "\tmag_out-mag=" << mag_out-ham3.mag 
                      << ", rhoI_out-rhoI=" << rhoI_out-ham3.rhoI 
                      << std::endl;
        
        // Test for convergence
        converged = (std::abs(mag_out-ham3.mag)<ham3.tol) && (std::abs(rhoI_out-ham3.rhoI)<ham3.tol);
        fail = (!converged) && (counter>ham3.loops_lim); // Must come after converged line
    } while (!converged && !fail);
    
    if ((converged==true) && (fail==true) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (A)\n";
    if ((converged==false) && (fail==false) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (B)\n";
    
    return fail;
}


// ######################################################################################
int ParameterStudy()
{
    const bool with_output = false; // Show output for diagnostics
    
    // Format display output
    std::cout << std::scientific << std::showpos;
    
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    std::string GlobalAttr; // String to hold the metadata
    
    // Declare object of type pspace (parameter space)
    pspace_t pspace(alpha_pts, alpha_bounds, U_pts, U_bounds);
    
    // Loop over values of the parameter space
    /* PARALLELIZATION:, note that different threads do not write to the same parts of 
    pspace. For some reason, kx_bounds and ky_bounds need to be declared as shared even 
    though they are const. In any event, they are only read and cannot be written to. */
    #pragma omp parallel default(none) shared(pspace,GlobalAttr,std::cout) reduction(+:numfails)
    {
    
    // Declare and construct an instance of ham3_t (local to each thread)
    ham3_t ham3;
    
    #pragma omp single
    {
      GlobalAttr = ham3.GetAttributes();//Single thread assigns attributes to GlobalAttr
      std::cout << "\ntol = " << ham3.tol << "\n\n";//Print tolerance for equality of MFs.
    }
    
    /* Declare (and construct) an instance of kspace_t (local to each thread). */
    kspace_t kspace(ham3.kx_pts, ham3.kx_bounds, ham3.ky_pts, ham3.ky_bounds, 
                                                                         ham3.bands_num);
    /* Declare an array to hold the evecs (local to each thread). */
    std::complex<double>*const*const evecs = Alloc2D_z(ham3.bands_num, ham3.bands_num);
    ValInitArray(ham3.bands_num*ham3.bands_num, &(evecs[0][0])); // Initialize to zero
    
    
    #pragma omp for schedule(dynamic,1)
    for (int g=0; g<alpha_pts; ++g)
      for (int h=0; h<U_pts;  ++h)
      {
        if (with_output)
            std::cout << "alpha = " << pspace.alpha_grid[g] << ", "
                      << "U = "  << pspace.U_grid[h] << std::endl; //Print current params
        
        // Adjust phase space parameters
        ham3.parsII.SetScaling(pspace.alpha_grid[g]); // parsII is scaled by alpha
        ham3.tperp = ham3.tperp_0 * pspace.alpha_grid[g];//tperp is also scaled by alpha
        ham3.U = pspace.U_grid[h];
        
        // Set the OPs to their startvals
        double mag = ham3.mag_startval;
        double rhoI = ham3.rhoI_startval;
        
        const bool fail = IterativeSearch(mag, rhoI, ham3, kspace, evecs, with_output);
        
        if (fail)
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "alpha = " << pspace.alpha_grid[g] << ", "
                      << "U = "     << pspace.U_grid[h] << "\n"; //Print current params
            ++numfails;
        }
        else
        {
            // We save the converged OP values to the pspace arrays.
            pspace.M_grid[g][h] = mag;
            pspace.rhoI_grid[g][h] = rhoI;
        }
        
        if (with_output) std::cout << std::endl;
      }
    
    // Deallocate the memory
    Dealloc2D(evecs);
    }
    
    
    // Print out M and rhoI arrays
    if (with_output)
    {
        std::cout << std::endl << "pspace.M_grid = " << std::endl;
        PrintMatrix(alpha_pts, U_pts, pspace.M_grid, std::cout);
        std::cout << std::endl << "pspace.rhoI_grid = " << std::endl;
        PrintMatrix(alpha_pts, U_pts, pspace.rhoI_grid, std::cout);
        std::cout << std::endl;
    }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    const std::string path="data/ham3/"; //Choose the path for saving (include final '/')
    pspace.SaveData(GlobalAttr, path); // Call saving method
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int main(int argc, char* argv[])
{
    int info = ParameterStudy();
    
    return info;
}