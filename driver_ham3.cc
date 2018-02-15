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


// Class that defines the parameter space for this Hamiltonian
class pspace_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspace_t(const pspace_t&);
    // Private assignment operator (prohibits assignment)
    const pspace_t& operator=(const pspace_t&);
    
  public:
    // scaling factor
    const int alpha_pts = 6; const double alpha_bounds [2] = {1., 3.};
    // Hubbard interaction strength
    const int U_pts = 6;     const double U_bounds [2] = {0., 10.};
    
    // Coordinate variables
    double*const alpha_grid; // alpha coordinate variable
    double*const U_grid; // U coordinate variable
    
    // MF parameters
    double*const*const rhoI_s_grid; // rhoI_s variable
    double*const*const rhoI_a_grid; // rhoI_a variable
    double*const*const mag_s_grid; // mag_s variable
    double*const*const mag_a_grid; // mag_a variable
    double*const*const loops_grid; // holds the number of loops done at each point
    
    // Constructor declaration
    pspace_t()
        :alpha_grid(new double [alpha_pts]), U_grid(new double [U_pts]),
         rhoI_s_grid(Alloc2D_d(alpha_pts, U_pts)),rhoI_a_grid(Alloc2D_d(alpha_pts, U_pts)),
         mag_s_grid(Alloc2D_d(alpha_pts, U_pts)), mag_a_grid(Alloc2D_d(alpha_pts, U_pts)),
         loops_grid(Alloc2D_d(alpha_pts, U_pts)) // Important -- Note order: alpha, U
    {
        /* The initialization list initializes the parameters and allocates memory for 
        the arrays. */
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(alpha_bounds[0], alpha_bounds[1], alpha_pts, alpha_grid, endpoint);
        LinInitArray(U_bounds[0],  U_bounds[1],  U_pts,  U_grid,  endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(alpha_pts*U_pts, &(rhoI_s_grid[0][0]), -99.);//Init to unlikely value
        ValInitArray(alpha_pts*U_pts, &(rhoI_a_grid[0][0]), -99.); // Likewise
        ValInitArray(alpha_pts*U_pts, &(mag_s_grid[0][0]), -99.); // Likewise
        ValInitArray(alpha_pts*U_pts, &(mag_a_grid[0][0]), -99.); // Likewise
        ValInitArray(alpha_pts*U_pts, &(loops_grid[0][0]), -99.); // Likewise
        
        
        std::cout << "pspace_t instance created.\n";
    }
    // Destructor declaration
    ~pspace_t() {
        delete [] alpha_grid;
        delete [] U_grid;
        Dealloc2D(rhoI_s_grid);
        Dealloc2D(rhoI_a_grid);
        Dealloc2D(mag_s_grid);
        Dealloc2D(mag_a_grid);
        Dealloc2D(loops_grid);
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
        const int vars_num = 5; // Variables other than coord variables
        const std::string var_names [vars_num] = {"rho_s", "rho_a", "mag_s", "mag_a", 
                                                  "loops"}; // Use same order below
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(GlobalAttr_, dims_num, dim_names, dim_lengths,
                      vars_num, var_names, path_);
        
        const double*const coord_vars [dims_num] = {alpha_grid, U_grid};
        newDS.WriteCoordVars(coord_vars); // Write the coordinate variables
        
        const double*const vars [vars_num] = {&(rhoI_s_grid[0][0]), &(rhoI_a_grid[0][0]),
                                              &(mag_s_grid[0][0]), &(mag_a_grid[0][0]),
                                              &(loops_grid[0][0])};
        newDS.WriteVars(vars); // Write the variables
    }
};

bool IterativeSearch(double& rhoI_s, double& rhoI_a, double& mag_s, double& mag_a, 
                     ham3_t& ham3, kspace_t& kspace, std::complex<double>*const*const evecs, 
                     int*const num_loops_p=NULL, const bool with_output=false)
{
    /* Performs the iterative self-consistent search using the parameters from ham3 and 
    the arrays kspace and evecs. The initial values of mag and rhoI are used as the 
    starting values for the search; the end values are also output to mag and rhoI. */
    
    /* For clarity, we define references for the MF parameter arguments, which are used 
    to output to. */
    double& rhoI_s_out = rhoI_s;
    double& rhoI_a_out = rhoI_a;
    double& mag_s_out  = mag_s;
    double& mag_a_out  = mag_a;
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        ham3.rhoI_s = rhoI_s_out;
        ham3.rhoI_a = rhoI_a_out; 
        ham3.mag_s  = mag_s_out;
        ham3.mag_a  = mag_a_out; // Update mean-field values
        if (with_output) std::cout << "rhoIs=" << ham3.rhoI_s
                                   << " rhoIa=" << ham3.rhoI_a
                                   << " ms="  << ham3.mag_s 
                                   << " ma="  << ham3.mag_a << "\t";
        
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
        double rhoI_s_accumulator = 0.;
        double rhoI_a_accumulator = 0.;
        double mag_s_accumulator = 0.;
        double mag_a_accumulator = 0.;
        
        for (int i=0; i<ham3.kx_pts; ++i)
          for (int j=0; j<ham3.ky_pts; ++j)
          {
            ham3.Assign_ham(kspace.kx_grid[i], kspace.ky_grid[j]);
            simple_zheev(ham3.bands_num, &(ham3.ham_array[0][0]), 
                                      &(kspace.energies[i][j][0]), true, &(evecs[0][0]));
            rhoI_s_accumulator += ham3.ComputeTerm_rhoI_s(mu,&(kspace.energies[i][j][0]),evecs);
            rhoI_a_accumulator += ham3.ComputeTerm_rhoI_a(mu,&(kspace.energies[i][j][0]),evecs);
            mag_s_accumulator  += ham3.ComputeTerm_mag_s(mu,&(kspace.energies[i][j][0]),evecs);
            mag_a_accumulator  += ham3.ComputeTerm_mag_a(mu,&(kspace.energies[i][j][0]),evecs);
          }
        rhoI_s_out = rhoI_s_accumulator;
        rhoI_a_out = rhoI_a_accumulator;
        mag_s_out  = mag_s_accumulator;
        mag_a_out  = mag_a_accumulator;
        
        // Print out final OP values
        if (with_output) std::cout << "drhoIs=" << rhoI_s_out - ham3.rhoI_s 
                                   << " drhoIa=" << rhoI_a_out - ham3.rhoI_a
                                   << " dms="  << mag_s_out - ham3.mag_s
                                   << " dma="  << mag_a_out - ham3.mag_a
                                   << std::endl;
        
        // Test for convergence
        converged =    (std::abs(rhoI_s_out - ham3.rhoI_s)<ham3.tol)
                    && (std::abs(rhoI_a_out - ham3.rhoI_a)<ham3.tol)
                    && (std::abs(mag_s_out - ham3.mag_s)<ham3.tol) 
                    && (std::abs(mag_a_out - ham3.mag_a)<ham3.tol);
        fail = (!converged) && (counter>ham3.loops_lim); // Must come after converged line
    } while (!converged && !fail);
    
    // Unless num_loops_p is the null pointer, assign the number of loops to its location
    if (num_loops_p!=NULL) *num_loops_p = counter;
    
    // We make sure that either converged or fail is true.
    if ((converged==true) && (fail==true) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (A)\n";
    if ((converged==false) && (fail==false) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (B)\n";
    
    return fail;
}


// ######################################################################################
int ParameterStudy()
{
    /* This routine performs the mean-field iterative search at every point in the 
    parameter space defined above. */
    
    const bool with_output = false; // Show output for diagnostics
    
    std::cout << std::scientific << std::showpos; // Format display output
    
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    std::string GlobalAttr; // String to hold the metadata
    
    pspace_t pspace; // Declare object of type pspace (parameter space)
    
    // Loop over values of the parameter space
    /* PARALLELIZATION:, note that different threads do not write to the same parts of 
    pspace. For some reason, kx_bounds and ky_bounds need to be declared as shared even 
    though they are const. In any event, they are only read and cannot be written to. */
    #pragma omp parallel default(none) shared(pspace,GlobalAttr,std::cout) reduction(+:numfails)
    {
    
    ham3_t ham3; // Declare and construct an instance of ham3_t (local to each thread)
    
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
    for (int g=0; g<pspace.alpha_pts; ++g)
      for (int h=0; h<pspace.U_pts;  ++h)
      {
        if (with_output)
          std::cout << "alpha = " << pspace.alpha_grid[g] << ", "
                    << "U = "     << pspace.U_grid[h] << "\n"; // Print current params
        
        // Adjust phase space parameters
        ham3.parsII.SetScaling(pspace.alpha_grid[g]); // parsII is scaled by alpha
        ham3.tperp = ham3.tperp_0 * pspace.alpha_grid[g];//tperp is also scaled by alpha
        ham3.U = pspace.U_grid[h];
        
        // Set the OPs to their startvals
        double rhoI_s = ham3.rhoI_s_startval;
        double rhoI_a = ham3.rhoI_a_startval;
        double mag_s = ham3.mag_s_startval;
        double mag_a = ham3.mag_a_startval;
        
        int loops=0; // Will receive the number of loops performed
        const bool fail = IterativeSearch(rhoI_s, rhoI_a, mag_s, mag_a, ham3, kspace, 
                                                             evecs, &loops, with_output);
        
        pspace.loops_grid[g][h] = loops; // Save the number of loops to pspace array.
        
        if (fail)
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "alpha = " << pspace.alpha_grid[g] << ", "
                      << "U = "     << pspace.U_grid[h] << "\n"; //Print current params
            ++numfails;
        }
        else
        {
            // We save the converged MF parameters to the pspace arrays.
            pspace.rhoI_s_grid[g][h] = rhoI_s;
            pspace.rhoI_a_grid[g][h] = rhoI_a;
            pspace.mag_s_grid[g][h] = mag_s;
            pspace.mag_a_grid[g][h] = mag_a;
        }
        
        if (with_output) std::cout << std::endl;
      }
    
    // Deallocate the memory
    Dealloc2D(evecs);
    }
    
    
    // Print out MF parameters
    if (with_output)
    {
        std::cout << std::endl << "pspace.rhoI_s_grid = " << std::endl;
        PrintMatrix(pspace.alpha_pts, pspace.U_pts, pspace.rhoI_s_grid, std::cout);
        std::cout << std::endl << "pspace.rhoI_a_grid = " << std::endl;
        PrintMatrix(pspace.alpha_pts, pspace.U_pts, pspace.rhoI_a_grid, std::cout);
        std::cout << std::endl << "pspace.mag_s_grid = " << std::endl;
        PrintMatrix(pspace.alpha_pts, pspace.U_pts, pspace.mag_s_grid, std::cout);
        std::cout << std::endl << "pspace.mag_a_grid = " << std::endl;
        PrintMatrix(pspace.alpha_pts, pspace.U_pts, pspace.mag_a_grid, std::cout);
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