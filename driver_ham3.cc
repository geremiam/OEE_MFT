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
#include "nc_IO.h" // Class for creating simple NetCDF datasets
#include "IO.h" // For PrintMatrix()
#include "ham3.h" // Source code for ham3
using std::complex;
using std::string;


// Class that defines a parameter space for this Hamiltonian
class pspaceA_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceA_t(const pspaceA_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceA_t& operator=(const pspaceA_t&);
    
  public:
    // rho
    const double rho_bounds [2] = {0.3, 0.7}; const size_t rho_pts = 3;
    // V1
    const double V1_bounds [2] = {0., 1.};    const size_t V1_pts = 5;
    
    // Coordinate variables
    double*const rho_grid; // rho coordinate variable
    double*const V1_grid; // V1 coordinate variable
    
    /* Variables for MF parameters. We store them in 1D arrays, but we must keep track of 
    the order of the coordinate variable "indices". */
            double *const rho_a_grid;
    complex<double>*const u1_grid;
    complex<double>*const u1p_s_grid;
    complex<double>*const u1p_a_grid;
    complex<double>*const u2A_grid;
    complex<double>*const u2B_grid;
    complex<double>*const u3_s_grid;
    complex<double>*const u3_a_grid;
                int*const loops_grid; // holds the number of loops done at each point
            double *const energy_grid; // Holds the MF energy for later comparison
    
    // Constructor declaration
    pspaceA_t()
        :rho_grid(new double [rho_pts]), // Coord vars
         V1_grid (new double [V1_pts]), 
         rho_a_grid(new         double  [rho_pts*V1_pts]), // Vars
         u1_grid   (new complex<double> [rho_pts*V1_pts]),
         u1p_s_grid(new complex<double> [rho_pts*V1_pts]),
         u1p_a_grid(new complex<double> [rho_pts*V1_pts]),
         u2A_grid  (new complex<double> [rho_pts*V1_pts]),
         u2B_grid  (new complex<double> [rho_pts*V1_pts]),
         u3_s_grid (new complex<double> [rho_pts*V1_pts]),
         u3_a_grid (new complex<double> [rho_pts*V1_pts]),
         loops_grid(new         int     [rho_pts*V1_pts]),
         energy_grid(new        double  [rho_pts*V1_pts]) // Important -- Note order
    {
        /* The initialization list initializes the parameters and allocates memory for 
        the arrays. */
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(rho_bounds[0], rho_bounds[1], rho_pts, rho_grid, endpoint);
        LinInitArray(V1_bounds[0],  V1_bounds[1],  V1_pts,  V1_grid,  endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(rho_pts*V1_pts, rho_a_grid, -99.);//Init to unlikely value
        ValInitArray(rho_pts*V1_pts, u1_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, u1p_s_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, u1p_a_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, u2A_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, u2B_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, u3_s_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, u3_a_grid, -99.); // Likewise
        ValInitArray(rho_pts*V1_pts, loops_grid, 0); // Initialize to 0
        
        
        std::cout << "pspaceA_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceA_t()
    {
        delete [] rho_grid;
        delete [] V1_grid;
        delete [] rho_a_grid;
        delete [] u1_grid;
        delete [] u1p_s_grid;
        delete [] u1p_a_grid;
        delete [] u2A_grid;
        delete [] u2B_grid;
        delete [] u3_s_grid;
        delete [] u3_a_grid;
        delete [] loops_grid;
        delete [] energy_grid;
        
        std::cout << "pspaceA_t instance deleted.\n";
    }
    
    int Index(const int i, const int j)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. We pick the order rho, V1 so that V1 varies faster than rho.
        return i*V1_pts + j;
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        const size_t dims_num = 2;
        const string dim_names [dims_num] = {"rho", "V1"};
        const size_t dim_lengths [dims_num] = {rho_pts, V1_pts};
        const size_t vars_num = 8; // Variables other than coord variables
        const string var_names [vars_num] = {"rho_a", "u1", "u1p_s", "u1p_a", "u2A", "u2B", "u3_s", "u3_a"}; // Use same order below
        const bool var_complex [vars_num] = {false,   true, true,    true,    true,  true,  true,   true};
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
        
        newDS.DefCoordVar(0, "rho"); // Define "default" coordinate variables
        newDS.DefCoordVar(1, "V1");
        
        newDS.EndDef(); // Exit definition mode
        
        newDS.WriteCoordVar(0, rho_grid); // Write coordinate variables
        newDS.WriteCoordVar(1, V1_grid);
        
        const double*const vars [vars_num] = {rho_a_grid, 
                                              reinterpret_cast<double*const>(u1_grid), 
                                              reinterpret_cast<double*const>(u1p_s_grid), 
                                              reinterpret_cast<double*const>(u1p_a_grid), 
                                              reinterpret_cast<double*const>(u2A_grid), 
                                              reinterpret_cast<double*const>(u2B_grid), 
                                              reinterpret_cast<double*const>(u3_s_grid), 
                                              reinterpret_cast<double*const>(u3_a_grid)};
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid);
        newDS.WriteEnergy(energy_grid);
    }
};



// ######################################################################################
int pstudyA()
{
    /* This routine performs the mean-field iterative search at every point in the 
    parameter space defined above. */
    
    const bool with_output = true; // Show output for diagnostics
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    string GlobalAttr; // String to hold the metadata
    
    pspaceA_t pspaceA; // Declare object of type pspaceA (parameter space)
    
    // Declare and construct an instance of ham3_t
    ham3_t ham3(20,20,20);
    // Make any initial adjustment to the parameters
    ham3.set_zerotemp();
    
    GlobalAttr = ham3.GetAttributes(); // assign attributes to GlobalAttr
    std::cout << "\n\n\ttol = " << ham3.tol_ << "\n";//Print tolerance for equality of MFs.
    
    // Loop over values of the parameter space
    for (int g=0; g<pspaceA.rho_pts; ++g)
      for (int h=0; h<pspaceA.V1_pts; ++h)
      {
        if (with_output)
          std::cout << "\n\nrho = " << pspaceA.rho_grid[g] << ", "
                    << "V1 = "     << pspaceA.V1_grid[h] << "\n"; // Print current params
        
        // Adjust phase space parameters
        ham3.assign_rho(pspaceA.rho_grid[g]); // Assign value of rho
        ham3.V1_ = pspaceA.V1_grid[h]; // Assign value of V1
        
        int loops=0; // Will be assigned the number of loops performed
        const bool fail = ham3.FixedPoint(&loops, with_output);
        
        
        if (fail)
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "rho = " << pspaceA.rho_grid[g] << ", "
                      << "V1 = "     << pspaceA.V1_grid[h] << "\n"; //Print current params
            ++numfails;
        }
        
        // We save the MF parameters to the pspaceA arrays.
        pspaceA.rho_a_grid[pspaceA.Index(g,h)] = ham3.rho_a_;
        pspaceA.u1_grid[pspaceA.Index(g,h)]    = ham3.u1_;
        pspaceA.u1p_s_grid[pspaceA.Index(g,h)] = ham3.u1p_s_;
        pspaceA.u1p_a_grid[pspaceA.Index(g,h)] = ham3.u1p_a_;
        pspaceA.u2A_grid[pspaceA.Index(g,h)]   = ham3.u2A_;
        pspaceA.u2B_grid[pspaceA.Index(g,h)]   = ham3.u2B_;
        pspaceA.u3_s_grid[pspaceA.Index(g,h)]  = ham3.u3_s_;
        pspaceA.u3_a_grid[pspaceA.Index(g,h)]  = ham3.u3_a_;
        // Save the number of loops to pspaceA array.
        pspaceA.loops_grid[pspaceA.Index(g,h)] = loops;
        
        if (with_output) std::cout << std::endl;
      }
    
    
    // Print out MF parameters
    if (with_output)
    {
        std::cout << std::endl << "pspaceA.rho_a_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.rho_a_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u1_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u1_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u1p_s_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u1p_s_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u1p_a_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u1p_a_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u2A_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u2A_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u2B_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u2B_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u3_s_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u3_s_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.u3_a_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.u3_a_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.energy_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.energy_grid, std::cout);
        
        std::cout << std::endl << "pspaceA.loops_grid = " << std::endl;
        PrintMatrix(pspaceA.rho_pts, pspaceA.V1_pts, pspaceA.loops_grid, std::cout);
        std::cout << std::endl;
    }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    const string path="data/ham3/"; //Choose the path for saving (include final '/')
    pspaceA.SaveData(GlobalAttr, path); // Call saving method
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int main(int argc, char* argv[])
{
    int info = pstudyA();
    
    return info;
}