// driver_ham4.cc
/* Mean-field theory of antiferromagnetism in the Haldane bilayer, where layer II also 
has a Hubbard interaction. This driver defines model-specific parameters and functions, 
and performs the iteration until self-consistency is achieved. The data is saved as a 
NetCDF dataset. */
#include <iostream>
#include <sstream> // For stringstreams
#include <complex> // For complex numbers
#include <cmath> // For many math functions
#include <string>
#include "alloc.h"
#include "init_routines.h" // Array initialization
#include "nc_IO.h" // Class for creating simple NetCDF datasets
#include "IO.h" // For PrintMatrix()
#include "ham4.h" // Source code for ham4
using std::complex;
using std::string;

class pspaceA_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceA_t(const pspaceA_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceA_t& operator=(const pspaceA_t&);
    
  public:
    // rho
    const size_t rho_pts = 5;
    const double rho_bounds [2] = {0.4, 0.6};
    // g: determines all the interaction strengths
    const size_t g_pts = 11;
    const double V1_bounds  [2] = {0.,5.0};
    const double V1p_bounds [2] = {0.,2.5};
    const double V2_bounds  [2] = {0.,2.5};
    const double V3_bounds  [2] = {0.,2.5};
    
    const int parspace_pts = rho_pts*g_pts;
    
    
    // Coordinate variables
    double*const rho_grid;// rho coordinate variable
    double*const V1_grid; //  V1 coordinate variable (depends on g)
    double*const V1p_grid;// V1p coordinate variable (depends on g)
    double*const V2_grid; //  V2 coordinate variable (depends on g)
    double*const V3_grid; //  V3 coordinate variable (depends on g)
    
    // Variables for MF parameters. We store them in 2D arrays. The first index is the 
    // harmonic, while the next keeps track of the point in parameter space (using the 
    // method index()). 
    const size_t num_harmonics_; // We need the number of harmonics being considered.
    double*const*const rho_s_grid; // 2D array
    double*const*const rho_a_grid; // 2D array
    int*const loops_grid; // holds the number of loops done at each point
    double*const energy_grid; // Holds the MF energy for later comparison
    double*const mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspaceA_t(const int num_harmonics)
        :num_harmonics_(num_harmonics),
         rho_grid(new double [rho_pts]), // Coord vars
         V1_grid (new double [g_pts]),
         V1p_grid(new double [g_pts]),
         V2_grid (new double [g_pts]),
         V3_grid (new double [g_pts]),
         rho_s_grid (Alloc2D_d(num_harmonics, parspace_pts)), // Vars
         rho_a_grid (Alloc2D_d(num_harmonics, parspace_pts)),
         loops_grid (new int    [parspace_pts]),
         energy_grid(new double [parspace_pts]),
         mu_grid    (new double [parspace_pts]) // Important -- Note order
    {
        // The initialization list initializes the parameters and allocates memory for 
        // the arrays.
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(rho_bounds[0], rho_bounds[1], rho_pts, rho_grid, endpoint);
        
        LinInitArray( V1_bounds[0],  V1_bounds[1], g_pts,  V1_grid, endpoint);
        LinInitArray(V1p_bounds[0], V1p_bounds[1], g_pts, V1p_grid, endpoint);
        LinInitArray( V2_bounds[0],  V2_bounds[1], g_pts,  V2_grid, endpoint);
        LinInitArray( V3_bounds[0],  V3_bounds[1], g_pts,  V3_grid, endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(num_harmonics * parspace_pts, &(rho_s_grid[0][0]), -99.);//Init to unlikely value
        ValInitArray(num_harmonics * parspace_pts, &(rho_a_grid[0][0]), -99.);//Init to unlikely value
        
        ValInitArray(parspace_pts, loops_grid, 0); // Initialize to 0
        ValInitArray(parspace_pts, energy_grid, -99.);
        ValInitArray(parspace_pts, mu_grid,      -9.);
        
        std::cout << "pspaceA_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceA_t()
    {
        delete [] rho_grid; // Coordinate variables
        delete [] V1_grid; delete [] V1p_grid;
        delete [] V2_grid; delete [] V3_grid;
        
        Dealloc2D(rho_s_grid); // MF variables
        Dealloc2D(rho_a_grid); // MF variables
        delete [] loops_grid; // Other variables
        delete [] energy_grid;
        delete [] mu_grid;
        
        std::cout << "pspaceA_t instance deleted.\n";
    }
    
    int Index(const int f, const int g)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. We pick the order rho, V1 so that V1 varies faster than rho.
        return g + 
               g_pts*f;
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        const size_t dims_num = 2;
        const string dim_names [dims_num] = {"rho", "g"};
        const size_t dim_lengths [dims_num] = {rho_pts, g_pts};
        const size_t vars_num = 2*num_harmonics_; // Variables other than coord variables
        
        string var_names [2*num_harmonics_]; // List for the variable names
        bool var_complex [2*num_harmonics_]; // List for indicating whether vars are complex
        for (int Q=0; Q<num_harmonics_; ++Q) { // Loop over harmonics to fill up lists
          std::ostringstream strs_s, strs_a;
          strs_s << "rho_s[" << Q << "]";
          strs_a << "rho_a[" << Q << "]";
          var_names[2*Q]   = strs_s.str(); // Use same order below
          var_names[2*Q+1] = strs_a.str();
          
          var_complex[2*Q]   = false;
          var_complex[2*Q+1] = false;
        }
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
        
        newDS.DefCoordVar(0, "rho"); // Define "default" coordinate variables
        // Define "custom" coordinate variables
        const int varid_V1 = newDS.DefCoordVar_custom(1, "V1"); // Vary with g
        const int varid_V1p= newDS.DefCoordVar_custom(1, "V1p");
        const int varid_V2 = newDS.DefCoordVar_custom(1, "V2");
        const int varid_V3 = newDS.DefCoordVar_custom(1, "V3");
        
        newDS.EndDef(); // Exit definition mode
        
        newDS.WriteCoordVar(0, rho_grid); // Write "default" coordinate variables
        newDS.WriteCoordVar_custom(varid_V1,  V1_grid); // Write "custom" coordinate variables
        newDS.WriteCoordVar_custom(varid_V1p, V1p_grid);
        newDS.WriteCoordVar_custom(varid_V2,  V2_grid);
        newDS.WriteCoordVar_custom(varid_V3,  V3_grid);
        
        double* vars [2*num_harmonics_]; // List for holding the pointers to the vars
        for (int Q=0; Q<num_harmonics_; ++Q) { // Fill up that list
          vars[2*Q]   = rho_s_grid[Q];
          vars[2*Q+1] = rho_a_grid[Q];
        }
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid); // Write loops variable
        newDS.WriteEnergy(energy_grid); // Write energy variable
        newDS.Writemu(mu_grid); // Write mu variable
    }
};

class pspaceAA_t { // Interaction strength and temperature varied at constant filling
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceAA_t(const pspaceAA_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceAA_t& operator=(const pspaceAA_t&);
    
  public:
    // T (temperature)
    const size_t T_pts = 5;
    // g: determines all the interaction strengths
    const size_t g_pts = 11;
    const double V1_bounds  [2] = {0.,5.0};
    const double V1p_bounds [2] = {0.,2.5};
    const double V2_bounds  [2] = {0.,2.5};
    const double V3_bounds  [2] = {0.,2.5};
    
    const int parspace_pts = T_pts*g_pts;
    
    
    // Coordinate variables
    double*const T_grid;// T (temperature) coordinate variable
    double*const V1_grid; //  V1 coordinate variable (depends on g)
    double*const V1p_grid;// V1p coordinate variable (depends on g)
    double*const V2_grid; //  V2 coordinate variable (depends on g)
    double*const V3_grid; //  V3 coordinate variable (depends on g)
    
    // Variables for MF parameters. We store them in 2D arrays. The first index is the 
    // harmonic, while the next keeps track of the point in parameter space (using the 
    // method index()). 
    const size_t num_harmonics_; // We need the number of harmonics being considered.
    double*const*const rho_s_grid; // 2D array
    double*const*const rho_a_grid; // 2D array
    
    int*const loops_grid; // holds the number of loops done at each point
    double*const energy_grid; // Holds the MF energy for later comparison
    double*const mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspaceAA_t(const int num_harmonics)
        :num_harmonics_(num_harmonics),
         T_grid  (new double [T_pts]), // Coord vars
         V1_grid (new double [g_pts]),
         V1p_grid(new double [g_pts]),
         V2_grid (new double [g_pts]),
         V3_grid (new double [g_pts]),
         rho_s_grid (Alloc2D_d(num_harmonics, parspace_pts)), // Vars
         rho_a_grid (Alloc2D_d(num_harmonics, parspace_pts)),
         loops_grid (new int    [parspace_pts]),
         energy_grid(new double [parspace_pts]),
         mu_grid    (new double [parspace_pts]) // Important -- Note order
    {
        // The initialization list initializes the parameters and allocates memory for 
        // the arrays.
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        //LinInitArray(rho_bounds[0], rho_bounds[1], T_pts, rho_grid, endpoint);
        T_grid[0] = 1.e-3;
        T_grid[1] = 1.e-2;
        T_grid[2] = 0.1;
        T_grid[3] = 1.0;
        T_grid[4] = 10.;
        
        LinInitArray( V1_bounds[0],  V1_bounds[1], g_pts,  V1_grid, endpoint);
        LinInitArray(V1p_bounds[0], V1p_bounds[1], g_pts, V1p_grid, endpoint);
        LinInitArray( V2_bounds[0],  V2_bounds[1], g_pts,  V2_grid, endpoint);
        LinInitArray( V3_bounds[0],  V3_bounds[1], g_pts,  V3_grid, endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(num_harmonics * parspace_pts, &(rho_s_grid[0][0]), -99.);//Init to unlikely value
        ValInitArray(num_harmonics * parspace_pts, &(rho_a_grid[0][0]), -99.);//Init to unlikely value
        
        ValInitArray(parspace_pts, loops_grid, 0); // Initialize to 0
        ValInitArray(parspace_pts, energy_grid, -99.);
        ValInitArray(parspace_pts, mu_grid,      -9.);
        
        std::cout << "pspaceAA_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceAA_t()
    {
        delete [] T_grid; // Coordinate variables
        delete [] V1_grid; delete [] V1p_grid;
        delete [] V2_grid; delete [] V3_grid;
        
        Dealloc2D(rho_s_grid); // MF variables
        Dealloc2D(rho_a_grid); // MF variables
        delete [] loops_grid; // Other variables
        delete [] energy_grid;
        delete [] mu_grid;
        
        std::cout << "pspaceAA_t instance deleted.\n";
    }
    
    int Index(const int f, const int g)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. We pick the order T, V1 so that V1 varies faster than rho.
        return g + 
               g_pts*f;
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        const size_t dims_num = 2;
        const string dim_names [dims_num] = {"T", "g"};
        const size_t dim_lengths [dims_num] = {T_pts, g_pts};
        const size_t vars_num = 2*num_harmonics_; // Variables other than coord variables
        
        string var_names [2*num_harmonics_]; // List for the variable names
        bool var_complex [2*num_harmonics_]; // List for indicating whether vars are complex
        for (int Q=0; Q<num_harmonics_; ++Q) { // Loop over harmonics to fill up lists
          std::ostringstream strs_s, strs_a;
          strs_s << "rho_s[" << Q << "]";
          strs_a << "rho_a[" << Q << "]";
          var_names[2*Q]   = strs_s.str(); // Use same order below
          var_names[2*Q+1] = strs_a.str();
          
          var_complex[2*Q]   = false;
          var_complex[2*Q+1] = false;
        }
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
        
        newDS.DefCoordVar(0, "T"); // Define "default" coordinate variables
        // Define "custom" coordinate variables
        const int varid_V1 = newDS.DefCoordVar_custom(1, "V1"); // Vary with g
        const int varid_V1p= newDS.DefCoordVar_custom(1, "V1p");
        const int varid_V2 = newDS.DefCoordVar_custom(1, "V2");
        const int varid_V3 = newDS.DefCoordVar_custom(1, "V3");
        
        newDS.EndDef(); // Exit definition mode
        
        newDS.WriteCoordVar(0, T_grid); // Write "default" coordinate variables
        newDS.WriteCoordVar_custom(varid_V1,  V1_grid); // Write "custom" coordinate variables
        newDS.WriteCoordVar_custom(varid_V1p, V1p_grid);
        newDS.WriteCoordVar_custom(varid_V2,  V2_grid);
        newDS.WriteCoordVar_custom(varid_V3,  V3_grid);
        
        double* vars [2*num_harmonics_]; // List for holding the pointers to the vars
        for (int Q=0; Q<num_harmonics_; ++Q) { // Fill up that list
          vars[2*Q]   = rho_s_grid[Q];
          vars[2*Q+1] = rho_a_grid[Q];
        }
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid); // Write loops variable
        newDS.WriteEnergy(energy_grid); // Write energy variable
        newDS.Writemu(mu_grid); // Write mu variable
    }
};

// Class that defines a parameter space for this Hamiltonian
class pspaceB_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceB_t(const pspaceB_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceB_t& operator=(const pspaceB_t&);
    
  public:
    // rho
    const size_t rho_pts = 3; const double rho_bounds [2] = {0.5, 0.7};
    // g: determines V1 and V1p
    const size_t g_pts = 20; const double V1_bounds [2] = {0.,40.}; const double V1p_bounds [2] = {0.,30.}; // Vary with g
    // h: determines V2 and V3
    const size_t h_pts = 15; const double V2_bounds [2] = {0.,6.}; const double V3_bounds  [2] = {0.,8.}; // Vary with h
    
    const int parspace_pts = rho_pts*g_pts*h_pts;
    
    
    // Coordinate variables
    double*const rho_grid; // rho coordinate variable
    double*const V1_grid; // V1 coordinate variable (depends on g)
    double*const V1p_grid; // V1p coordinate variable (depends on g)
    double*const V2_grid; // V2 coordinate variable (depends on h)
    double*const V3_grid; // V3 coordinate variable (depends on h)
    
    /* Variables for MF parameters. We store them in 1D arrays, but we must keep track of 
    the order of the coordinate variable "indices". */
    const size_t num_harmonics_; // We need the number of harmonics being considered.
    double*const*const rho_s_grid;
    double*const*const rho_a_grid;
    int*const loops_grid; // holds the number of loops done at each point
    double*const energy_grid; // Holds the MF energy for later comparison
    double*const mu_grid; // Holds the chemical potential for the converged system
    
    // Constructor declaration
    pspaceB_t(const int num_harmonics)
        :num_harmonics_(num_harmonics),
         rho_grid(new double [rho_pts]), // Coord vars
         V1_grid (new double [g_pts]), V1p_grid(new double [g_pts]),
         V2_grid (new double [h_pts]), V3_grid (new double [h_pts]),
         rho_s_grid (Alloc2D_d(num_harmonics, parspace_pts)),
         rho_a_grid (Alloc2D_d(num_harmonics, parspace_pts)), // Vars
         loops_grid (new int    [parspace_pts]),
         energy_grid(new double [parspace_pts]),
         mu_grid    (new double [parspace_pts]) // Important -- Note order
    {
        // The initialization list initializes the parameters and allocates memory for 
        // the arrays.
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(rho_bounds[0], rho_bounds[1], rho_pts, rho_grid, endpoint);
        
        LinInitArray(V1_bounds[0],  V1_bounds [1], g_pts, V1_grid,  endpoint);
        LinInitArray(V1p_bounds[0], V1p_bounds[1], g_pts, V1p_grid, endpoint);
        
        LinInitArray(V2_bounds[0], V2_bounds[1], h_pts, V2_grid, endpoint);
        LinInitArray(V3_bounds[0], V3_bounds[1], h_pts, V3_grid, endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(num_harmonics * parspace_pts, &(rho_s_grid[0][0]), -99.);//Init to unlikely value
        ValInitArray(num_harmonics * parspace_pts, &(rho_a_grid[0][0]), -99.);//Init to unlikely value
        
        ValInitArray(parspace_pts, loops_grid, 0);     // Initialize to 0
        ValInitArray(parspace_pts, energy_grid, -99.);
        ValInitArray(parspace_pts, mu_grid,      -9.);
        
        std::cout << "pspaceB_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceB_t()
    {
        delete [] rho_grid; // Coordinate variables
        delete [] V1_grid; delete [] V1p_grid;
        delete [] V2_grid; delete [] V3_grid;
        
        Dealloc2D(rho_s_grid); // MF variables
        Dealloc2D(rho_a_grid); // MF variables
        delete [] loops_grid; // Other variables
        delete [] energy_grid;
        delete [] mu_grid;
        
        std::cout << "pspaceB_t instance deleted.\n";
    }
    
    int Index(const int f, const int g, const int h)
    {
        // Returns the index of the var arrays corresponding to indices i, j of the coord
        // var arrays. We pick the order rho, V1 so that V1 varies faster than rho.
        return h + 
               h_pts*g + 
               h_pts*g_pts*f;
    }
    
    void SaveData(const string GlobalAttr, const string path)
    {
        // Method for saving the data of this class. This method uses the class from the 
        // module nc_IO that creates a simple NetCDF class and allows writing of variables.
        // We define parameters required to create the dataset. Don't forget to adjust these depending on the parameter space defined above. 
        // Important: Note that the order of the variables must be kept consistent.
        const size_t dims_num = 3;
        const string dim_names [dims_num] = {"rho", "g", "h"};
        const size_t dim_lengths [dims_num] = {rho_pts, g_pts, h_pts};
        const size_t vars_num = 2*num_harmonics_; // Variables other than coord variables
        
        string var_names [2*num_harmonics_]; // List for the variable names
        bool var_complex [2*num_harmonics_]; // List for indicating whether vars are complex
        for (int Q=0; Q<num_harmonics_; ++Q) { // Loop over harmonics to fill up lists
          std::ostringstream strs_s, strs_a;
          strs_s << "rho_s[" << Q << "]";
          strs_a << "rho_a[" << Q << "]";
          var_names[2*Q]   = strs_s.str(); // Use same order below
          var_names[2*Q+1] = strs_a.str();
          
          var_complex[2*Q]   = false;
          var_complex[2*Q+1] = false;
        }
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(dims_num, dim_names, dim_lengths, vars_num, var_names, var_complex, GlobalAttr, path);
        
        newDS.DefCoordVar(0, "rho"); // Define "default" coordinate variables
        // Define "custom" coordinate variables
        const int varid_V1 = newDS.DefCoordVar_custom(1, "V1"); // Vary with g
        const int varid_V1p= newDS.DefCoordVar_custom(1, "V1p");
        
        const int varid_V2 = newDS.DefCoordVar_custom(2, "V2"); // Vary with h
        const int varid_V3 = newDS.DefCoordVar_custom(2, "V3");
        
        newDS.EndDef(); // Exit definition mode
        
        newDS.WriteCoordVar(0, rho_grid); // Write "default" coordinate variables
        newDS.WriteCoordVar_custom(varid_V1,  V1_grid); // Write "custom" coordinate variables
        newDS.WriteCoordVar_custom(varid_V1p, V1p_grid);
        newDS.WriteCoordVar_custom(varid_V2,  V2_grid);
        newDS.WriteCoordVar_custom(varid_V3,  V3_grid);
        
        double* vars [2*num_harmonics_]; // List for holding the pointers to the vars
        for (int Q=0; Q<num_harmonics_; ++Q) { // Fill up that list
          vars[2*Q]   = rho_s_grid[Q];
          vars[2*Q+1] = rho_a_grid[Q];
        }
        newDS.WriteVars(vars); // Write the variables
        
        newDS.WriteLoops(loops_grid); // Write loops variable
        newDS.WriteEnergy(energy_grid); // Write energy variable
        newDS.Writemu(mu_grid); // Write mu variable
    }
};

// ######################################################################################

int pstudyA()
{
    // This routine performs the mean-field iterative search at every point in the 
    // parameter space defined above.
    
    const bool with_output = true; // Show output for diagnostics
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    // Declare and construct an instance of ham4_t
    ham4_t ham4(62,62,62); // Choose twice a prime number for grid resolution
    
    pspaceA_t pspaceA(ham4.num_harmonics); // Declare object of type pspaceA (parameter space)
    
    // Make any initial adjustment to the parameters
    ham4.set_nonzerotemp(1.e-3);
    //ham4.set_zerotemp();
    
    const string GlobalAttr = ham4.GetAttributes(); // assign attributes to GlobalAttr
    std::cout << "\n\n\ttol = " << ham4.tol_ << "\n";//Print tolerance for equality of MFs.
    
    // Loop over values of the parameter space
    for (int f=0; f<pspaceA.rho_pts; ++f)
      for (int g=0; g<pspaceA.g_pts; ++g)
      {
        if (with_output) // Print current params
          std::cout << "\n\nrho = " << pspaceA.rho_grid[f] << "; " << "g = " << g 
                    << " (V1 = " << pspaceA.V1_grid[g] << ", V1p = " << pspaceA.V1p_grid[g]
                    << ", V2 = " << pspaceA.V2_grid[g] << ", V3 = " << pspaceA.V3_grid[g] << ")\n";
        
        // Adjust phase space parameters
        ham4.assign_rho(pspaceA.rho_grid[f]); // Assign value of rho
        ham4.V1_  = pspaceA.V1_grid[g]; // g-dependent params
        ham4.V1p_ = pspaceA.V1p_grid[g];
        ham4.V2_  = pspaceA.V2_grid[g];
        ham4.V3_  = pspaceA.V3_grid[g];
        
        ham4.resetMFs(); // Resets MFs to default starting values.
        
        int loops=0; // Will be assigned the number of loops performed
        const bool fail = ham4.FixedPoint(&loops, with_output);
        
        
        if (fail) // Print current params
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "rho = " << pspaceA.rho_grid[f] << ", " << "g = " << g
                      << " (V1 = " << pspaceA.V1_grid[g] << ", V1p = " << pspaceA.V1p_grid[g]
                      << ", V2 = " << pspaceA.V2_grid[g] << ", V3 = "  << pspaceA.V3_grid[g]  << ")\n";
            ++numfails;
        }
        
        for (int Q=0; Q<ham4.num_harmonics; ++Q)
        { // We save the MF parameters to the pspaceA arrays.
          pspaceA.rho_s_grid[Q][pspaceA.Index(f,g)] = ham4.rho_s_[Q];
          pspaceA.rho_a_grid[Q][pspaceA.Index(f,g)] = ham4.rho_a_[Q];
        }
        pspaceA.loops_grid[pspaceA.Index(f,g)] = loops; // Save the number of loops to pspaceA array.
        pspaceA.energy_grid[pspaceA.Index(f,g)]= ham4.HFE_;
        
        if (with_output) std::cout << std::endl;
      }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO. Call saving method.
    pspaceA.SaveData(GlobalAttr, "data/ham4/"); // Include final '/' in path for saving
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int pstudyAA()
{
    // This routine performs the mean-field iterative search at every point in the 
    // parameter space defined above.
    
    const bool with_output = true; // Show output for diagnostics
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    // Declare and construct an instance of ham4_t
    ham4_t ham4(62,62,62); // Choose twice a prime number for grid resolution
    
    pspaceAA_t pspaceAA(ham4.num_harmonics); // Declare object of type pspaceAA (parameter space)
    
    // Make any initial adjustment to the parameters
    ham4.assign_rho(0.5); // Assign value of rho
    
    const string GlobalAttr = ham4.GetAttributes(); // assign attributes to GlobalAttr
    std::cout << "\n\n\ttol = " << ham4.tol_ << "\n";//Print tolerance for equality of MFs.
    
    // Loop over values of the parameter space
    for (int f=0; f<pspaceAA.T_pts; ++f)
      for (int g=0; g<pspaceAA.g_pts; ++g)
      {
        if (with_output) // Print current params
          std::cout << "\n\nT = " << pspaceAA.T_grid[f] << "; " << "g = " << g 
                    << " (V1 = " << pspaceAA.V1_grid[g] << ", V1p = " << pspaceAA.V1p_grid[g]
                    << ", V2 = " << pspaceAA.V2_grid[g] << ", V3 = " << pspaceAA.V3_grid[g] << ")\n";
        
        // Adjust phase space parameters
        ham4.set_nonzerotemp(pspaceAA.T_grid[f]); // Assign temperature
        ham4.V1_  = pspaceAA.V1_grid[g]; // g-dependent params
        ham4.V1p_ = pspaceAA.V1p_grid[g];
        ham4.V2_  = pspaceAA.V2_grid[g];
        ham4.V3_  = pspaceAA.V3_grid[g];
        
        ham4.resetMFs(); // Resets MFs to default starting values.
        
        int loops=0; // Will be assigned the number of loops performed
        const bool fail = ham4.FixedPoint(&loops, with_output);
        
        
        if (fail) // Print current params
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "T = " << pspaceAA.T_grid[f] << ", " << "g = " << g
                      << " (V1 = " << pspaceAA.V1_grid[g] << ", V1p = " << pspaceAA.V1p_grid[g]
                      << ", V2 = " << pspaceAA.V2_grid[g] << ", V3 = "  << pspaceAA.V3_grid[g]  << ")\n";
            ++numfails;
        }
        
        for (int Q=0; Q<ham4.num_harmonics; ++Q)
        { // We save the MF parameters to the pspaceAA arrays.
          pspaceAA.rho_s_grid[Q][pspaceAA.Index(f,g)] = ham4.rho_s_[Q];
          pspaceAA.rho_a_grid[Q][pspaceAA.Index(f,g)] = ham4.rho_a_[Q];
        }
        pspaceAA.loops_grid[pspaceAA.Index(f,g)] = loops; // Save the number of loops to pspaceAA array.
        pspaceAA.energy_grid[pspaceAA.Index(f,g)]= ham4.HFE_;
        
        if (with_output) std::cout << std::endl;
      }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO. Call saving method.
    pspaceAA.SaveData(GlobalAttr, "data/ham4/"); // Include final '/' in path for saving
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int pstudyB()
{
    // This routine performs the mean-field iterative search at every point in the 
    // parameter space defined above.
    
    const bool with_output = true; // Show output for diagnostics
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    // Declare and construct an instance of ham4_t
    ham4_t ham4(62,62,62); // Choose twice a prime number for grid resolution
    
    pspaceB_t pspaceB(ham4.num_harmonics); // Declare object of type pspaceB (parameter space)
    
    // Make any initial adjustment to the parameters
    ham4.set_nonzerotemp(1.e-3);
    
    const string GlobalAttr = ham4.GetAttributes(); // assign attributes to GlobalAttr
    std::cout << "\n\n\ttol = " << ham4.tol_ << "\n";//Print tolerance for equality of MFs.
    
    // Loop over values of the parameter space
    for (int f=0; f<pspaceB.rho_pts; ++f)
     for (int g=0; g<pspaceB.g_pts; ++g)
      for (int h=0; h<pspaceB.h_pts; ++h)
      {
        if (with_output) // Print current params
          std::cout << "\n\nrho = " << pspaceB.rho_grid[f] << "; "
                    << "g = " << g << " (V1 = " << pspaceB.V1_grid[g] << ", V1p = " << pspaceB.V1p_grid[g] << "); "
                    << "h = " << h << " (V2 = " << pspaceB.V2_grid[h] << ", V3 = " << pspaceB.V3_grid[h] << ")\n";
        
        // Adjust phase space parameters
        ham4.assign_rho(pspaceB.rho_grid[f]); // Assign value of rho
        ham4.V1_ = pspaceB.V1_grid[g];  ham4.V1p_ = pspaceB.V1p_grid[g]; // g-dependent params
        ham4.V2_ = pspaceB.V2_grid[h];  ham4.V3_ = pspaceB.V3_grid[h]; // h-dependent params
        
        ham4.resetMFs(); // Resets MFs to default starting values.
        
        int loops=0; // Will be assigned the number of loops performed
        const bool fail = ham4.FixedPoint(&loops, with_output);
        
        
        if (fail) // Print current params
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t" << "rho = " << pspaceB.rho_grid[f] << ", "
                      << "(V1 = " << pspaceB.V1_grid[g] << ", V1p = " << pspaceB.V1p_grid[g] << "), "
                      << "(V2 = " << pspaceB.V2_grid[h] << ", V3 = "  << pspaceB.V3_grid[h]  << ")\n";
            ++numfails;
        }
        
        for (int Q=0; Q<ham4.num_harmonics; ++Q)
        { // We save the MF parameters to the pspaceB arrays.
          pspaceB.rho_s_grid[Q][pspaceB.Index(f,g,h)] = ham4.rho_s_[Q];
          pspaceB.rho_a_grid[Q][pspaceB.Index(f,g,h)] = ham4.rho_a_[Q];
        }
        pspaceB.loops_grid[pspaceB.Index(f,g,h)] = loops; // Save the number of loops to pspaceB array.
        pspaceB.energy_grid[pspaceB.Index(f,g,h)]= ham4.HFE_;
        
        if (with_output) std::cout << std::endl;
      }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO. Call saving method.
    pspaceB.SaveData(GlobalAttr, "data/ham4/"); // Include final '/' in path for saving
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int main(int argc, char* argv[])
{
    //int info = pstudyA();
    
    int info = pstudyAA();
    
    //int info = pstudyB();
    
    return info;
}