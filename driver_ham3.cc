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
#include "nc_IO.h" // Class for creating simple NetCDF datasets
#include "IO.h" // For PrintMatrix()
#include "ham3.h" // Source code for ham3


// Class that defines a parameter space for this Hamiltonian
class pspaceA_t {
  private:
    // Private copy constructor (prohibits copy creation)
    pspaceA_t(const pspaceA_t&);
    // Private assignment operator (prohibits assignment)
    const pspaceA_t& operator=(const pspaceA_t&);
    
  public:
    // scaling factor
    const int lambda_pts = 4; const double lambda_bounds [2] = {1./3., 1.};
    // Hubbard interaction strength
    const int U_pts = 6;     const double U_bounds [2] = {0., 10.};
    
    // Coordinate variables
    double*const lambda_grid; // lambda coordinate variable
    double*const U_grid; // U coordinate variable
    
    // MF parameters
    double*const*const rhoI_s_grid; // rhoI_s variable
    double*const*const rhoI_a_grid; // rhoI_a variable
    double*const*const mag_s_grid; // mag_s variable
    double*const*const mag_a_grid; // mag_a variable
    double*const*const loops_grid; // holds the number of loops done at each point
    
    // Constructor declaration
    pspaceA_t()
        :lambda_grid(new double [lambda_pts]), U_grid(new double [U_pts]),
         rhoI_s_grid(Alloc2D_d(lambda_pts, U_pts)),rhoI_a_grid(Alloc2D_d(lambda_pts, U_pts)),
         mag_s_grid(Alloc2D_d(lambda_pts, U_pts)), mag_a_grid(Alloc2D_d(lambda_pts, U_pts)),
         loops_grid(Alloc2D_d(lambda_pts, U_pts)) // Important -- Note order: lambda, U
    {
        /* The initialization list initializes the parameters and allocates memory for 
        the arrays. */
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(lambda_bounds[0], lambda_bounds[1], lambda_pts, lambda_grid, endpoint);
        LinInitArray(U_bounds[0],  U_bounds[1],  U_pts,  U_grid,  endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(lambda_pts*U_pts, &(rhoI_s_grid[0][0]), -99.);//Init to unlikely value
        ValInitArray(lambda_pts*U_pts, &(rhoI_a_grid[0][0]), -99.); // Likewise
        ValInitArray(lambda_pts*U_pts, &(mag_s_grid[0][0]), -99.); // Likewise
        ValInitArray(lambda_pts*U_pts, &(mag_a_grid[0][0]), -99.); // Likewise
        ValInitArray(lambda_pts*U_pts, &(loops_grid[0][0]), -99.); // Likewise
        
        
        std::cout << "pspaceA_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceA_t() {
        delete [] lambda_grid;
        delete [] U_grid;
        Dealloc2D(rhoI_s_grid);
        Dealloc2D(rhoI_a_grid);
        Dealloc2D(mag_s_grid);
        Dealloc2D(mag_a_grid);
        Dealloc2D(loops_grid);
        std::cout << "pspaceA_t instance deleted.\n";
    }
    
    void SaveData(const std::string GlobalAttr_, const std::string path_) {
        /* Method for saving the data of this class. This method uses the class from the 
        module nc_IO that creates a simple NetCDF class and allows writing of variables.*/
        /* We define parameters required to create the dataset. Don't forget to adjust 
        these depending on the parameter space defined above. 
        Important: Note that the order of the variables must be kept consistent. */
        const int dims_num = 2;
        const std::string dim_names [dims_num] = {"lambda", "U"};
        const int dim_lengths [dims_num] = {lambda_pts, U_pts};
        const int vars_num = 5; // Variables other than coord variables
        const std::string var_names [vars_num] = {"rho_s", "rho_a", "mag_s", "mag_a", 
                                                  "loops"}; // Use same order below
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(GlobalAttr_, dims_num, dim_names, dim_lengths,
                      vars_num, var_names, path_);
        
        const double*const coord_vars [dims_num] = {lambda_grid, U_grid};
        newDS.WriteCoordVars(coord_vars); // Write the coordinate variables
        
        const double*const vars [vars_num] = {&(rhoI_s_grid[0][0]), &(rhoI_a_grid[0][0]),
                                              &(mag_s_grid[0][0]), &(mag_a_grid[0][0]),
                                              &(loops_grid[0][0])};
        newDS.WriteVars(vars); // Write the variables
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
    // Energy bias between the layers
    const int L_pts = 4; const double L_bounds [2] = {0., 1.5};
    // Hubbard interaction strength
    const int U_pts = 4; const double U_bounds [2] = {0., 8.};
    
    // Coordinate variables
    double*const L_grid; // L coordinate variable
    double*const U_grid; // U coordinate variable
    
    // MF parameters
    double*const*const rhoI_s_grid; // rhoI_s variable
    double*const*const rhoI_a_grid; // rhoI_a variable
    double*const*const mag_s_grid; // mag_s variable
    double*const*const mag_a_grid; // mag_a variable
    double*const*const loops_grid; // holds the number of loops done at each point
    
    // Constructor declaration
    pspaceB_t()
        :L_grid(new double [L_pts]), U_grid(new double [U_pts]),
         rhoI_s_grid(Alloc2D_d(L_pts, U_pts)), rhoI_a_grid(Alloc2D_d(L_pts, U_pts)),
         mag_s_grid(Alloc2D_d(L_pts, U_pts)), mag_a_grid(Alloc2D_d(L_pts, U_pts)),
         loops_grid(Alloc2D_d(L_pts, U_pts)) // Important -- Note order: L, U
    {
        /* The initialization list initializes the parameters and allocates memory for 
        the arrays. */
        
        // Coordinate variables are initialized
        const bool endpoint = true; // Include endpoint (not an important choice)
        LinInitArray(L_bounds[0], L_bounds[1], L_pts, L_grid, endpoint);
        LinInitArray(U_bounds[0], U_bounds[1], U_pts, U_grid, endpoint);
        
        // The other variables (the order parameters) span the coordinate space
        ValInitArray(L_pts*U_pts, &(rhoI_s_grid[0][0]), -99.); // Init to unlikely value
        ValInitArray(L_pts*U_pts, &(rhoI_a_grid[0][0]), -99.); // Likewise
        ValInitArray(L_pts*U_pts, &(mag_s_grid[0][0]),  -99.); // Likewise
        ValInitArray(L_pts*U_pts, &(mag_a_grid[0][0]),  -99.); // Likewise
        ValInitArray(L_pts*U_pts, &(loops_grid[0][0]),  -99.); // Likewise
        
        
        std::cout << "pspaceB_t instance created.\n";
    }
    // Destructor declaration
    ~pspaceB_t() {
        delete [] L_grid;
        delete [] U_grid;
        Dealloc2D(rhoI_s_grid);
        Dealloc2D(rhoI_a_grid);
        Dealloc2D(mag_s_grid);
        Dealloc2D(mag_a_grid);
        Dealloc2D(loops_grid);
        std::cout << "pspaceB_t instance deleted.\n";
    }
    
    void SaveData(const std::string GlobalAttr_, const std::string path_) {
        /* Method for saving the data of this class. This method uses the class from the 
        module nc_IO that creates a simple NetCDF class and allows writing of variables.*/
        /* We define parameters required to create the dataset. Don't forget to adjust 
        these depending on the parameter space defined above. 
        Important: Note that the order of the variables must be kept consistent. */
        const int dims_num = 2;
        const std::string dim_names [dims_num] = {"L", "U"};
        const int dim_lengths [dims_num] = {L_pts, U_pts};
        const int vars_num = 5; // Variables other than coord variables
        const std::string var_names [vars_num] = {"rho_s", "rho_a", "mag_s", "mag_a", 
                                                  "loops"}; // Use same order below
        
        // Constructor for the dataset class creates a dataset
        newDS_t newDS(GlobalAttr_, dims_num, dim_names, dim_lengths,
                      vars_num, var_names, path_);
        
        const double*const coord_vars [dims_num] = {L_grid, U_grid};
        newDS.WriteCoordVars(coord_vars); // Write the coordinate variables
        
        const double*const vars [vars_num] = {&(rhoI_s_grid[0][0]), &(rhoI_a_grid[0][0]),
                                              &(mag_s_grid[0][0]), &(mag_a_grid[0][0]),
                                              &(loops_grid[0][0])};
        newDS.WriteVars(vars); // Write the variables
    }
};


// ######################################################################################
int pstudyA()
{
    /* This routine performs the mean-field iterative search at every point in the 
    parameter space defined above. */
    
    const bool with_output = false; // Show output for diagnostics
    
    std::cout << std::scientific << std::showpos; // Format display output
    
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    std::string GlobalAttr; // String to hold the metadata
    
    pspaceA_t pspaceA; // Declare object of type pspaceA (parameter space)
    
    // Loop over values of the parameter space
    /* PARALLELIZATION:, note that different threads do not write to the same parts of 
    pspaceA. */
    #pragma omp parallel default(none) shared(pspaceA,GlobalAttr,std::cout) reduction(+:numfails)
    {
    
    /* Declare and construct an instance of ham3_t (local to each thread). This works 
    because each instance of ham3 internally has its own arrays for the Hamiltonian, the 
    evals, and the evecs. */
    ham3_t ham3;
    
    #pragma omp single
    {
      GlobalAttr = ham3.GetAttributes();//Single thread assigns attributes to GlobalAttr
      std::cout << "\ntol = " << ham3.tol << "\n\n";//Print tolerance for equality of MFs.
    }
    
    
    #pragma omp for schedule(dynamic,1)
    for (int g=0; g<pspaceA.lambda_pts; ++g)
      for (int h=0; h<pspaceA.U_pts;  ++h)
      {
        if (with_output)
          std::cout << "lambda = " << pspaceA.lambda_grid[g] << ", "
                    << "U = "     << pspaceA.U_grid[h] << "\n"; // Print current params
        
        // Adjust phase space parameters
        ham3.parsI.SetScaling(pspaceA.lambda_grid[g]); // parsI is scaled by lambda
        ham3.U = pspaceA.U_grid[h];
        
        // Set the OPs to their startvals
        double rhoI_s = ham3.rhoI_s_startval;
        double rhoI_a = ham3.rhoI_a_startval;
        double mag_s = ham3.mag_s_startval;
        double mag_a = ham3.mag_a_startval;
        
        int loops=0; // Will receive the number of loops performed
        const bool fail = FixedPoint(rhoI_s, rhoI_a, mag_s, mag_a, ham3, &loops, 
                                                                            with_output);
        
        pspaceA.loops_grid[g][h] = loops; // Save the number of loops to pspaceA array.
        
        if (fail)
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "lambda = " << pspaceA.lambda_grid[g] << ", "
                      << "U = "     << pspaceA.U_grid[h] << "\n"; //Print current params
            ++numfails;
        }
        else
        {
            // We save the converged MF parameters to the pspaceA arrays.
            pspaceA.rhoI_s_grid[g][h] = rhoI_s;
            pspaceA.rhoI_a_grid[g][h] = rhoI_a;
            pspaceA.mag_s_grid[g][h] = mag_s;
            pspaceA.mag_a_grid[g][h] = mag_a;
        }
        
        if (with_output) std::cout << std::endl;
      }
    
    }
    
    
    // Print out MF parameters
    if (with_output)
    {
        std::cout << std::endl << "pspaceA.rhoI_s_grid = " << std::endl;
        PrintMatrix(pspaceA.lambda_pts, pspaceA.U_pts, pspaceA.rhoI_s_grid, std::cout);
        std::cout << std::endl << "pspaceA.rhoI_a_grid = " << std::endl;
        PrintMatrix(pspaceA.lambda_pts, pspaceA.U_pts, pspaceA.rhoI_a_grid, std::cout);
        std::cout << std::endl << "pspaceA.mag_s_grid = " << std::endl;
        PrintMatrix(pspaceA.lambda_pts, pspaceA.U_pts, pspaceA.mag_s_grid, std::cout);
        std::cout << std::endl << "pspaceA.mag_a_grid = " << std::endl;
        PrintMatrix(pspaceA.lambda_pts, pspaceA.U_pts, pspaceA.mag_a_grid, std::cout);
        std::cout << std::endl << "pspaceA.loops_grid = " << std::endl;
        PrintMatrix(pspaceA.lambda_pts, pspaceA.U_pts, pspaceA.loops_grid, std::cout);
        std::cout << std::endl;
    }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    const std::string path="data/ham3/"; //Choose the path for saving (include final '/')
    pspaceA.SaveData(GlobalAttr, path); // Call saving method
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int pstudyB()
{
    /* This routine performs the mean-field iterative search at every point in the 
    parameter space defined above. */
    
    const bool with_output = false; // Show output for diagnostics
    
    std::cout << std::scientific << std::showpos; // Format display output
    
    int numfails = 0; // Tracks number of points which failed to converge after loops_lim
    
    std::string GlobalAttr; // String to hold the metadata
    
    pspaceB_t pspaceB; // Declare object of type pspaceB (parameter space)
    
    // Loop over values of the parameter space
    /* PARALLELIZATION:, note that different threads do not write to the same parts of 
    pspaceB. */
    #pragma omp parallel default(none) shared(pspaceB,GlobalAttr,std::cout) reduction(+:numfails)
    {
    
    /* Declare and construct an instance of ham3_t (local to each thread). This works 
    because each instance of ham3 internally has its own arrays for the Hamiltonian, the 
    evals, and the evecs. */
    ham3_t ham3;
    
    #pragma omp single
    {
      GlobalAttr = ham3.GetAttributes();//Single thread assigns attributes to GlobalAttr
      std::cout << "\ntol = " << ham3.tol << "\n\n";//Print tolerance for equality of MFs.
    }
    
    ham3.parsI.SetScaling(ham3.lambda); // We manually set the parsI scaling
    
    
    #pragma omp for schedule(dynamic,1)
    for (int g=0; g<pspaceB.L_pts; ++g)
      for (int h=0; h<pspaceB.U_pts;  ++h)
      {
        if (with_output)
          std::cout << "L = " << pspaceB.L_grid[g] << ", "
                    << "U = " << pspaceB.U_grid[h] << "\n"; // Print current params
        
        // Adjust phase space parameters
        ham3.L = pspaceB.L_grid[g]; // L value is set
        ham3.U = pspaceB.U_grid[h] * ham3.lambda; // SCALED U value is set
        
        // Set the OPs to their startvals
        double rhoI_s = ham3.rhoI_s_startval;
        double rhoI_a = ham3.rhoI_a_startval;
        double mag_s = ham3.mag_s_startval;
        double mag_a = ham3.mag_a_startval;
        
        int loops=0; // Will receive the number of loops performed
        const bool fail = FixedPoint(rhoI_s, rhoI_a, mag_s, mag_a, ham3, &loops, 
                                                                            with_output);
        
        pspaceB.loops_grid[g][h] = loops; // Save the number of loops to pspaceB array.
        
        if (fail)
        {
            std::cout << "\tWARNING: failure to converge after limit reached\t"
                      << "L = " << pspaceB.L_grid[g] << ", "
                      << "U = " << pspaceB.U_grid[h] << "\n"; //Print current params
            ++numfails;
        }
        else
        {
            // We save the converged MF parameters to the pspaceB arrays.
            pspaceB.rhoI_s_grid[g][h] = rhoI_s;
            pspaceB.rhoI_a_grid[g][h] = rhoI_a;
            pspaceB.mag_s_grid[g][h] = mag_s;
            pspaceB.mag_a_grid[g][h] = mag_a;
        }
        
        if (with_output) std::cout << std::endl;
      }
    
    }
    
    
    // Print out MF parameters
    if (with_output)
    {
        std::cout << std::endl << "pspaceB.rhoI_s_grid = " << std::endl;
        PrintMatrix(pspaceB.L_pts, pspaceB.U_pts, pspaceB.rhoI_s_grid, std::cout);
        std::cout << std::endl << "pspaceB.rhoI_a_grid = " << std::endl;
        PrintMatrix(pspaceB.L_pts, pspaceB.U_pts, pspaceB.rhoI_a_grid, std::cout);
        std::cout << std::endl << "pspaceB.mag_s_grid = " << std::endl;
        PrintMatrix(pspaceB.L_pts, pspaceB.U_pts, pspaceB.mag_s_grid, std::cout);
        std::cout << std::endl << "pspaceB.mag_a_grid = " << std::endl;
        PrintMatrix(pspaceB.L_pts, pspaceB.U_pts, pspaceB.mag_a_grid, std::cout);
        std::cout << std::endl << "pspaceB.loops_grid = " << std::endl;
        PrintMatrix(pspaceB.L_pts, pspaceB.U_pts, pspaceB.loops_grid, std::cout);
        std::cout << std::endl;
    }
    
    
    // We save to a NetCDF dataset using the class defined in nc_IO
    const std::string path="data/ham3/"; //Choose the path for saving (include final '/')
    pspaceB.SaveData(GlobalAttr, path); // Call saving method
    
    // Print out numfails
    std::cout << "\n\tNUMBER OF FAILURES TO CONVERGE: " << numfails << "\n\n";
    
    return numfails;
}

int main(int argc, char* argv[])
{
    //int info = pstudyA();
    int info = pstudyB();
    
    return info;
}