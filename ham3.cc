// ham3.cc
/* Description */

#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <complex> // Needed to import alloc.h
#include "alloc.h" // Allocation/deallocation of arrays
#include "init_routines.h" // Initialization of arrays
#include "math_routines.h" // Various math functions
#include "chempot.h"
#include "kspace.h" // Defines a class for holding a band structure
#include "diag_routines.h" // Routines for finding evals and evecs
#include "ham3.h" // Include header file for consistency check
//using std::cos; using std::sin; using std::conj; // Not sure if these are necessary
using std::polar;
using std::to_string;
using std::complex;
using std::abs;

double Set_chi(const int counter, const int*const counter_vals, const double*const chi_vals, const int len, const bool with_output=false)
{
    /* When 'counter' reaches one of the values in 'counter_values', the corresponding 
    entry of 'chi_values' is returned. */
    double chi = 1.; // Starting value should be 1.
    for (int i=0; i<len; ++i)
        if (counter==counter_vals[i]) 
        {
            chi = chi_vals[i];
            if (with_output)
                std::cout << "\n\t Counter has reached " << counter
                          << ".\tchi = " << chi << "\n\n";
        }
    return chi;
}


/* Useful functions of momentum. */
double ham3_t::zeta(double ka, double kb)
{
    return 4.*cos(a_*ka/2.)*cos(a_*kb/2.);
}
complex<double> ham3_t::chi(double ka, double kb)
{
    return polar(1.,-a_*ka) + polar(1.,-a_*kb);
}
complex<double> ham3_t::f(double ka, double kb, double kc)
{
    const complex<double> ans( 4.*cos(a_*ka/2.)*cos(a_*kb/2.)*cos(c_*kc), 
                                    4.*sin(a_*ka/2.)*sin(a_*kb/2.)*sin(c_*kc) );
    return ans;
}
double ham3_t::ftilde(double ka, double kb, double kc)
{
    return 8. * cos(a_*ka/2.)*cos(a_*kb/2.)*cos(c_*kc);
}


/* Methods that set parameters with interdependencies. */
void ham3_t::assign_rho(const double rho)
{
    // Assign rho_ and dependent variables. Gets called in the constructor.
    rho_ = rho;
    filled_states = (int)( rho * (double)(num_states) );
}

void ham3_t::set_zerotemp()
{
    zerotemp_ = true; // Make sure this is true.
    T_ = -1.; // Irrelevant when zerotemp_ is true
}

void ham3_t::set_nonzerotemp(const double T)
{
    zerotemp_ = false; // Make sure this is false.
    T_ = T; // Set the private member to the right value
}


/* Constructor and destructor */
ham3_t::ham3_t(const int ka_pts, const int kb_pts, const int kc_pts)
    :ka_pts_(ka_pts), kb_pts_(kb_pts), kc_pts_(kc_pts)
{
    /* Constructor implementation */
    assign_rho(rho_); // Sets initial value of 'filled_states'
    std::cout << "ham3_t instance created.\n";
}

ham3_t::~ham3_t()
{
    /* Destructor implementation */
    std::cout << "ham3_t instance deleted.\n";
}

/* Assigns Hamiltonian */
void ham3_t::Assign_ham(const double ka, const double kb, const double kc, complex<double>*const*const ham_array)
{
    /* Given the arguments momentum as well as member parameters, calculate the 2*2
    k-space Hamiltonian and assign it to ham_array in full storage layout. The diag 
    routine only uses the lower triangle, so we only assign that part. */
    
    // For convenience, compute the A and B parts of the MFs
    const              double  rho_A = rho_   + rho_a_;
    const              double  rho_B = rho_   - rho_a_;
    const complex<double> u1p_A = u1p_s_ + u1p_a_;
    const complex<double> u1p_B = u1p_s_ - u1p_a_;
    const complex<double> u3_A  = u3_s_  + u3_a_;
    const complex<double> u3_B  = u3_s_  - u3_a_;
    
    // The "total" hopping amplitude includes the MF contribution
    const complex<double> t1_tot    = t1_ + V1_*u1_;
    const complex<double> t2_tot    = (t2A_ + V2_*u2A_) + (t2B_ + V2_*u2B_);
    const complex<double> t1p_A_tot = t1p_ + V1p_*u1p_A;
    const complex<double> t1p_B_tot = t1p_ + V1p_*u1p_B;
    const complex<double> t3_A_tot  = t3_ + V3_*u3_A;
    const complex<double> t3_B_tot  = t3_ + V3_*u3_B;
    
    
    complex<double>*const*const& H = ham_array; // For convenience, define reference
    
    // Initialize elements to zero
    H[0][0] = 0.;    H[0][1] = 0.;
    H[1][0] = 0.;    H[1][1] = 0.;
    
    // Add the density terms
    H[0][0] += 2.*(2.*V1p_+V3_)*rho_A + 4.*(V1_+2.*V2_)*rho_B;
    H[1][1] += 2.*(2.*V1p_+V3_)*rho_B + 4.*(V1_+2.*V2_)*rho_A;
    
    // Add the x- and y-direction hopping
    H[0][1] += - t1_tot * zeta(ka,kb);
    
    // Add the solenoid-like hopping
    H[0][1] += - t2_tot * f(ka,kb,kc);
    
    // Add the 'a'- and 'b'-direction hopping
    H[0][0] += -2. * std::real( t1p_A_tot * chi(ka,kb) );
    H[1][1] += -2. * std::real( t1p_B_tot * chi(ka,kb) );
    
    // Add the c-direction hopping
    H[0][0] += -2. * std::real( t3_A_tot * polar(1.,-c_*kc) );
    H[1][1] += -2. * std::real( t3_B_tot * polar(1.,-c_*kc) );
    
    
    // Assign the other off-diagonal element
    H[1][0] = std::conj(H[0][1]);
}

/* Methods for computing mean fields */
double ham3_t::ComputeTerm_rho_a(const double*const occs, const complex<double>*const*const evecs)
{
    // Evaluates the contribution to rho_a from a single k (see notes)
    const complex<double> mat[num_bands*num_bands] = {1., 0., // Square mat, row-major storage
                                                      0., -1.}; 
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace /= (double)(num_states); // Normalize appropriately (see notes)
    
    const double imag_part = std::imag(trace); // Check for zero imaginary part
    if ( abs(imag_part) > 1.e-15 )
      std::cerr << "WARNING: rho_a has nonzero imaginary part: " << imag_part<<std::endl;
    
    return std::real(trace);
}
complex<double> ham3_t::ComputeTerm_u1(const double ka, const double kb, const double*const occs, const complex<double>*const*const evecs)
{
    // Evaluates the contribution to u3_s from a single k (see notes)
    const complex<double> mat[num_bands*num_bands] = {0., 0., // Square mat, row-major storage
                                                      1., 0.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,a_*(ka+kb)/2.) 
           / (double)(ka_pts_*kb_pts_*kc_pts_); // Multiply by appropriate factors
    
    return trace;
}
complex<double> ham3_t::ComputeTerm_u1p_s(const double ka, const double*const occs, const complex<double>*const*const evecs)
{
    // Evaluates the contribution to u3_s from a single k (see notes)
    const complex<double> mat[num_bands*num_bands] = {1., 0., // Square mat, row-major storage
                                                      0., 1.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,a_*ka) / (double)(num_states); // Multiply by appropriate factors
    
    return trace;
}
complex<double> ham3_t::ComputeTerm_u1p_a(const double ka, const double*const occs, const complex<double>*const*const evecs)
{
    // Evaluates the contribution to u3_s from a single k (see notes)
    const complex<double> mat[num_bands*num_bands] = {1., 0., // Square mat, row-major storage
                                                      0., -1.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,a_*ka) / (double)(num_states); // Multiply by appropriate factors
    
    return trace;
}
complex<double> ham3_t::ComputeTerm_u2A(const double ka, const double kb, const double kc, const double*const occs, const complex<double>*const*const evecs)
{
    const complex<double> mat[num_bands*num_bands] = {0., 0., // Square mat, row-major storage
                                                      1., 0.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,-a_*(ka+kb)/2.+c_*kc) 
           / (double)(ka_pts_*kb_pts_*kc_pts_); // Multiply by appropriate factors
    
    return trace;
}
complex<double> ham3_t::ComputeTerm_u2B(const double ka, const double kb, const double kc, const double*const occs, const complex<double>*const*const evecs)
{
    const complex<double> mat[num_bands*num_bands] = {0., 0., // Square mat, row-major storage
                                                      1., 0.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,+a_*(ka+kb)/2.-c_*kc) 
           / (double)(ka_pts_*kb_pts_*kc_pts_); // Multiply by appropriate factors
    
    return trace;
}
complex<double> ham3_t::ComputeTerm_u3_s(const double kc, const double*const occs, const complex<double>*const*const evecs)
{
    // Evaluates the contribution to u3_s from a single k (see notes)
    const complex<double> mat[num_bands*num_bands] = {1., 0., // Square mat, row-major storage
                                                      0., 1.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,c_*kc) / (double)(num_states); // Multiply by appropriate factors
    
    return trace;
}
complex<double> ham3_t::ComputeTerm_u3_a(const double kc, const double*const occs, const complex<double>*const*const evecs)
{
    // Evaluates the contribution to u3_s from a single k (see notes)
    const complex<double> mat[num_bands*num_bands] = {1., 0., // Square mat, row-major storage
                                                      0., -1.};
    
    complex<double> trace = TraceMF(num_bands, evecs, mat, occs);
    trace *= polar(1.,c_*kc) / (double)(num_states); // Multiply by appropriate factors
    
    return trace;
}




void ham3_t::ComputeMFs(double& rho_a_out, complex<double>& u1_out,
                        complex<double>& u1p_s_out, complex<double>& u1p_a_out,
                        complex<double>& u2A_out, complex<double>& u2B_out,
                        complex<double>& u3_s_out, complex<double>& u3_a_out)
{
    // Declare (and construct) an instance of kspace_t.
    kspace_t kspace(a_, a_, c_, ka_pts_, kb_pts_, kc_pts_, num_bands);
    // array to hold Ham (local to thread)
    complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); //Initialize to zero
    
    // Given the parameters, diagonalize the Hamiltonian at each grid point
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          const int index = kspace.index(i, j, k, 0); // Last index is the band index
          simple_zheev(num_bands, &(ham_array[0][0]), &(kspace.energies[index]));
        }
    
    // Use all energies to compute chemical potential
    double mu = 666.;
    if (zerotemp_) // (elements get reordered)
        mu = FermiEnerg(num_states, filled_states, kspace.energies);
    else
        mu = ChemPotBisec(num_states, filled_states, kspace.energies, T_, 1.e-13);
    
    // Array to hold evecs (local to each thread)
    complex<double>*const*const evecs = Alloc2D_z(num_bands, num_bands);
    // Array to hold evals in second loop (local to thread)
    double*const evals = new double [num_bands];
    // Array to hold occupations (local to thread)
    double*const occs = new double [num_bands];
    ValInitArray(num_bands*num_bands, &(evecs[0][0])); // Initialize to zero
    ValInitArray(num_bands, evals); // Initialize to zero
    ValInitArray(num_bands, occs); // Initialize to zero
    
    // Use all the occupation numbers and the evecs to find the order parameter
    // Probably best to diagonalize a second time to avoid storing the evecs
            double  rho_a_accum = 0.;
    complex<double> u1_accum    = {0.,0.};
    complex<double> u1p_s_accum = {0.,0.};
    complex<double> u1p_a_accum = {0.,0.};
    complex<double> u2A_accum   = {0.,0.};
    complex<double> u2B_accum   = {0.,0.};
    complex<double> u3_s_accum  = {0.,0.};
    complex<double> u3_a_accum  = {0.,0.};
    
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          simple_zheev(num_bands, &(ham_array[0][0]), evals, true, &(evecs[0][0]));
          // Calculate occupations from energies, mu, and temperature
          Occupations(num_bands, mu, evals, occs, zerotemp_, T_);
          
          rho_a_accum += ComputeTerm_rho_a(occs, evecs);
          u1_accum    += ComputeTerm_u1(kspace.ka_grid[i], kspace.kb_grid[j], occs, evecs);
          u1p_s_accum += ComputeTerm_u1p_s(kspace.ka_grid[i], occs, evecs);
          u1p_a_accum += ComputeTerm_u1p_a(kspace.ka_grid[i], occs, evecs);
          u2A_accum   += ComputeTerm_u2A(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], occs, evecs);
          u2B_accum   += ComputeTerm_u2B(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], occs, evecs);
          u3_s_accum  += ComputeTerm_u3_s(kspace.kc_grid[k], occs, evecs);
          u3_a_accum  += ComputeTerm_u3_a(kspace.kc_grid[k], occs, evecs);
        }
    
    // This is where we should check that quantities are real, etc.
    
    rho_a_out = rho_a_accum;
    u1_out    = u1_accum;
    u1p_s_out = u1p_s_accum;
    u1p_a_out = u1p_a_accum;
    u2A_out   = u2A_accum;
    u2B_out   = u2B_accum;
    u3_s_out  = u3_s_accum;
    u3_a_out  = u3_a_accum;
    
    // Deallocate memory for arrays. The kspace_t destructor is called automatically.
    delete [] occs;
    delete [] evals;
    Dealloc2D(evecs);
    Dealloc2D(ham_array);
}


/*
std::string ham3_t::GetAttributes()
{
    // Define a string of metadata
    const std::string Attributes = "Haldane Hubbard bilayer (ham3)"
        ": NN distance a = "+to_string(a)+"; kx_pts = "+to_string(kx_pts)+
        "; kx_bounds = "+to_string(kx_bounds[0])+", "+to_string(kx_bounds[1])+
        "; ky_pts = " + to_string(ky_pts)+
        "; ky_bounds = "+to_string(ky_bounds[0])+", "+to_string(ky_bounds[1])+
        "; num_bands = "+to_string(num_bands)+
        "; parsI = ("+to_string(parsI.t1) +","+to_string(parsI.t2)+","+
                      to_string(parsI.eps)+","+to_string(parsI.phi)+")"
        "; parsII = ("+to_string(parsII.t1) +","+to_string(parsII.t2)+","+
                       to_string(parsII.eps)+","+to_string(parsII.phi)+")"
        "; tperp = "+to_string(tperp)+"; L = "+to_string(L)+"; rho = "+to_string(rho)+
        "; U = "+to_string(U)+"; rhoI_s_startval = "+to_string(rhoI_s_startval)+
        "; rhoI_a_startval = "+to_string(rhoI_a_startval)+
        "; mag_s_startval = "+to_string(mag_s_startval)+
        "; mag_a_startval = "+to_string(mag_a_startval)+"; tol = "+to_string(tol)+
        "; loops_lim_ = "+to_string(loops_lim_);
    
    return Attributes;
}
*/

bool ham3_t::FixedPoint(int*const num_loops_p, const bool with_output)
{
    // Performs the iterative self-consistent search using the parameters from ham3. 
    // The initial values of mag and rhoI are used as the starting values for the search; 
    // the end values are also output to mag and rhoI.
    if (with_output) // Format display output and set precision
        std::cout << std::scientific << std::showpos << std::setprecision(2);
    
    // Declare output variables and *initialize them to input values*.
            double  rho_a_out = rho_a_;
    complex<double> u1_out    = u1_;
    complex<double> u1p_s_out = u1p_s_;
    complex<double> u1p_a_out = u1p_a_;
    complex<double> u2A_out   = u2A_;
    complex<double> u2B_out   = u2B_;
    complex<double> u3_s_out  = u3_s_;
    complex<double> u3_a_out  = u3_a_;
    
    if (with_output)
      std::cout << "\t" "rho_a_" "\t" "u1_" "\t" "u1p_s_" "\t" "u1p_a_" "\t" "u2A_" "\t" "u2B_" "\t" "u3_s_" "\t" "u3_a_" << std::endl;
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        
        // Past a certain number of loops, we mix in part of the previous input vals
        const int len = 12;
        const int counter_vals [len] = {25,  50, 75, 100, 150, 200, 250, 300, 350, 400,  450,  500};
        const double  chi_vals [len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01};
        // Mixing fraction chi (chi=1 corresponds to using fully new value)
        const double chi = Set_chi(counter, counter_vals, chi_vals, len, with_output);
        
        rho_a_ = (1.-chi)*rho_a_ + chi*rho_a_out;
        u1_    = (1.-chi)*u1_    + chi*u1_out;
        u1p_s_ = (1.-chi)*u1p_s_ + chi*u1p_s_out;
        u1p_a_ = (1.-chi)*u1p_a_ + chi*u1p_a_out;
        u2A_   = (1.-chi)*u2A_   + chi*u2A_out;
        u2B_   = (1.-chi)*u2B_   + chi*u2B_out;
        u3_s_  = (1.-chi)*u3_s_  + chi*u3_s_out;
        u3_a_  = (1.-chi)*u3_a_  + chi*u3_a_out; // Update mean-field values
        
        
        if (with_output)
          std::cout << "\ninput\t"
                    << rho_a_ << "\t" << u1_ << "\t" << u1p_s_ << "\t" << u1p_a_ << "\t" 
                    << u2A_ << "\t" << u2B_ << "\t" << u3_s_ << "\t" << u3_a_ << std::endl;
        
        // Evaluate function. Output is assigned to the 'out' arguments.
        ComputeMFs(rho_a_out, u1_out, u1p_s_out, u1p_a_out, u2A_out, u2B_out, u3_s_out, u3_a_out);
        
                double  rho_a_diff = rho_a_out - rho_a_;
        complex<double> u1_diff    = u1_out    - u1_;
        complex<double> u1p_s_diff = u1p_s_out - u1p_s_;
        complex<double> u1p_a_diff = u1p_a_out - u1p_a_;
        complex<double> u2A_diff   = u2A_out   - u2A_;
        complex<double> u2B_diff   = u2B_out   - u2B_;
        complex<double> u3_s_diff  = u3_s_out  - u3_s_;
        complex<double> u3_a_diff  = u3_a_out  - u3_a_; // Differences between outputs and inputs
        
        if (with_output) // Print the differences
          std::cout << "diff\t"
                    << rho_a_diff << "\t" << u1_diff << "\t" << u1p_s_diff << "\t" << u1p_a_diff << "\t" 
                    << u2A_diff << "\t" << u2B_diff << "\t" << u3_s_diff << "\t" << u3_a_diff << std::endl;
        
        // Test for convergence
        converged = (abs(rho_a_diff)<tol_) && (abs(u1_diff)<tol_) && (abs(u1p_s_diff)<tol_) 
                 && (abs(u1p_a_diff)<tol_) && (abs(u2A_diff)<tol_) && (abs(u2B_diff)<tol_)
                 && (abs(u3_s_diff)<tol_) && (abs(u3_a_diff)<tol_);
        
        fail = (!converged) && (counter>loops_lim_); // Must come after converged line
        
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
