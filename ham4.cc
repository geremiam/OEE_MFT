// ham4.cc
/* Description */

#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <sstream> // For stringstreams
#include <complex> // Needed to import alloc.h
#include "alloc.h" // Allocation/deallocation of arrays
#include "init_routines.h" // Initialization of arrays
#include "math_routines.h" // Various math functions
#include "misc_routines.h"
#include "chempot.h"
#include "kspace.h" // Defines a class for holding a band structure
#include "diag_routines.h" // Routines for finding evals and evecs
#include "ham4.h" // Include header file for consistency check
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
double ham4_t::zeta(double ka, double kb) const
{
    return 4.*cos(a_*ka/2.)*cos(a_*kb/2.);
}
complex<double> ham4_t::f(double ka, double kb, double kc) const
{
    const complex<double> ans( 4.*cos(a_*ka/2.)*cos(a_*kb/2.)*cos(c_*kc), 
                               4.*sin(a_*ka/2.)*sin(a_*kb/2.)*sin(c_*kc) );
    return ans;
}


int ham4_t::idx(const int Q, const int alpha) const
{
    // Gives the appropriate index in the direct product between Q and alpha.
    // The direct product is defined so that alpha varies faster than Q.
    return idx_composite(num_harmonics, states_per_cell, Q, alpha);
}


/* Methods that set parameters with interdependencies. */
void ham4_t::assign_rho(const double rho)
{
    // Assign rho_ and dependent variables. Gets called in the constructor.
    rho_ = rho;
    filled_states = (int)( rho * (double)(num_states) );
}

void ham4_t::set_zerotemp()
{
    zerotemp_ = true; // Make sure this is true.
    T_ = -1.; // Irrelevant when zerotemp_ is true
}

void ham4_t::set_nonzerotemp(const double T)
{
    zerotemp_ = false; // Make sure this is false.
    T_ = T; // Set the private member to the right value
}


/* Constructor and destructor */
ham4_t::ham4_t(const int ka_pts, const int kb_pts, const int kc_pts)
    :ka_pts_(ka_pts), kb_pts_(kb_pts), kc_pts_(kc_pts)
{
    /* Constructor implementation */
    assign_rho(rho_); // Sets initial value of 'filled_states'
    //resetMFs(); // Sets the mean fields to default value
    std::cout << "ham4_t instance created.\n";
}

ham4_t::~ham4_t()
{
    /* Destructor implementation */
    std::cout << "ham4_t instance deleted.\n";
}

// Methods for assigning the Hamiltonian
void ham4_t::Assign_h(double ka, double kb, double kc, complex<double>*const h) const
{
    // With this convention, it is not necessary to use momenta inside the FBZ.
    //ka = FBZ(ka, a_);
    //kb = FBZ(kb, a_);
    //kc = FBZ(kc, c_);
    // Assign to h
    h[0] = -2.*t1p_ * (cos(a_*ka) + cos(a_*kb)) - 2.*t3_ * cos(c_*kc);
    h[1] = //(- t1_ * zeta(ka,kb) - t2A_ * f(ka,kb,kc) - t2B_ * std::conj(f(ka,kb,kc))) * polar(1.,-a_*(ka+kb)/2.);
           - t1_  * ( 1. + polar(1.,-a_*ka) + polar(1.,-a_*kb) + polar(1.,-a_*ka - a_*kb) )
           - t2A_ * ( polar(1.,-c_*kc) + polar(1.,+c_*kc - a_*ka) + polar(1.,+c_*kc - a_*kb) + polar(1.,-a_*ka - a_*kb - c_*kc) )
           - t2B_ * ( polar(1.,+c_*kc) + polar(1.,-c_*kc - a_*ka) + polar(1.,-c_*kc - a_*kb) + polar(1.,-a_*ka - a_*kb + c_*kc) );
    h[2] = std::conj(h[1]);
    h[3] = h[0];
}
void ham4_t::Assign_V(const int Q, complex<double>*const V) const
{
    // Assigns the matrix to V with qa = Qa[Q], qb = Qb[Q], etc.
    const double qa = Qa[Q];
    const double qb = Qb[Q];
    const double qc = Qc[Q];
    // Assign to V
    V[0] = 2.*V1p_ * (cos(a_*qa) + cos(a_*qb)) + 2.*V3_ * cos(c_*qc);
    V[1] = (V1_ + 2.*V2_*cos(c_*qc)) * (1. + polar(1.,-a_*qa) + polar(1.,-a_*qb) + polar(1.,-a_*qa-a_*qb));
    V[2] = std::conj(V[1]);
    V[3] = V[0];
}
void ham4_t::Assign_V_manual(const int Q, complex<double>*const V) const
{
    // Manually computed assignment formulas for V. Assigns the matrix to V with 
    // qa = Qa[Q], qb = Qb[Q], etc.
    // NEEDS TO BE REWRITTEN FOR DIFFERENT SETS OF HARMONICS
    if (Q==0) {
      V[0] = 4.*V1p_ + 2.*V3_;
      V[1] = 4.*(V1_ + 2.*V2_);
    }
    
    else if (Q==1) {
      V[0] = 2.*V3_;
      V[1] = 0.;
    }
    
    else if (Q==2) {
      V[0] = 2.*V3_;
      V[1] = 0.;
    }
    
    else if (Q==3) {
      V[0] = -4.*V1p_ + 2.*V3_;
      V[1] = 0.;
    }
    else
      std::cout << "WARNING: In routine Assign_V_manual, index Q is out of range." << std::endl;
    
    V[2] = std::conj(V[1]);
    V[3] = V[0];
}

// Assigns Hamiltonian
void ham4_t::Assign_ham(const double ka, const double kb, const double kc, complex<double>*const*const ham_array) const
{
    /* Given the arguments momentum as well as member parameters, calculate the 2*2
    k-space Hamiltonian and assign it to ham_array in full storage layout. The diag 
    routine only uses the lower triangle, but we assign the whole thing. */
    
    complex<double>*const*const& H = ham_array; // For convenience, define reference
    // Initialize elements to zero
    for (int i=0; i<num_bands*num_bands; ++i)
        H[0][i] = 0.;
    // **********************************************************************************
    
    // Add the kinetic terms
    
    // Array for holding the 2*2 kinetic part of the Hamiltonian
    complex<double> h[states_per_cell*states_per_cell]={0.};
    for (int Q=0; Q<num_harmonics; ++Q)
    {
      Assign_h(ka+Qa[Q], kb+Qb[Q], kc+Qc[Q], h); // Assign the 2*2 kinetic part of the Hamiltonian
      H[idx(Q,0)][idx(Q,0)] += h[0];
      H[idx(Q,0)][idx(Q,1)] += h[1];
      H[idx(Q,1)][idx(Q,0)] += h[2];
      H[idx(Q,1)][idx(Q,1)] += h[3];
    }
    // **********************************************************************************
    
    
    // Declare an array rhotilde for each alpha. These have one index for the Harmonic.
    complex<double> rhotilde_A [num_harmonics]={0.};
    complex<double> rhotilde_B [num_harmonics]={0.};
    // Declare the array V, to be evaluated for every Q, needed to compte rhotilde.
    complex<double> V[states_per_cell*states_per_cell]={0.};
    for (int Q=0; Q<num_harmonics; ++Q)
    {
        // For convenience, assign 
        const double rho_A = rho_s_[Q] + rho_a_ [Q];
        const double rho_B = rho_s_[Q] - rho_a_ [Q];
        Assign_V(Q, V);// The matrix V, which depends on Q, gets assigned
        rhotilde_A[Q] = V[0]*rho_A + V[1]*rho_B;// V[0,0]*rho_A_[Q] + V[0,1]*rho_B_[Q];
        rhotilde_B[Q] = V[2]*rho_A + V[3]*rho_B;// V[1,0]*rho_A_[Q] + V[1,1]*rho_B_[Q];
    }
    
    // Add the density terms to the Hamiltonian
    for (int P=0; P<num_harmonics; ++P)
      for (int Q=0; Q<num_harmonics; ++Q)
      {
          const int Q_sum = addition_table[P][Q];
          // Loop over the momentum double sum; use the addition table.
          H[idx(P,0)][idx(Q_sum,0)] += std::conj(rhotilde_A[Q]);
          H[idx(P,1)][idx(Q_sum,1)] += std::conj(rhotilde_B[Q]);
      }
}


// Methods for computing mean fields
void ham4_t::AddContribution_rho(const double*const occs, const complex<double>*const*const evecs, double*const rho_A, double*const rho_B) const
{
    // Evaluates the contribution to rho_A and rho_B from a single k (see notes) and adds 
    // it to the arrays "rho_A" and "rho_B".
    
    for (int Q=0; Q<num_harmonics; ++Q) // Loop over the harmonics i.e. the components of rho_A and rho_B
    {
      complex<double> temp_A = {0.,0.}; // These will be used in the summation.
      complex<double> temp_B = {0.,0.}; // Initialize to zero
      
      for (int beta=0; beta<states_per_cell; ++beta) // summation loop
        for (int P=0; P<num_harmonics; ++P)          // summation loop
          for (int R=0; R<num_harmonics; ++R)        // summation loop
          {
            const int QplusR = addition_table[Q][R];
            
            temp_A += conj(evecs[idx(R,0)][idx(P,beta)]) * evecs[idx(QplusR,0)][idx(P,beta)] * occs[idx(P,beta)];
            temp_B += conj(evecs[idx(R,1)][idx(P,beta)]) * evecs[idx(QplusR,1)][idx(P,beta)] * occs[idx(P,beta)];
          }
      if (abs(std::imag(temp_A))>1.e-16)
        std::cout << "WARNING: nonzero imaginary part: temp_A = " << temp_A << std::endl;
      if (abs(std::imag(temp_B))>1.e-16)
        std::cout << "WARNING: nonzero imaginary part: temp_B = " << temp_B << std::endl;
      rho_A[Q] += std::real(temp_A) / (double)(num_unit_cells);
      rho_B[Q] += std::real(temp_B) / (double)(num_unit_cells);
    }
}


double ham4_t::ComputeMFs(double*const rho_s_out, double*const rho_a_out) const
{
    // Declare (and construct) an instance of kspace_t.
    kspace_t kspace(a_, a_, c_, ka_pts_, kb_pts_, kc_pts_, num_bands);
    
    // Step 1: diagonalize to find all the energy evals and store them in kspace
    #pragma omp parallel default(none) shared(kspace)
    {
    // array to hold Ham (local to thread)
    complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols);
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); //Initialize to zero
    #pragma omp for collapse(3)
    // Given the parameters, diagonalize the Hamiltonian at each grid point
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          const int index = kspace.index(i, j, k, 0); // Last index is the band index
          simple_zheev(num_bands, &(ham_array[0][0]), &(kspace.energies[index]));
        }
    Dealloc2D(ham_array);
    }
    
    // Step 2: Use all energies to compute chemical potential
    double mu = 666.;
    if (zerotemp_) // (ELEMENTS GET REORDERED)
        mu = FermiEnerg(num_states, filled_states, kspace.energies);
    else // USES OPENMP PARALLELIZATION
    {
        const bool show_output = false;
        const bool usethreads = true;
        mu = ChemPotBisec(num_states, filled_states, kspace.energies, T_, show_output, usethreads);
    }
    
    std::cout << "\tCheckpoint 1" << std::endl;
    // Step 3: Use all the occupation numbers and the evecs to find the order parameter
    // Probably best to diagonalize a second time to avoid storing the evecs
    double rho_A_accum [num_harmonics] = {0.}; // IMPORTANT: MUST BE INITIALIZED TO ZERO
    double rho_B_accum [num_harmonics] = {0.}; // IMPORTANT: MUST BE INITIALIZED TO ZERO
    
    #pragma omp parallel default(none) firstprivate(mu) shared(kspace, std::cout) reduction(+:rho_A_accum,rho_B_accum)
    {
    complex<double>*const*const ham_array = Alloc2D_z(ham_array_rows, ham_array_cols); // array to hold Ham (local to thread)
    complex<double>*const*const     evecs = Alloc2D_z(num_bands, num_bands); // Array to hold evecs (local to each thread)
    double*const evals = new double [num_bands]; // Array to hold evals in second loop (local to thread)
    double*const  occs = new double [num_bands]; // Array to hold occupations (local to thread)
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0])); // Initialize to zero
    ValInitArray(num_bands*num_bands, &(evecs[0][0])); // Initialize to zero
    ValInitArray(num_bands, evals); // Initialize to zero
    ValInitArray(num_bands, occs); // Initialize to zero
    
    #pragma omp for collapse(3)
    for (int i=0; i<ka_pts_; ++i)
      for (int j=0; j<kb_pts_; ++j)
        for (int k=0; k<kc_pts_; ++k)
        {
          //std::cout << "i,j,k = " << i << "," << j << "," << k << std::endl;
          Assign_ham(kspace.ka_grid[i], kspace.kb_grid[j], kspace.kc_grid[k], ham_array);
          simple_zheev(num_bands, &(ham_array[0][0]), evals, true, &(evecs[0][0]));
          // Calculate occupations from energies, mu, and temperature
          Occupations(num_bands, mu, evals, occs, zerotemp_, T_);
          
          AddContribution_rho(occs, evecs, rho_A_accum, rho_B_accum);
        }
    delete [] occs; // Deallocate memory for arrays. 
    delete [] evals;
    Dealloc2D(evecs);
    Dealloc2D(ham_array);
    }
    
    // Assign the results to the output arrays
    
    for (int Q=0; Q<num_harmonics; ++Q)
    {
        rho_s_out[Q] = 0.5 * (rho_A_accum[Q] + rho_B_accum[Q]);
        rho_a_out[Q] = 0.5 * (rho_A_accum[Q] - rho_B_accum[Q]);
    }
    
    // The kspace_t destructor is called automatically.
    
    return 0.; // Should return the Helmholtz free energy
}


/*
std::string ham4_t::GetAttributes()
{
    // Define a string of metadata
    // Watch out: this method resets the values of the MFs to their starting values.
    resetMFs();
        
    std::ostringstream strs; // Declare a stringstream to which to write the attributes
    strs << "Gyrotropic SSB (ham4)"
         << ": a = " << a_ << ", c = " << c_
         << "; k-space points: " << ka_pts_ << ", " << kb_pts_ << ", " << kc_pts_
         << "; num_bands = " << num_bands << "; rho = " << rho_
         << "; zerotemp = " << zerotemp_ << ", T = " << T_ << "; tol = " << tol_
         << "; t1 = " << t1_ << ", t1p = " << t1p_ << ", t2A = " << t2A_ << "; t2B = " << t2B_ << ", t3 = " << t3_ 
         << ", V1 = " << V1_ << "; V1p = " << V1p_ << ", V2 = " << V2_ << ", V3 = " << V3_
         << "; loops_lim = " << loops_lim_
         << "; Starting values of the MFs: rho_a_ = " << rho_a_ << ", u1_ = " << u1_
         << ", u1p_s_ = " << u1p_s_ << ", u1p_a_ = " << u1p_a_
         << ", u2A_ = "   << u2A_   << ", u2B_ = "   << u2B_
         << ", u3_s_ = "  << u3_s_  << ", u3_a_ = "  << u3_a_;
    
    const std::string Attributes = strs.str(); // Write the stream contents to a string
    
    return Attributes;
}
*/

bool ham4_t::FixedPoint(int*const num_loops_p, const bool with_output)
{
    // Performs the iterative self-consistent search using the parameters from ham4. 
    // The initial values of the MF attributes are used as the starting values for the 
    // search; the end values are also stored in these same MF attributes.
    if (with_output) // Format display output and set precision
        std::cout << std::scientific << std::showpos << std::setprecision(2);
    
    // Declare output variables and INITIALIZE THEM TO INPUT VALUES
    double rho_s_out [num_harmonics];
    double rho_a_out [num_harmonics];
    copy_array(num_harmonics, rho_s_, rho_s_out); // Copy rho_s_ into rho_s_out
    copy_array(num_harmonics, rho_a_, rho_a_out); // Copy rho_a_ into rho_a_out
    double HFE_prev = 0.; // For keeping track of free energy at previous step
    
    if (with_output)
      std::cout << "\t" "rho_s_" "\t" "rho_a_" << std::endl;
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        
        // Past a certain number of loops, we mix in part of the previous input vals
        const int len = 16;
        const int counter_vals [len] = {50,  100, 150, 200, 250, 300, 450, 500, 550, 600,  650,  700,  800,   900,   1000,  1100};
        const double  chi_vals [len] = {0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.03, 0.01, 5.e-3, 1.e-3, 5.e-4, 1.e-5};
        // Mixing fraction chi (chi=1 corresponds to using fully new value)
        const double chi = Set_chi(counter, counter_vals, chi_vals, len, with_output);
        
        for (int Q=0; Q<num_harmonics; ++Q) // Update mean-field values
        {
            rho_s_[Q] = (1.-chi)*rho_s_[Q] + chi*rho_s_out[Q];
            rho_a_[Q] = (1.-chi)*rho_a_[Q] + chi*rho_a_out[Q];
        }
        HFE_prev = HFE_; // Store previous free energy in HFE_prev
        
        
        if (with_output)
        {
          std::cout << "\ninput\t";
          for (int Q=0; Q<num_harmonics; ++Q)
            std::cout << rho_s_[Q] << "\t";
          for (int Q=0; Q<num_harmonics; ++Q)
            std::cout << rho_a_[Q] << "\t";
          std::cout << std::endl;
        }
        
        // Compute MFs. Output is assigned to the 'out' arguments.
        // The HFE for the final set of MF values is left in the attribute HFE_.
        HFE_ = ComputeMFs(rho_s_out, rho_a_out);
        
        double rho_s_diff [num_harmonics]; // Declare arrays to store the changes in MF
        double rho_a_diff [num_harmonics];
        for (int Q=0; Q<num_harmonics; ++Q) // Differences between outputs and inputs
        {
            rho_s_diff[Q] = rho_s_out[Q] - rho_s_[Q];
            rho_a_diff[Q] = rho_a_out[Q] - rho_a_[Q];
        }
        
        if (with_output) // Print the differences
        {
          std::cout << "diff\t";
          for (int Q=0; Q<num_harmonics; ++Q)
            std::cout << rho_s_diff[Q] << "\t";
          for (int Q=0; Q<num_harmonics; ++Q)
            std::cout << rho_a_diff[Q] << "\t";
          std::cout << std::endl;
        }
        
        // Test for convergence
        converged = check_bound_array(num_harmonics, tol_, rho_s_diff) && check_bound_array(num_harmonics, tol_, rho_a_diff);
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

/*

// **************************************************************************************
// ROUTINES FOR CALCULATING THE FREE ENERGY
double ham4_t::Helmholtz(const double*const energies, const double mu,
                         const double rho_a_out, const complex<double> u1_out, const complex<double> u1p_s_out, const complex<double> u1p_a_out,
                         const complex<double> u2A_out, const complex<double> u2B_out, const complex<double> u3_s_out, const complex<double> u3_a_out) const
{
    // The Helmholtz FE normalized according to the NUMBER OF ATOMS
    // Should be use to compare states of a fixed particle number
    // Second term is + mu N (properly normalized).
    const double ans = Omega_trial(energies,mu,rho_a_out,u1_out,u1p_s_out,u1p_a_out,u2A_out,u2B_out,u3_s_out,u3_a_out) + mu*rho_;
    return ans;
}

double ham4_t::Omega_trial(const double*const energies, const double mu,
                           const double rho_a_out, const complex<double> u1_out, const complex<double> u1p_s_out, const complex<double> u1p_a_out,
                           const complex<double> u2A_out, const complex<double> u2B_out, const complex<double> u3_s_out, const complex<double> u3_a_out) const
{
    // The (GC) free energy normalized according to the NUMBER OF ATOMS
    // Each component is already normalized.
    const double ans = Omega_MF(energies, mu)
                     + mean_Hint(rho_a_out,u1_out,u1p_s_out,u1p_a_out,u2A_out,u2B_out,u3_s_out,u3_a_out)
                     - mean_Hint_MF(rho_a_out,u1_out,u1p_s_out,u1p_a_out,u2A_out,u2B_out,u3_s_out,u3_a_out);
    return ans;
}

double ham4_t::Omega_MF(const double*const energies, const double mu) const
{
    // (GC) free energy of the MF Hamiltonian. 
    // Normalized according to the NUMBER OF ATOMS
    // It's calculated differently at zero and nonzero temperatures.
    double accumulator = 0.;
    
    if (zerotemp_)
    {
      for (int i=0; i<num_states; ++i)
        if (energies[i]<=mu)
          accumulator += (energies[i] - mu);
    }
    
    else
      for (int i=0; i<num_states; ++i)
        accumulator += - T_ * log_1p_exp(-(energies[i] - mu)/T_);
    
    // Normalize by the number of atoms
    accumulator /= (double)(2*ka_pts_*kb_pts_*kc_pts_);
    
    return accumulator;
}

double ham4_t::mean_Hint(const double rho_a_out, const complex<double> u1_out, const complex<double> u1p_s_out, const complex<double> u1p_a_out,
                         const complex<double> u2A_out, const complex<double> u2B_out, const complex<double> u3_s_out, const complex<double> u3_a_out) const
{
    // Expectation value of the full interaction term in the Hamiltonian
    // Normalized according to the NUMBER OF ATOMS
    
    // For convenience, define
    const double rho_A_out = rho_ + rho_a_out;
    const double rho_B_out = rho_ - rho_a_out;
    
    const complex<double> u1p_A_out = u1p_s_out + u1p_a_out;
    const complex<double> u1p_B_out = u1p_s_out - u1p_a_out;
    
    const complex<double> u3_A_out = u3_s_out + u3_a_out;
    const complex<double> u3_B_out = u3_s_out - u3_a_out;
    
    double ans = 0.;
    
    // NOTE: norm() is the MODULUS SQUARED, for some reason.
    ans += 0.5*(V3_+2.*V1p_) * ( std::norm(rho_A_out) + std::norm(rho_B_out) ); // Density terms
    ans += 2. *(V1_+2.*V2_)  * rho_A_out * rho_B_out;
    
    ans += - 2.*V1_ * std::norm(u1_out);
    ans += - 2.*V2_ * ( std::norm(u2A_out) + std::norm(u2B_out) );
    ans += -    V1p_* ( std::norm(u1p_A_out) + std::norm(u1p_B_out) );
    ans += - 0.5*V3_* ( std::norm(u3_A_out)  + std::norm(u3_B_out)  );
    
    return ans;
}

double ham4_t::mean_Hint_MF(const double rho_a_out, const complex<double> u1_out, const complex<double> u1p_s_out, const complex<double> u1p_a_out,
                            const complex<double> u2A_out, const complex<double> u2B_out, const complex<double> u3_s_out, const complex<double> u3_a_out) const
{
    // Expectation value of the M-F interaction term in the Hamiltonian.
    // Normalized according to the NUMBER OF ATOMS
    
    // For convenience, define
    const double rho_A     = rho_ + rho_a_;
    const double rho_B     = rho_ - rho_a_;
    const double rho_A_out = rho_ + rho_a_out;
    const double rho_B_out = rho_ - rho_a_out;
    
    const complex<double> u1p_A     = u1p_s_    + u1p_a_;
    const complex<double> u1p_B     = u1p_s_    - u1p_a_;
    const complex<double> u1p_A_out = u1p_s_out + u1p_a_out;
    const complex<double> u1p_B_out = u1p_s_out - u1p_a_out;
    
    const complex<double> u3_A     = u3_s_    + u3_a_;
    const complex<double> u3_B     = u3_s_    - u3_a_;
    const complex<double> u3_A_out = u3_s_out + u3_a_out;
    const complex<double> u3_B_out = u3_s_out - u3_a_out;
    
    double ans = 0.;
    
    ans += 2.*(V1_+2.*V2_) * (rho_A*rho_B_out + rho_B*rho_A_out); // Density terms
    ans +=   (2.*V1p_+V3_) * (rho_A*rho_A_out + rho_B*rho_B_out);
    
    ans += -4.*V1_ * real( u1_ * conj(u1_out) ); // Exchange terms
    ans += -4.*V2_ * real( u2A_  * conj(u2A_out)    +  u2B_  * conj(u2B_out)   );
    ans += -2.*V1p_* real( u1p_A * conj(u1p_A_out)  +  u1p_B * conj(u1p_B_out) );
    ans += -   V3_ * real( u3_A  * conj(u3_A_out)   +  u3_B  * conj(u3_B_out)  );
    
    return ans;
}
*/