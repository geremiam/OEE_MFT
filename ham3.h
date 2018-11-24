// ham3.h
/* Description */
#ifndef HAM3SOURCE_H
#define HAM3SOURCE_H

#include <complex>
using std::complex;

class ham3_t
{
  /* Class for holding the Hamiltonian parameters. Parameters than are not likely to be 
  changed in a parameter study are made private (and const). Parameters that may be 
  varied during a parameter study are public, except for those on which other parameters 
  depend; such parameters are modifiable with a public method. */
  
  private:
    
    // Private copy constructor (prohibits copy creation)
    ham3_t(const ham3_t&);
    // Private assignment operator (prohibits assignment)
    const ham3_t& operator=(const ham3_t&);
    
  //public:
    
    /* Parameters that the user doesn't need to modify after instantiation. */
    
    const double a_ = 1.; // a- and b-axis length
    const double c_ = 1.; // c-axis length
    
    /* Parameters for the momentum space grid. Assigned in the constructor. */
    const int ka_pts_; // Useful to allow the user to set these for convergence studies
    const int kb_pts_;
    const int kc_pts_;
    
    /* Size of the Hamiltonian, or equivalently number of bands */
    const int num_bands = 2; // The number of bands, i.e. the order of the matrix for each k
    const int ham_array_rows = num_bands; // Same as matrix order for full storage
    const int ham_array_cols = num_bands; // Same as matrix order for full storage
    const int num_states = ka_pts_*kb_pts_*kc_pts_*num_bands;
    
    const int loops_lim_ = 1000; // Limit to the number of iteration loops
    
    
    /* Parameters that have dependencies or are dependent on other parameters. These 
    should be modified with an appropriate method. */
    double rho_ = 0.4; // Average (global) electron density, between 0 and 1
    int filled_states=0; // Value assigned in constructor or upon calling 'assign_rho()'
    
    bool zerotemp_ = true; // Sets whether the temperature is zero or not
    double T_ = 0.; // Sets the temperature. Irrelevant when zerotemp_ is true.
    
    
    // Functions useful for computing the Hamiltonian
    double              zeta(double ka, double kb);
    complex<double> chi(double ka, double kb);
    complex<double>   f(double ka, double kb, double kc);
    double            ftilde(double ka, double kb, double kc);
    
    
    // These functions are used in the method ComputeMFs()
    double ComputeTerm_rho_a(const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u1   (const double ka, const double kb, const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u1p_s(const double ka, const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u1p_a(const double ka, const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u2A  (const double ka, const double kb, const double kc, const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u2B  (const double ka, const double kb, const double kc, const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u3_s (const double kc, const double*const occs, const complex<double>*const*const evecs);
    complex<double> ComputeTerm_u3_a (const double kc, const double*const occs, const complex<double>*const*const evecs);
    
    void Assign_ham(const double ka, const double kb, const double kc, complex<double>*const*const ham_array);
    
    void ComputeMFs(double& rho_a_out, complex<double>& u1_out, complex<double>& u1p_s_out, complex<double>& u1p_a_out,
                    complex<double>& u2A_out, complex<double>& u2B_out,complex<double>& u3_s_out, complex<double>& u3_a_out);
    
  public:
    
    /* Settings for the iterative search */
    
    // The MF values that are used in the Hamiltonian.
    // This is where the user sets the initial values of the MFs.
    
            double  rho_a_ = -99.;
    complex<double> u1_    = {-99.,0.};
    complex<double> u1p_s_ = {-99.,0.};
    complex<double> u1p_a_ = {-99.,0.};
    complex<double> u2A_   = {-99.,0.};
    complex<double> u2B_   = {-99.,0.};
    complex<double> u3_s_  = {-99.,0.};
    complex<double> u3_a_  = {-99.,0.};
    
    const double tol_ = 1.e-9; // Tolerance for the equality of the mean fields
    
    
    // Hamiltonian parameters that the user may want to change
    double t1_  = 1.; // x- and y-direction hopping
    double t1p_ = 0.7; // a- and b-direction hopping
    double t2A_ = 0.1; // "A-emanating" solenoid-like hopping
    double t2B_ = 0.4; // "B-emanating" solenoid-like hopping
    double t3_  = 0.5; // c-direction hopping
    double V1_  = 1.0; // Repulsion between neighbours in x and y directions
    double V1p_ = 0.7; // Repulsion between neighbours in a and b directions
    double V2_  = 0.4; // We choose to use the same V2 for A- and B-emanating
    double V3_  = 0.5; // Repulsion between neighbours in c direction
    
    
    void resetMFs(); // Resets MFs to default starting values.
    
    void assign_rho(const double rho); // Assign rho_ and dependent variables
    void set_zerotemp();
    void set_nonzerotemp(const double T);
    
    ham3_t(const int ka_pts=100, const int kb_pts=100, const int kc_pts=100); // Constructor declaration
    ~ham3_t(); // Destructor declaration
    
    bool FixedPoint(int*const num_loops_p=NULL, const bool with_output=false);
    
    std::string GetAttributes();
};

#endif
