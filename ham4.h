// ham4.h
/* Description */
#ifndef HAM4SOURCE_H
#define HAM4SOURCE_H

#include <complex>
#include <cmath> // For the constant M_PI
using std::complex;

class ham4_t
{
  /* Class for holding the Hamiltonian parameters. Parameters than are not likely to be 
  changed in a parameter study are made private (and const). Parameters that may be 
  varied during a parameter study are public, except for those on which other parameters 
  depend; such parameters are modifiable with a public method. */
  
  private:
    
    // Private copy constructor (prohibits copy creation)
    ham4_t(const ham4_t&);
    // Private assignment operator (prohibits assignment)
    const ham4_t& operator=(const ham4_t&);
    
  public:
    
    /* Parameters that the user doesn't need to modify after instantiation. */
    
    const double a_ = 1.; // a- and b-axis length
    const double c_ = 1.; // c-axis length
    
    // Parameters for the momentum space grid. Assigned in the constructor.
    // Sets the number of momentum points in the REDUCED BZ
    const int ka_pts_; // Useful to allow the user to set these for convergence studies
    const int kb_pts_;
    const int kc_pts_;
    const int states_per_cell = 2; // Number of states in an (original) unit cell
    static const int num_harmonics = 4; // Number of harmonics considered.
    
    const double Qa [num_harmonics] = {0., M_PI/a_, 0.,      M_PI/a_};
    const double Qb [num_harmonics] = {0., 0.,      M_PI/a_, M_PI/a_};
    const double Qc [num_harmonics] = {0., 0.,      0.,      0.};
    
    const int addition_table [num_harmonics][num_harmonics]={{0, 1, 2, 3},
                                                             {1, 0, 3, 2},
                                                             {2, 3, 0, 1},
                                                             {3, 2, 1, 0}};
    
    // Size of the Hamiltonian, or equivalently number of bands
    const int num_bands = states_per_cell * num_harmonics; // The number of bands, i.e. the order of the matrix for each k
    const int ham_array_rows = num_bands; // Same as matrix order for full storage
    const int ham_array_cols = num_bands; // Same as matrix order for full storage
    const int num_states = ka_pts_*kb_pts_*kc_pts_*num_bands;
    
    const int loops_lim_ = 1400; // Limit to the number of iteration loops
    
    
    /* Parameters that have dependencies or are dependent on other parameters. These 
    should be modified with an appropriate method. */
    double rho_ = 0.4; // Average (global) electron density, between 0 and 1
    int filled_states=0; // Value assigned in constructor or upon calling 'assign_rho()'
    
    bool zerotemp_ = true; // Sets whether the temperature is zero or not
    double T_ = 0.; // Sets the temperature. Irrelevant when zerotemp_ is true.
    
    // Function for performing the direct product in the indices Q and alpha.
    int idx(const int Q, const int alpha) const;
    
    // Functions useful for computing the Hamiltonian
    double          zeta(double ka, double kb) const;
    complex<double> f(double ka, double kb, double kc) const;
    
    
    // These functions are used in the method ComputeMFs()
    void AddContribution_rho_s(const double*const occs, const complex<double>*const*const evecs, double*const rho_s) const;
    void AddContribution_rho_a(const double*const occs, const complex<double>*const*const evecs, double*const rho_a) const;
    
    void Assign_h(double ka, double kb, double kc, complex<double>*const h) const;
    void Assign_V(const int Q, complex<double>*const V) const;
    void Assign_V_manual(const int Q, complex<double>*const V) const;
    void Assign_ham(const double ka, const double kb, const double kc, complex<double>*const*const ham_array) const;
    
    double ComputeMFs(double*const rho_s_out, double*const rho_a_out) const;
    
    
  //public:
    
    /* Settings for the iterative search */
    
    // The MF values that are used in the Hamiltonian.
    // This is where the user sets the initial values of the MFs.
    
    double rho_A_ [num_harmonics] = {1.,1.,1.,1.};//{0.};
    double rho_B_ [num_harmonics] = {0.};
    
    double HFE_ = -99.; // For storing the free energy
    
    const double tol_ = 1.e-6; // Tolerance for the equality of the mean fields
    
    
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
    
    
    //void resetMFs(); // Resets MFs to default starting values.
    
    void assign_rho(const double rho); // Assign rho_ and dependent variables
    void set_zerotemp();
    void set_nonzerotemp(const double T);
    
    ham4_t(const int ka_pts=62, const int kb_pts=62, const int kc_pts=62); // Constructor declaration
    ~ham4_t(); // Destructor declaration
    
    bool FixedPoint(int*const num_loops_p=NULL, const bool with_output=false);
    
    std::string GetAttributes();
};

#endif
