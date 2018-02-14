// ham3.h
/* Description */
#ifndef HAM3SOURCE_H
#define HAM3SOURCE_H

const double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

// Class for holding the parameters proper to each layer
class pars_t
{
  private:
    // Default values of the parameters
    const double t1_0 = 1.; // NN hopping
    const double t2_0 = 0.25; // NNN hopping
    const double eps_0 = 0.; // Potential difference between A and B sublattices
  public:
    // Default values of the parameters
    const double phi = pi/2.; // Flux phase in the Haldane model
    double t1=t1_0; // NN hopping
    double t2=t2_0; // NNN hopping
    double eps=eps_0; // Potential difference between A and B sublattices
    
    pars_t(const double lambda=1.);
    void SetScaling(const double lambda);
};

class ham3_t
{
  private:
    
    // Private copy constructor (prohibits copy creation)
    ham3_t(const ham3_t&);
    // Private assignment operator (prohibits assignment)
    const ham3_t& operator=(const ham3_t&);
    
    const double a = 1.; // We take a to be the NN distance
    
  public:
    
    /* Settings for the iterative search */
    const double mag_startval = 0.1; // Choose a starting value
    const double rhoI_startval = 0.5; // Choose starting value
    const int loops_lim = 3000; // Limit to the number of iteration loops
    const double tol = 1.e-6; // Tolerance for the equality of the mean fields
    
    /* Size of the Hamiltonian, or equivalently number of bands (constant) */
    const int bands_num = 8; // The number of bands, i.e. the order of the matrix for each k
    const int ham_array_rows = bands_num; // Same as matrix order for full storage
    const int ham_array_cols = bands_num; // Same as matrix order for full storage
    
    /* Make sure the rectangular zone used here is equivalent to the first Brillouin zone. */
    const int kx_pts = 173;
    const double kx_bounds [2] = {-(2.*pi)/(3.*sqrt(3.)*a), (4.*pi)/(3.*sqrt(3.)*a)};
    const int ky_pts = 200;
    const double ky_bounds [2] = {-(2.*pi)/(3.*a), (2.*pi)/(3.*a)};
    
    
    
    const double tperp_0 = 0.3; // Base value of tperp (gets scaled)
    double tperp = tperp_0;
    double L = 0.; // bias voltage between layers I and II
    double rho = 1.; // Average (global) electron density, between 0 and 2
    double U = 0.; // Hubbard interaction strength
    
    pars_t parsI, parsII;
    
    double mag=-88.;
    double rhoI=-88.;
    
    const int num_states = kx_pts*ky_pts*bands_num;
    int filled_states = (int)( rho * (double)(4*kx_pts*ky_pts) );
    
    
    std::complex<double>*const*const ham_array; // array for storing the Hamiltonian
    
    
    ham3_t(); // Constructor declaration
    ~ham3_t(); // Destructor declaration
    
    void Assign_ham(const double kx, const double ky);
    
    double ComputeTerm_mag(const double mu, const double*const evals, 
                           const std::complex<double>*const*const evecs);
    double ComputeTerm_rhoI(const double mu, const double*const evals, 
                            const std::complex<double>*const*const evecs);
    
    std::string GetAttributes();
};

#endif
