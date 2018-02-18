// ham3.cc
/* Description */

#include <iostream>
#include <complex> // Needed to import alloc_dealloc.h
#include "alloc_dealloc.h" // Allocation/deallocation of arrays
#include "init_routines.h" // Initialization of arrays
#include "math_routines.h" // Various math functions
#include "kspace.h" // Defines a class for holding a band structure
#include "diag_routines.h" // Routines for finding evals and evecs
#include "ham3.h" // Include header file for consistency check
//using std::cos; using std::sin; using std::conj; // Not sure if these are necessary
using std::polar;
using std::to_string;

void Assign_h(const double kx, const double ky, const double a, const pars_t pars,
              std::complex<double>*const h)
{
    // Define local variables for convenience
    const double& t1 = pars.t1; const double& t2 = pars.t2;
    const double& eps = pars.eps; const double& phi = pars.phi;
    /* Assigns to h (a 1D array in row-major layout) the values of the momentum-dependant 
    2*2 Haldane Hamiltonian. */
    h[0] = +eps - 2.*t2*( 2.*cos(sqrt(3.)/2.*a*kx-phi)*cos(3./2.*a*ky) + cos(sqrt(3.)*a*kx+phi) );
    h[1] = -t1*( polar(1.,a*ky) + polar(1.,-a*ky/2.)*2.*cos(sqrt(3.)/2.*a*kx) );
    h[2] = conj(h[1]);
    h[3] = -eps - 2.*t2*( 2.*cos(sqrt(3.)/2.*a*kx+phi)*cos(3./2.*a*ky) + cos(sqrt(3.)*a*kx-phi) );
}

double del2(const double p0, const double p1, const double d)
{
    return p0 - (p1-p0)*(p1-p0)/d;
}



pars_t::pars_t(const double lambda)
{
    /* Upon construction, the scaling is set to lambda (default value is 1). */
    SetScaling(lambda);
}

void pars_t::SetScaling(const double lambda)
{
    // The energy parameters get scaled by lambda.
    t1 = t1_0 * lambda; // NN hopping
    t2 = t2_0 * lambda; // NNN hopping
    eps = eps_0 * lambda; // Potential difference between A and B sublattices
}



ham3_t::ham3_t()
    :ham_array_(Alloc2D_z(ham_array_rows, ham_array_cols)),
    evecs_(Alloc2D_z(bands_num, bands_num)), 
    kspace_(kx_pts, kx_bounds, ky_pts, ky_bounds, bands_num)
{
    /* Constructor implementation. ham_array is allocated and initialized. */
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array_[0][0]));//Initialize to zero
    ValInitArray(bands_num*bands_num, &(evecs_[0][0])); // Initialize to zero
    std::cout << "ham3_t instance created.\n";
}

ham3_t::~ham3_t()
{
    /* Destructor implementation. ham_array is deallocated. */
    Dealloc2D(ham_array_);
    Dealloc2D(evecs_);
    std::cout << "ham3_t instance deleted.\n";
}

void ham3_t::Assign_ham(const double kx, const double ky)
{
    /* Given the arguments kx and ky as well as member parameters, calculate the 8*8 
    k-space Hamiltonian and assign it to ham_array in full storage layout. The diag 
    routine only uses the lower triangle, so we only assign that part. */
    
    // We calculate the kinetic energy parts of the Hamiltonian
    std::complex<double> hI [4] = {0.,0.};
    Assign_h(kx, ky, a, parsI, hI);
    std::complex<double> hII [4] = {0.,0.};
    Assign_h(kx, ky, a, parsII, hII);
    
    std::complex<double>*const*const& H = ham_array_; //For convenience, define reference
    
    H[0][0]=hI[0]+U*(rhoI_s_+rhoI_a_)/2.+L/2.; 
    H[1][0]=hI[2]; H[1][1]=hI[3]+U*(rhoI_s_-rhoI_a_)/2.+L/2.; 
    
    H[2][0]=-U*(mag_s_+mag_a_); H[2][1]=0.;    H[2][2]=hI[0]+U*(rhoI_s_+rhoI_a_)/2.+L/2.; 
    H[3][0]=0.; H[3][1]=-U*(mag_s_-mag_a_);    H[3][2]=hI[2]; H[3][3]=hI[3]+U*(rhoI_s_-rhoI_a_)/2.+L/2.; 
    
    H[4][0]=tperp; H[4][1]=0.;    H[4][2]=0.; H[4][3]=0.;    H[4][4]=hII[0]-L/2.; 
    H[5][0]=0.; H[5][1]=tperp;    H[5][2]=0.; H[5][3]=0.;    H[5][4]=hII[2]; H[5][5]=hII[3]-L/2.; 
    
    H[6][0]=0.; H[6][1]=0.;    H[6][2]=tperp; H[6][3]=0.;    H[6][4]=0.; H[6][5]=0.;    H[6][6]=hII[0]-L/2.; 
    H[7][0]=0.; H[7][1]=0.;    H[7][2]=0.; H[7][3]=tperp;    H[7][4]=0.; H[7][5]=0.;    H[7][6]=hII[2]; H[7][7]=hII[3]-L/2.;
    
}

double ham3_t::ComputeTerm_rhoI_s(const double mu, const double*const evals, 
                                const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to rhoI from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double B [bands_num][bands_num] = {{1., 0., 0., 0.,    0., 0., 0., 0.}, 
                                             {0., 1., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 1., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 1.,    0., 0., 0., 0.},
                                             
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int b=0; b<bands_num; ++b)
        for (int c=0; c<bands_num; ++c)
            for (int d=0; d<bands_num; ++d)
                accumulator += conj(evecs[c][b])*B[c][d]*evecs[d][b]*nF0(evals[b]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(double)(2*kx_pts*ky_pts));
    if (std::abs(imag_part)>1.e-15)
        std::cerr << "WARNING: rhoI_s has nonzero imaginary part: " << imag_part << "\n";
    
    return std::real(accumulator/(double)(2*kx_pts*ky_pts));
}

double ham3_t::ComputeTerm_rhoI_a(const double mu, const double*const evals, 
                                const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to rhoI from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double B [bands_num][bands_num] = {{+1., 0., 0., 0.,    0., 0., 0., 0.}, 
                                             {0., -1., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., +1., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., -1.,    0., 0., 0., 0.},
                                             
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int b=0; b<bands_num; ++b)
        for (int c=0; c<bands_num; ++c)
            for (int d=0; d<bands_num; ++d)
                accumulator += conj(evecs[c][b])*B[c][d]*evecs[d][b]*nF0(evals[b]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(double)(2*kx_pts*ky_pts));
    if (std::abs(imag_part)>1.e-15)
        std::cerr << "WARNING: rhoI_a has nonzero imaginary part: " << imag_part << "\n";
    
    return std::real(accumulator/(double)(2*kx_pts*ky_pts));
}

double ham3_t::ComputeTerm_mag_s(const double mu, const double*const evals, 
                               const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to mag from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double C [bands_num][bands_num] = {{0., 0., 1., 0.,   0., 0., 0., 0.}, 
                                             {0., 0., 0., 1.,   0., 0., 0., 0.},
                                             {1., 0., 0., 0.,   0., 0., 0., 0.},
                                             {0., 1., 0., 0.,   0., 0., 0., 0.},
                                             
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int b=0; b<bands_num; ++b)
        for (int c=0; c<bands_num; ++c)
            for (int d=0; d<bands_num; ++d)
                accumulator += conj(evecs[c][b])*C[c][d]*evecs[d][b]*nF0(evals[b]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(double)(4*kx_pts*ky_pts));
    if (std::abs(imag_part)>1.e-15)
        std::cerr << "WARNING: mag_s has nonzero imaginary part: " << imag_part << "\n";
    
    return std::real(accumulator/(double)(4*kx_pts*ky_pts));
}

double ham3_t::ComputeTerm_mag_a(const double mu, const double*const evals, 
                               const std::complex<double>*const*const evecs)
{
    // Evaluates the contribution to mag from a single k (see notes)
    // Not good to implement matrix mult. by hand... but we will for simplicity
    const double C [bands_num][bands_num] = {{0., 0., +1., 0.,   0., 0., 0., 0.}, 
                                             {0., 0., 0., -1.,   0., 0., 0., 0.},
                                             {+1., 0., 0., 0.,   0., 0., 0., 0.},
                                             {0., -1., 0., 0.,   0., 0., 0., 0.},
                                             
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.},
                                             {0., 0., 0., 0.,    0., 0., 0., 0.}};
    // Calculate the trace (see notes)
    std::complex<double> accumulator = {0.,0.};
    for (int b=0; b<bands_num; ++b)
        for (int c=0; c<bands_num; ++c)
            for (int d=0; d<bands_num; ++d)
                accumulator += conj(evecs[c][b])*C[c][d]*evecs[d][b]*nF0(evals[b]-mu);
    
    // Test for zero imaginary part
    const double imag_part = std::imag(accumulator/(double)(4*kx_pts*ky_pts));
    if (std::abs(imag_part)>1.e-15)
        std::cerr << "WARNING: mag_a has nonzero imaginary part: " << imag_part << "\n";
    
    return std::real(accumulator/(double)(4*kx_pts*ky_pts));
}

void ham3_t::ComputeMFs(const double rhoI_s_in, const double rhoI_a_in, 
                        const double mag_s_in, const double mag_a_in, 
                        double& rhoI_s_out, double& rhoI_a_out, 
                        double& mag_s_out, double& mag_a_out)
{
    // Update (private) class members for mean-field values (used by Assign_ham())
    rhoI_s_ = rhoI_s_in;
    rhoI_a_ = rhoI_a_in; 
    mag_s_  = mag_s_in;
    mag_a_  = mag_a_in;
    
    // Given the parameters, diagonalize the Hamiltonian at each grid point
    for (int i=0; i<kx_pts; ++i)
      for (int j=0; j<ky_pts; ++j)
      {
        Assign_ham(kspace_.kx_grid[i], kspace_.ky_grid[j]);
        simple_zheev(bands_num, &(ham_array_[0][0]), &(kspace_.energies[i][j][0]));
      }
    
    // Use all energies to compute chemical potential (elements get reordered)
    double mu = FermiEnerg(num_states, filled_states, &(kspace_.energies[0][0][0]));
    
    // Use all the occupation numbers and the evecs to find the order parameter
    // Probably best to diagonalize a second time to avoid storing the evecs
    double rhoI_s_accumulator = 0.;
    double rhoI_a_accumulator = 0.;
    double mag_s_accumulator = 0.;
    double mag_a_accumulator = 0.;
    
    for (int i=0; i<kx_pts; ++i)
      for (int j=0; j<ky_pts; ++j)
      {
        Assign_ham(kspace_.kx_grid[i], kspace_.ky_grid[j]);
        simple_zheev(bands_num, &(ham_array_[0][0]), 
                                    &(kspace_.energies[i][j][0]), true, &(evecs_[0][0]));
        rhoI_s_accumulator += ComputeTerm_rhoI_s(mu,&(kspace_.energies[i][j][0]),evecs_);
        rhoI_a_accumulator += ComputeTerm_rhoI_a(mu,&(kspace_.energies[i][j][0]),evecs_);
        mag_s_accumulator  += ComputeTerm_mag_s(mu,&(kspace_.energies[i][j][0]),evecs_);
        mag_a_accumulator  += ComputeTerm_mag_a(mu,&(kspace_.energies[i][j][0]),evecs_);
      }
    rhoI_s_out = rhoI_s_accumulator;
    rhoI_a_out = rhoI_a_accumulator;
    mag_s_out  = mag_s_accumulator;
    mag_a_out  = mag_a_accumulator;
}

std::string ham3_t::GetAttributes()
{
    // Define a string of metadata
    const std::string Attributes = "Haldane Hubbard bilayer (ham3)"
        ": NN distance a = "+to_string(a)+"; kx_pts = "+to_string(kx_pts)+
        "; kx_bounds = "+to_string(kx_bounds[0])+", "+to_string(kx_bounds[1])+
        "; ky_pts = " + to_string(ky_pts)+
        "; ky_bounds = "+to_string(ky_bounds[0])+", "+to_string(ky_bounds[1])+
        "; bands_num = "+to_string(bands_num)+
        "; parsI = ("+to_string(parsI.t1) +","+to_string(parsI.t2)+","+
                      to_string(parsI.eps)+","+to_string(parsI.phi)+")"
        "; parsII = ("+to_string(parsII.t1) +","+to_string(parsII.t2)+","+
                       to_string(parsII.eps)+","+to_string(parsII.phi)+")"
        "; tperp = "+to_string(tperp)+"; L = "+to_string(L)+"; rho = "+to_string(rho)+
        "; U = "+to_string(U)+"; rhoI_s_startval = "+to_string(rhoI_s_startval)+
        "; rhoI_a_startval = "+to_string(rhoI_a_startval)+
        "; mag_s_startval = "+to_string(mag_s_startval)+
        "; mag_a_startval = "+to_string(mag_a_startval)+"; tol = "+to_string(tol)+
        "; loops_lim = "+to_string(loops_lim);
    
    return Attributes;
}

bool FixedPoint(double& rhoI_s, double& rhoI_a, double& mag_s, double& mag_a, 
                ham3_t& ham3, int*const num_loops_p, const bool with_output)
{
    /* Performs the iterative self-consistent search using the parameters from ham3. 
    The initial values of mag and rhoI are used as the starting values for the search; 
    the end values are also output to mag and rhoI. */
    std::cout << std::scientific << std::showpos; // Format display output
    
    /* Declare output variables and *initialize them to input values*. */
    double rhoI_s_out = rhoI_s;
    double rhoI_a_out = rhoI_a;
    double mag_s_out  = mag_s;
    double mag_a_out  = mag_a;
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        rhoI_s = rhoI_s_out;
        rhoI_a = rhoI_a_out; 
        mag_s  = mag_s_out;
        mag_a  = mag_a_out; // Update mean-field values
        if (with_output) std::cout << "rhoIs="  << rhoI_s
                                   << " rhoIa=" << rhoI_a
                                   << " ms="    << mag_s 
                                   << " ma="    << mag_a << "\t";
        
        // Evaluate function. Output is assigned to the 'out' arguments.
        ham3.ComputeMFs(rhoI_s, rhoI_a, mag_s, mag_a, 
                        rhoI_s_out, rhoI_a_out, mag_s_out, mag_a_out);
        
        // Print out final OP values
        if (with_output) std::cout << "drhoIs=" << rhoI_s_out - rhoI_s 
                                   << " drhoIa=" << rhoI_a_out - rhoI_a
                                   << " dms="  << mag_s_out - mag_s
                                   << " dma="  << mag_a_out - mag_a
                                   << std::endl;
        
        // Test for convergence
        converged =    (std::abs(rhoI_s_out - rhoI_s)<ham3.tol)
                    && (std::abs(rhoI_a_out - rhoI_a)<ham3.tol)
                    && (std::abs(mag_s_out - mag_s)<ham3.tol) 
                    && (std::abs(mag_a_out - mag_a)<ham3.tol);
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

bool Steffensen(double& rhoI_s, double& rhoI_a, double& mag_s, double& mag_a, 
                ham3_t& ham3, int*const num_loops_p, const bool with_output)
{
    /* Performs the iterative Steffensen search using the parameters from ham3. The 
    initial values of mag and rhoI are used as the starting values for the search; the 
    end values are also output to mag and rhoI. The algorithm does not seem to work... */
    /*https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fixed_point.html
    (see source code) */
    std::cout << std::scientific << std::showpos; // Format display output
    
    const int num_vars = 4;
    
    // Declare arrays for intermediate values
    double p  [num_vars] = {rhoI_s, rhoI_a, mag_s, mag_a};
    double p0 [num_vars] = {rhoI_s, rhoI_a, mag_s, mag_a};
    
    int counter = 0; // Define counter for number of loops
    bool converged=false, fail=false; // Used to stop the while looping
    do // Iterate until self-consistency is achieved
    {
        ++counter; // Increment counter
        
        for (int i=0; i< num_vars; ++i) p0[i] = p[i];
        
        if (with_output) std::cout << "rhoIs="  << p0[0]
                                   << " rhoIa=" << p0[1]
                                   << " ms="    << p0[2] 
                                   << " ma="    << p0[3] << "\t";
        
        double p1   [num_vars] = {0.};
        double p2   [num_vars] = {0.};
        double d    [num_vars] = {0.};
        double diff [num_vars] = {0.};
        
        ham3.ComputeMFs(p0[0], p0[1], p0[2], p0[3], 
                        p1[0], p1[1], p1[2], p1[3]);
        
        ham3.ComputeMFs(p1[0], p1[1], p1[2], p1[3], 
                        p2[0], p2[1], p2[2], p2[3]);
        
        for (int i=0; i<num_vars; ++i)
        {
            d[i] = p2[i] - 2.*p1[i] + p0[i];
            if (d[i]!=0.) p[i] = del2(p0[i], p1[i], d[i]);
            else p[i] = p2[i];
            diff[i] = std::abs(p[i] - p0[i]);
        }
        
        if (with_output) std::cout << "drhoIs="  << p[0] - p0[0]
                                   << " drhoIa=" << p[1] - p0[1]
                                   << " dms="    << p[2] - p0[2]
                                   << " dma="    << p[3] - p0[3]
                                   << std::endl;
        
        // Test for convergence
        converged = ( MaxArrayValue(diff, num_vars) < ham3.tol );
        fail = (!converged) && (counter>ham3.loops_lim);//Must come after converged line
        
    } while (!converged && !fail);
    
    // Assign the final values to the function arguments
    rhoI_s = p[0];
    rhoI_a = p[1];
    mag_s  = p[2];
    mag_a  = p[3];
    
    // Unless num_loops_p is the null pointer, assign the number of loops to its location
    if (num_loops_p!=NULL) *num_loops_p = counter;
    
    // We make sure that either converged or fail is true.
    if ((converged==true) && (fail==true) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (A)\n";
    if ((converged==false) && (fail==false) )
        std::cout << "OUPS 1: This option shouldn't have occurred! (B)\n";
    
    return fail;
}
