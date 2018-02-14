// ham3.cc
/* Description */

#include <iostream>
#include <complex> // Needed to import alloc_dealloc.h
#include "alloc_dealloc.h" // Allocation/deallocation of arrays
#include "init_routines.h" // Initialization of arrays
#include "math_routines.h" // Various math functions
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



pars_t::pars_t(const double lambda)
{
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
    :ham_array(Alloc2D_z(ham_array_rows, ham_array_cols))
{
    /* Constructor implementation. */
    ValInitArray(ham_array_rows*ham_array_cols, &(ham_array[0][0]));//Initialize to zero
    std::cout << "ham3_t instance created.\n";
}

ham3_t::~ham3_t()
{
    // Destructor implementation
    Dealloc2D(ham_array);
    std::cout << "ham3_t instance deleted.\n";
}

void ham3_t::Assign_ham(const double kx, const double ky)
{
    /* Given the parameters kx, ky, and M_, calculate the 8*8 k-space Hamiltonian and 
    assign it to H in full storage layout. The diag routine only uses the lower triangle, 
    so we only assign that part. */
    
    // We calculate the kinetic energy parts of the Hamiltonian
    std::complex<double> hI [4] = {0.,0.};
    Assign_h(kx, ky, a, parsI, hI);
    std::complex<double> hII [4] = {0.,0.};
    Assign_h(kx, ky, a, parsII, hII);
    
    std::complex<double>*const*const& H = ham_array;//For convenience, define reference
    
    H[0][0]=hI[0]+U*rhoI/2.+L/2.; 
    H[1][0]=hI[2]; H[1][1]=hI[3]+U*rhoI/2.+L/2.; 
    
    H[2][0]=-U*mag; H[2][1]=0.;    H[2][2]=hI[0]+U*rhoI/2.+L/2.; 
    H[3][0]=0.; H[3][1]=+U*mag;    H[3][2]=hI[2]; H[3][3]=hI[3]+U*rhoI/2.+L/2.; 
    
    H[4][0]=tperp; H[4][1]=0.;    H[4][2]=0.; H[4][3]=0.;    H[4][4]=hII[0]-L/2.; 
    H[5][0]=0.; H[5][1]=tperp;    H[5][2]=0.; H[5][3]=0.;    H[5][4]=hII[2]; H[5][5]=hII[3]-L/2.; 
    
    H[6][0]=0.; H[6][1]=0.;    H[6][2]=tperp; H[6][3]=0.;    H[6][4]=0.; H[6][5]=0.;    H[6][6]=hII[0]-L/2.; 
    H[7][0]=0.; H[7][1]=0.;    H[7][2]=0.; H[7][3]=tperp;    H[7][4]=0.; H[7][5]=0.;    H[7][6]=hII[2]; H[7][7]=hII[3]-L/2.;
    
}

double ham3_t::ComputeTerm_mag(const double mu, const double*const evals, 
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
        std::cerr << "WARNING: M has nonzero imaginary part: " << imag_part << std::endl;
    
    return std::real(accumulator/(double)(4*kx_pts*ky_pts));
}

double ham3_t::ComputeTerm_rhoI(const double mu, const double*const evals, 
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
        std::cerr << "WARNING: rhoI has nonzero imaginary part: " << imag_part << std::endl;
    
    return std::real(accumulator/(double)(2*kx_pts*ky_pts));
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
        "; U = "+to_string(U)+"; mag_startval = "+to_string(mag_startval)+
        "; rhoI_startval = "+to_string(rhoI_startval)+"; tol = "+to_string(tol)+
        "; loops_lim = "+to_string(loops_lim);
    
    return Attributes;
}
