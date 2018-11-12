// math_routines_test.cc
/* Testing suite for the module math_routines. */

#include <iostream> // For terminal input and output
#include <complex>
#include "math_routines.h" // Module's header

void test_nF0()
{
    /* Test for the routine nF0 */
    
    const int ArrLen = 4;
    const double E [ArrLen] = {-1., -0.1, 0.2, 0.7};
    const float  e [ArrLen] = {-1., -0.1, 0.2, 0.7};
    
    const std::complex<double> z = {2., 3.};
    
    for (int i=0; i<ArrLen; ++i)
        std::cout << "nF0(" << E[i] << ") = " << z*nF0(E[i]) << "\t"
                  << "nF0(" << e[i] << ") = " << nF0(e[i]) << std::endl;
}

void test_nF()
{
    /* Test for the routine nF */
    
    const int ArrLen = 4;
    const double E [ArrLen] = {-1., -0.1, 0.2, 0.7};
    const float  e [ArrLen] = {-1., -0.1, 0.2, 0.7};
    
    const double BETA = 1./0.1;
    const float  beta = 1./0.1;
    
    for (int i=0; i<ArrLen; ++i)
        std::cout << "nF(BETA, " << E[i] << ") = " << nF(BETA, E[i]) << "\t"
                  << "nF(beta, " << e[i] << ") = " << nF(beta, e[i]) << std::endl;
    
}

void test_FermiEnerg_cpy()
{
    const int num_states = 7;
    const double filled_states = num_states/2;
    double const energies [num_states] = {0.9, 0.7, 0.8, 0.4, 0.5, 0.3, 0.45};
    
    const double FE = FermiEnerg_cpy(num_states, (int)(filled_states), energies);
    
    std::cout << "Fermi energy is " << FE << std::endl;
}

void test_FermiEnerg()
{
    const int num_states = 7;
    const double filled_states = num_states/2;
    double energies [num_states] = {0.9, 0.7, 0.8, 0.4, 0.5, 0.3, 0.45};
    
    const double FE = FermiEnerg(num_states, (int)(filled_states), energies);
    
    std::cout << "Fermi energy is " << FE << std::endl;
}

void test_Occupations()
{
    const int arrlen = 50;
    
    double energies [arrlen] = {0.};
    for (int i=0; i<arrlen; ++i)
    {
        energies[i] = 0.1*(double)(i);
        std::cout << "E[" << i << "]=" << energies[i] << "\n";
    }
    std::cout << std::endl;
    
    const double mu = 3.55;
    const bool zerotemp = true;
    const double T = 0.001;
    
    double occs [arrlen] = {0.};
    Occupations(arrlen, mu, energies, occs, zerotemp, T);
    
    for (int i=0; i<arrlen; ++i)
    {
        std::cout << "occs[" << i << "]=" << occs[i] << "\n";
    }
}

int main()
{
    //test_nF0();
    
    //test_nF();
    
    //test_FermiEnerg_cpy();
    
    //test_FermiEnerg();
    
    //test_Occupations();
    
    return 0;
}
