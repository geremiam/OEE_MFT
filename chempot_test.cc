// chempot_test.cc
/* Testing suite for the module chempot. */

#include <iostream> // For terminal input and output
#include <iomanip> // For the function std::setprecision()
#include <complex>
#include <cmath>
#include <limits> // To get machine precision
#include "chempot.h" // Module's header

void test_ChemPot()
{
    /* Test for the routine ChemPotBisec() */
    
    // Get machine precision
    const double eps_double = std::numeric_limits<double>::epsilon();
    std::cout << "eps_double = " << eps_double << std::endl;
    // Seems the tolerance can safely be made about as small as epsilon times the root.
    
    const int ArrLen = 100*100*100;
    double energies1 [ArrLen] = {0.};
    for (int i=0; i<ArrLen; ++i)
      energies1[i] = std::sqrt((double)(i)) * 0.0012593;
    
    const int num_electrons = (int)((double)(ArrLen)*0.639);
    std::cout << "num_electrons = " << num_electrons << std::endl;
    const double T = 0.1;
    
    const bool show_output = true;
    const bool usethreads = true;  // Bisection method can use OpenMP parallelization
    double mu = ChemPotBisec(ArrLen, num_electrons, energies1, T, show_output, usethreads);
    //double mu = ChemPotNewton(ArrLen, num_electrons, energies1, T, tol, 20, show_output);
    
    std::cout << std::setprecision(20) << "mu = " << mu << std::endl;
}



int main()
{
    test_ChemPot();
    
    return 0;
}
