// chempot_test.cc
/* Testing suite for the module chempot. */

#include <iostream> // For terminal input and output
#include <iomanip> // For the function std::setprecision()
#include <complex>
#include <cmath>
#include "chempot.h" // Module's header

void test_ChemPot()
{
    /* Test for the routine ChemPotBisec() */
    
    const int ArrLen = 100*100*100;
    double energies1 [ArrLen] = {0.};
    for (int i=0; i<ArrLen; ++i)
      energies1[i] = std::sqrt((double)(i));
    
    const int num_electrons = (int)((double)(ArrLen)*0.639);
    std::cout << "num_electrons = " << num_electrons << std::endl;
    const double T = 0.1;
    const double tol = 1.e-13;
    
    double mu = ChemPotBisec(ArrLen, num_electrons, energies1, T, tol, true);
    //double mu = ChemPotNewton(ArrLen, num_electrons, energies1, T, tol, 20, true);
    
    std::cout << std::setprecision(40) << "mu = " << mu << std::endl;
}



int main()
{
    test_ChemPot();
    
    return 0;
}
