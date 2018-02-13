// ham3_test.cc
#include <iostream>
#include <complex> // For complex numbers
#include "IO.h" // For printing matrices
#include "ham3.h"

void pars_t_test()
{
    pars_t pars(1.5);
    std::cout << "pars.t1 = " << pars.t1 << std::endl;
    std::cout << "pars.t2 = " << pars.t2 << std::endl;
    std::cout << "pars.eps = " << pars.eps << std::endl;
    std::cout << "pars.phi = " << pars.phi << std::endl;
    std::cout << std::endl;
    
    pars.SetScaling(1.);
    std::cout << "pars.t1 = " << pars.t1 << std::endl;
    std::cout << "pars.t2 = " << pars.t2 << std::endl;
    std::cout << "pars.eps = " << pars.eps << std::endl;
    std::cout << "pars.phi = " << pars.phi << std::endl;
    std::cout << std::endl;
    
    pars.SetScaling(0.75);
    std::cout << "pars.t1 = " << pars.t1 << std::endl;
    std::cout << "pars.t2 = " << pars.t2 << std::endl;
    std::cout << "pars.eps = " << pars.eps << std::endl;
    std::cout << "pars.phi = " << pars.phi << std::endl;
}

void Assign_ham_test()
{
    ham3_t ham3;
    
    ham3.parsI.SetScaling(0.75);
    ham3.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham3.ham_array,std::cout);
    
    std::cout << std::endl;
    
    ham3.U = 0.5;
    ham3.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham3.ham_array,std::cout);
}



int main()
{
    //pars_t_test();
    Assign_ham_test();
    return 0;
}
