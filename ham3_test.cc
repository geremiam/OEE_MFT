// ham3_test.cc
#include <iostream>
#include <complex> // For complex numbers
#include <cmath> // For the constant M_PI
#include "IO.h" // For printing matrices
#include "alloc.h"
#include "ham3.h"

void Assign_ham_test()
{
    std::complex<double>** matrix = Alloc2D_z(2, 2);
    std::cout << "\nmatrix = " << std::endl;
    PrintMatrix(2,2,matrix, std::cout);
    std::cout << std::endl;
    
    const int ka_pts = 10;
    const int kb_pts = 10;
    const int kc_pts = 10;
    
    ham3_t ham3(ka_pts, kb_pts, kc_pts);
    
    std::cout << "num_states = " << ham3.num_states << std::endl;
    std::cout << "rho_ = " << ham3.rho_ << std::endl;
    std::cout << "filled_states = " << ham3.filled_states << std::endl;
    
    const double ka=M_PI, kb=M_PI, kc=M_PI;
    ham3.Assign_ham(ka, kb, kc, matrix);
    std::cout << "\nAfter assignment, matrix = " << std::endl;
    PrintMatrix(2,2,matrix, std::cout);
    std::cout << std::endl;
    
    /*
    ham3.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham3.ham_array,std::cout);
    
    std::cout << std::endl;
    
    ham3.U = 0.5;
    ham3.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham3.ham_array,std::cout);
    */
    Dealloc2D(matrix);
}

void test_FixedPoint()
{
    const int ka_pts = 100;
    const int kb_pts = 100;
    const int kc_pts = 100;
    
    ham3_t ham3(ka_pts, kb_pts, kc_pts);
    
    int num_loops;
    const bool with_output = true;
    
    ham3.FixedPoint(&num_loops, with_output);
    
    std::cout << "num_loops = " << num_loops << std::endl;
}



int main()
{
    //pars_t_test();
    //Assign_ham_test();
    
    test_FixedPoint();
    
    return 0;
}
