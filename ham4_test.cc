// ham4_test.cc
#include <iostream>
#include <iomanip> // For the function std::setprecision()
#include <complex> // For complex numbers
#include <cmath> // For the constant M_PI
#include "IO.h" // For printing matrices
#include "alloc.h"
#include "ham4.h"

void Assign_ham_test()
{
    std::complex<double>** matrix = Alloc2D_z(8, 8);
    std::cout << "\nmatrix = " << std::endl;
    PrintMatrix(8,8,matrix, std::cout);
    std::cout << std::endl;
    
    const int ka_pts = 10;
    const int kb_pts = 10;
    const int kc_pts = 10;
    
    ham4_t ham4(ka_pts, kb_pts, kc_pts);
    
    std::cout << "num_states = " << ham4.num_states << std::endl;
    std::cout << "rho_ = " << ham4.rho_ << std::endl;
    std::cout << "filled_states = " << ham4.filled_states << std::endl;
    
    const double ka=0.1*M_PI, kb=0.1*M_PI, kc=0.1*M_PI;
    ham4.Assign_ham(ka, kb, kc, matrix);
    std::cout << "\nAfter assignment, matrix = " << std::endl;
    PrintMatrix(8,8,matrix, std::cout);
    std::cout << std::endl;
    
    /*
    ham4.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham4.ham_array,std::cout);
    
    std::cout << std::endl;
    
    ham4.U = 0.5;
    ham4.Assign_ham(0.,0.);
    PrintMatrix(8,8,ham4.ham_array,std::cout);
    */
    Dealloc2D(matrix);
}

void test_ComputeMFs()
{
    const int ka_pts = 62;
    const int kb_pts = 62;
    const int kc_pts = 62;
    
    ham4_t ham4(ka_pts, kb_pts, kc_pts);
    ham4.set_nonzerotemp(1.e-1);
    
    std::cout << "num_states: " << ham4.num_states << std::endl;
    std::cout << "filled_states: " << ham4.filled_states << std::endl;
    
    double rho_s_out [4] = {0.};
    double rho_a_out [4] = {0.};
    
    ham4.ComputeMFs(rho_s_out, rho_a_out);
    
    for (int i=0; i<4; ++i)
    {
      std::cout << "rho_s_out[" << i << "] = " << rho_s_out[i] << "\t"
                << "rho_a_out[" << i << "] = " << rho_a_out[i] << std::endl;
    }
}

void test_FixedPoint()
{
    const int ka_pts = 62;
    const int kb_pts = 62;
    const int kc_pts = 62;
    
    ham4_t ham4(ka_pts, kb_pts, kc_pts);
    ham4.set_zerotemp(); // In this case, calculates EF from sorting
    //ham4.set_nonzerotemp(1.e-2); // In this case, calculates mu from bisection method
    
    int num_loops;
    const bool with_output = true;
    
    ham4.FixedPoint(&num_loops, with_output);
    
    std::cout << "num_loops = " << num_loops << std::endl;
}

int main()
{
    //Assign_ham_test();
    
    //test_ComputeMFs();
    
    test_FixedPoint();
    
    return 0;
}
