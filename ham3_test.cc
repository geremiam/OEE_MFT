// ham3_test.cc
#include <iostream>
#include <iomanip> // For the function std::setprecision()
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
    ham3.set_zerotemp(); // In this case, calculates EF from sorting
    //ham3.set_nonzerotemp(1.e-2); // In this case, calculates mu from bisection method
    
    int num_loops;
    const bool with_output = true;
    
    ham3.FixedPoint(&num_loops, with_output);
    
    std::cout << "num_loops = " << num_loops << std::endl;
}


void test_ComputeMFs()
{
    // Paste into the routine ComputeMFs() to inspect behaviour:
    //std::cout<<std::fixed<<std::setprecision(4)<<kspace.ka_grid[i]<<"  "<<kspace.kb_grid[j]<<"  "<<kspace.kc_grid[k]<<std::setprecision(10);
    //std::cout<<"  evals = "<<evals[0]<<" "<<evals[1]<<"\tmu = "<<mu<<"\toccs = "<<occs[0]<<" "<<occs[1]<<std::scientific;
    //std::cout << "  " << ComputeTerm_u1(kspace.ka_grid[i], kspace.kb_grid[j], occs, evecs) << std::endl;
    std::cout << std::scientific << std::showpos << std::setprecision(10);
    
    const int ka_pts = 4;
    const int kb_pts = 4;
    const int kc_pts = 4;
    
    ham3_t ham3(ka_pts, kb_pts, kc_pts);
    ham3.set_nonzerotemp(1.e-1); // In this case, calculates EF from sorting
    //ham3.assign_rho(0.6);
    //ham3.V1_ = 2.53e+01;
    //ham3.V1p_ = 1.89e+01;
    //ham3.V2_ = 3.57e-01;
    //ham3.V3_ = +1.79e-01;
    
    ham3.rho_a_ = 1.12051079202081405040e-01;
    ham3.u1_    = {+1.81809396905388331867e-01,+2.49800180540660221595e-16};
    ham3.u1p_s_ = {+6.62587689134149226966e-02,-1.83686399424232149613e-13};
    ham3.u1p_a_ = {+3.38161416511557379183e-03,-3.81639164714897560771e-16};
    ham3.u2A_   = {-1.72032695934012970496e-02,-1.28369537222283724986e-16};
    ham3.u2B_   = {+3.26383205594113143255e-02,+1.95156391047390798121e-16};
    ham3.u3_s_  = {+5.19325044381043557373e-02,-7.85728937696770352028e-18};
    ham3.u3_a_  = {+1.07683841914983473298e-03,-2.68456504720409813713e-19};
    
    
    double rho_a_out;
    complex<double> u1_out;
    complex<double> u1p_s_out;
    complex<double> u1p_a_out;
    complex<double> u2A_out;
    complex<double> u2B_out;
    complex<double> u3_s_out;
    complex<double> u3_a_out;
    
    ham3.ComputeMFs(rho_a_out, u1_out, u1p_s_out, u1p_a_out, u2A_out, u2B_out, u3_s_out, u3_a_out);
    
    std::cout << "\t" << rho_a_out << "\t" << u1_out << "\t" 
              << u1p_s_out << "\t" << u1p_a_out << "\t" 
              << u2A_out   << "\t" << u2B_out   << "\t" 
              << u3_s_out  << "\t" << u3_a_out  << std::endl;
}


int main()
{
    //pars_t_test();
    //Assign_ham_test();
    
    //test_ComputeMFs();
    
    test_FixedPoint();
    
    return 0;
}
