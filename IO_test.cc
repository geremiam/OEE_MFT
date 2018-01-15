// IO_test.cc
/* Testing suite for the module misc. */
#include <iostream>
#include <complex>
#include "alloc_dealloc.h" //Allocation and deallocation of arrays (must include complex)
#include "IO.h" // Include header file for consistency check

int test_PrintMatrix()
{
    /* Unit test for the (overloaded) function PrintMatrix() */
    const int num_rows = 2;
    const int num_cols = 2;
    
    double*const*const A = Alloc2D_d(num_rows, num_cols);
    for (int i=0; i<num_rows*num_cols; ++i) A[0][i] = (double)(i);
    PrintMatrix(num_rows, num_rows, A, std::cout);
    Dealloc2D(A);
    
    std::cout << std::endl;
    
    double*const B = new double [num_rows*num_cols];
    for (int i=0; i<num_rows*num_cols; ++i) B[i] = (double)(i*i);
    PrintMatrix(num_rows, num_rows, B, std::cout);
    delete [] B;
    
    std::cout << std::endl;
    
    std::complex<double>*const*const C = Alloc2D_z(num_rows, num_cols);
    for (int i=0; i<num_rows*num_cols; ++i) C[0][i] = {0.28*(double)(i), -1.23*(double)(i)};
    PrintMatrix(num_rows, num_rows, C, std::cout);
    Dealloc2D(C);
    
    std::cout << std::endl;
    
    std::complex<double>*const D = new std::complex<double> [num_rows*num_cols];
    for (int i=0; i<num_rows*num_cols; ++i) D[i] = {0.28*(double)(i), -1.23*(double)(i)};
    PrintMatrix(num_rows, num_rows, D, std::cout);
    delete [] D;
    
}

int main()
{
    test_PrintMatrix();
    
    return 0;
}
