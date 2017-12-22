// diag_routines_test.cc
// Testing suite for the module diag_routines.
#include <iostream> // Input/output to command line
#include <complex> // For complex numbers
#include "diag_routines.h"

void PrintMatrix(const int order, std::complex<double> const *const a);
void PrintVector(const int len, double const *const a);

void test_simple_zheev()
{
    const int order = 4;
    std::complex<double> a [order*order]= {{1,2},   {0,0},    {0,0},   {0,0},
                                           {3,4},   {5,6},    {0,0},   {0,0},
                                           {7,8},   {9,10},   {11,12}, {0,0},
                                           {13,14}, {15, 16}, {17,18}, {19,20}};
    std::complex<double> z [order*order];
    
    std::cout << "a =\n";
    PrintMatrix(order, a);
    
    double w [order]={0.};
    
    const bool OutputEvecs = true;
    simple_zheev(order, a, w, OutputEvecs, z);
    
    std::cout << "\nw = ";
    PrintVector(order, w);
    
    std::cout << "\na =\n";
    PrintMatrix(order, a);
    
    std::cout << "\nz =\n";
    PrintMatrix(order, z);
}

int main()
{
    test_simple_zheev();
    return 0;
}


void PrintMatrix(const int order, std::complex<double> const *const a)
{
    for (int i=0; i<order*order; ++i)
        {
            std::cout << a[i] << "\t";
            if (i%order==(order-1)) std::cout << std::endl;
        }
}
void PrintVector(const int len, double const *const w)
{
    for (int i=0; i<len; ++i)
        std::cout << w[i] << " ";
    std::cout << std::endl;
}
