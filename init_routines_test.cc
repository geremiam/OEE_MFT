// InitRoutines.cc

#include <iostream>
#include "init_routines.h"

void LinInitArray_TEST()
{
    const int NumPts = 20;
    double Array [NumPts] = {0.};
    const double Bounds [2] = {-3., 3.};
    
    std::cout << "dx = " << LinInitArray(Bounds[0], Bounds[1], NumPts, Array) 
              << std::endl;
    
    for (int i=0; i<NumPts; ++i)
        std::cout << "Array[" << i << "] = " << Array[i] << std::endl;
}

void ValInitArray_TEST()
{
    const int NumPts = 20;
    double Array [NumPts];
    
    for (int i=0; i<NumPts; ++i)
        std::cout << "Array[" << i << "] = " << Array[i] << std::endl;
    std::cout << std::endl;
    
    ValInitArray(NumPts, Array, 1.);
    
    for (int i=0; i<NumPts; ++i)
        std::cout << "Array[" << i << "] = " << Array[i] << std::endl;
}


int main()
{
    //LinInitArray_TEST();
    ValInitArray_TEST();
    
    return 0;
}