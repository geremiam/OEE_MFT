// InitRoutines.cc

#include <iostream>
#include "init_routines.h"

int InitArray_TEST()
{
    const int NumPts = 20;
    double Array [NumPts] = {0.};
    const double Bounds [2] = {-3., 3.};
    
    std::cout << "dx = " << InitArray(Bounds[0], Bounds[1], NumPts, Array) 
              << std::endl;
    
    for (int i=0; i<NumPts; ++i)
        std::cout << "Array[" << i << "] = " << Array[i] << std::endl;
}


int main()
{
    InitArray_TEST();
    
    return 0;
}