// misc_routines_test.cc
/* Testing suite for the module math_routines. */

#include <iostream> // For terminal input and output
#include "misc_routines.h" // Module's header

void test_check_bound_array()
{
    const double maxval = 0.31;
    const int len = 4;
    const double array1 [len] = {0.0,  0.1,  0.2,  0.15};
    const double array2 [len] = {-1.,  0.1,  0.2,  0.3};
    const double array3 [len] = {0.0,  0.145,0.2,  0.41};
    const double array4 [len] = {0.0,  0.1, -0.2,  0.3};
    const double array5 [len] = {0.0,  0.1, -0.2, -0.3};
    const double array6 [len] = {0.5,  0.4,  -1.,  0.5};
    
    std::cout << check_bound_array(len, maxval, array1) << std::endl;
    std::cout << check_bound_array(len, maxval, array2) << std::endl;
    std::cout << check_bound_array(len, maxval, array3) << std::endl;
    std::cout << check_bound_array(len, maxval, array4) << std::endl;
    std::cout << check_bound_array(len, maxval, array5) << std::endl;
    std::cout << check_bound_array(len, maxval, array6) << std::endl;
}

int main()
{
    test_check_bound_array();
    
    return 0;
}
