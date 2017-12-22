// nc_IO_test.cc

#include <iostream>
#include <string>
#include "nc_IO.h"



int main()
{
    const double Dim1 [4] = {1., 2., 3., 4.};
    const double Dim2 [5] = {10., 11., 12., 13., 14.};
    
    const std::string GlobalAttr = "The global attributes";
    
    const int dims_num_ = 2;
    const std::string dim_names [dims_num_] = {"Dim_1", "Dim_2"};
    const int dim_lengths [dims_num_] = {4, 6};
    
    const int vars_num_ = 2;
    const std::string var_names [vars_num_] = {"Var_1", "Var_2"};
    
    
    newDS_t newDS(GlobalAttr, dims_num_, dim_names, dim_lengths, vars_num_, var_names);
    
    const double*const coord_vars [2] = {Dim1, Dim2};
    newDS.WriteCoordVars(coord_vars);
    
    const double Var1 [4*6] = {11, 12, 13, 14, 15, 16,
                               21, 22, 23, 24, 25, 26,
                               31, 32, 33, 34, 35, 36,
                               41, 42, 43, 44, 45, 46};
    const double Var2 [4*6] = {11, 12, 13, 14, 15, 16,
                               21, 22, 23, 24, 25, 26,
                               31, 32, 33, 34, 35, 36,
                               41, 42, 43, 44, 45, 46};
    
    const double*const vars [vars_num_] = {Var1, Var2};
    newDS.WriteVars(vars);
    
    return 0;
}
