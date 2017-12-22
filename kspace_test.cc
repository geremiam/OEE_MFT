// kspace_test.cc
#include <iostream>
#include "kspace.h"

void kspace_t_test()
{
    const int kx_pts_ = 10;
    const double kx_bounds_ [2] = {-3., 3.};
    const int ky_pts_ = 20;
    const double ky_bounds_ [2] = {-3., 3.};
    const int bands_num_ = 2;
    
    kspace_t Inst1(kx_pts_, kx_bounds_, ky_pts_, ky_bounds_, bands_num_);
    
    /*
    std::cout << "Inst1.kx_pts = " << Inst1.kx_pts << std::endl;
    std::cout << "Inst1.kx_bounds = " << Inst1.kx_bounds[0] << " " 
                                      << Inst1.kx_bounds[1] << std::endl;
    std::cout << "Inst1.ky_pts = " << Inst1.ky_pts << std::endl;
    std::cout << "Inst1.ky_bounds = " << Inst1.ky_bounds[0] << " " 
                                      << Inst1.ky_bounds[1] << std::endl;
    */
    
    std::cout << "Inst1.kx_grid = " << std::endl;
    for (int i=0; i< kx_pts_; ++i)
        std::cout << Inst1.kx_grid[i] << " ";
    std::cout << std::endl;
    
    std::cout << "Inst1.ky_grid = " << std::endl;
    for (int i=0; i< ky_pts_; ++i)
        std::cout << Inst1.ky_grid[i] << " ";
    std::cout << std::endl;
    
    std::cout << "Inst1.energies = " << std::endl;
    for (int band=0; band<bands_num_; ++band)
    {
        for (int i=0; i<kx_pts_; ++i)
        {
            for (int j=0; j<ky_pts_; ++j)
                std::cout << Inst1.energies[i][j][band] << " ";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

int main()
{
    kspace_t_test();
    return 0;
}
