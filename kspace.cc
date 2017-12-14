// kspace.cc
/* Description */

#include <iostream>
#include <complex> // Needed to import alloc_dealloc.h
#include "alloc_dealloc.h" // Allocation/deallocation of arrays
#include "init_routines.h" // Initialization of arrays
#include "kspace.h" // Include header file for consistency check


/* Constructor implementation. To avoid naming ambiguities, we end the argument variables 
with an underscore. */
kspace_t::kspace_t(const int kx_pts_, const double*const kx_bounds_, 
                   const int ky_pts_, const double*const ky_bounds_, 
                   const int bands_num_)
    :kx_pts(kx_pts_), kx_bounds(kx_bounds_), kx_grid(new double [kx_pts_]),
    ky_pts(ky_pts_), ky_bounds(ky_bounds_), ky_grid(new double [ky_pts_]),
    bands_num(bands_num_), energies(Alloc3D_d(kx_pts_, ky_pts_, bands_num_))
{
    std::cout << "kspace_t instance created." << std::endl;
    // The momentum grids are initialized
    LinInitArray(kx_bounds[0], kx_bounds[1], kx_pts, kx_grid);
    LinInitArray(ky_bounds[0], ky_bounds[1], ky_pts, ky_grid);
    // The energies grid is initialized to zero
    ZeroInitArray(kx_pts*ky_pts*bands_num, energies[0][0]);
}

// Destructor implementation
kspace_t::~kspace_t()
{
    delete [] kx_grid;
    delete [] ky_grid;
    Dealloc3D(energies, kx_pts);
    std::cout << "kspace_t instance deleted." << std::endl;
}