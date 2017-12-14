// kspace.h
/* Description
Note that the object file must be linked to the alloc_dealloc.o and init_routines.o 
object files.*/
#ifndef KSPACE_H
#define KSPACE_H

class kspace_t
{
private:
    // Private copy constructor (prohibits copy creation)
    kspace_t(const kspace_t&);
    // Private assignment operator (prohibits assignment)
    const kspace_t& operator=(const kspace_t&);
public:
    // Declare quantities relevant to the momentum space
    const int kx_pts;
    const double*const kx_bounds; // Two-component array
    double*const kx_grid; // kx coordinate variable
    
    const int ky_pts;
    const double*const ky_bounds; // Two-component array
    double*const ky_grid; // ky coordinate variable
    
    // Declare quantities relevant to the dispersion
    const int bands_num;
    double*const*const*const energies; // energy variable
    
    // Constructor declaration
    kspace_t(const int kx_pts_, const double*const kx_bounds_, 
             const int ky_pts_, const double*const ky_bounds_, 
             const int bands_num_);
    ~kspace_t(); // Destructor declaration
};

#endif