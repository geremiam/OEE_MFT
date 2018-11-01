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
    
    // Points in momentum space grids. Choose even values to avoid the origin.
    const int ka_pts_;
    const int kb_pts_;
    const int kc_pts_;
    
    const int bands_num_; // Number of bands
    
    const double a_; // Length of translations in each orthogonal direction
    const double b_;
    const double c_;
public:
    double*const ka_grid; // ka coordinate variable
    double*const kb_grid; // kb coordinate variable
    double*const kc_grid; // kc coordinate variable
    /* The energies are held in a 1D array with the understanding that from slowest 
    varying to fastest varying, the indices are those for ka, kb, kc, and band.*/
    double*const energies; // energy variable
    
    // Constructor declaration
    kspace_t(const double a, const double b, const double c, 
                   const int ka_pts, const int kb_pts, const int kc_pts, 
                   const int bands_num);
    ~kspace_t(); // Destructor declaration
    
    int index(const int ka_ind, const int kb_ind, const int kc_ind, const int band_int);
};

#endif
