// ham1.h
#ifndef HAM1_H
#define HAM1_H

class params_t
{
private:
public:
    // Members (pointers in particular) are initialized in the constructor
    const int kx_pts;
    const double*const kx_bounds;
    const int ky_pts;
    const double*const ky_bounds;
    
    const int bands_num;
    const int ham_array_rows;
    const int ham_array_cols;
    
    params_t(); // Constructor declaration
};

void Set_params(params_t& params);

void HamEval(const double kx, const double ky, const double M, 
             std::complex<double>** ham_array);

double Evaluate_M(const int kx_pts, const int ky_pts, const int bands_num, 
                 const double mu, double*** evals);

#endif