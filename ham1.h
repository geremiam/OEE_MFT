// ham1.h
#ifndef HAM1_H
#define HAM1_H

struct params_t
{
    int kx_pts=0;
    const double* kx_bounds=NULL;
    int ky_pts=0;
    const double* ky_bounds=NULL;
    
    int bands_num = 0;
    int ham_array_rows = ham_array_rows=0;
    int ham_array_cols = ham_array_cols=0;
};

void Set_params(params_t& params);

void HamEval(const double kx, const double ky, const double M, 
             std::complex<double>** ham_array);

double Evaluate_M(const int kx_pts, const int ky_pts, const int bands_num, 
                 const double mu, double*** evals);

#endif