// ham1.h
#ifndef HAM1_H
#define HAM1_H

void HamEval(const double kx, const double ky, const double M, 
             std::complex<double>** ham_array);

double EvaluateM(const int kx_pts, const int ky_pts, const int bands_num, 
                 const double mu, double*** evals);

#endif