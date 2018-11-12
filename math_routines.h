// math_routines.h
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/
#ifndef MATH_ROUTINES_H
#define MATH_ROUTINES_H

#include <complex> // Includes complex numbers

void MinMaxArray(const double* const Array, const int ArrayLen, double& min, double& max);

double MaxArrayValue(const double* const Array, const int ArrayLen);

double nF0(double const energy);
float  nF0(float  const energy);
double nF(const double T, const double energy);
float  nF(const float  T, const float  energy);

double FermiEnerg_cpy(const int num_states, const int filled_states, 
                      double const *const energies, const bool print_bounds=false);
double FermiEnerg(const int num_states, const int filled_states, 
                  double*const energies, const bool print_diagnostics=false);

std::complex<double> TraceMF(const int num_bands, 
                             const std::complex<double>*const*const evecs, 
                             const std::complex<double>*const*const mat, 
                             const double*const occs);

void Occupations(const int arrlen, const double mu, const double*const energies, 
                 double*const occs, const bool zerotemp, const double T=-1.);

#endif
