// math_routines.h
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/
#ifndef MATH_ROUTINES_H
#define MATH_ROUTINES_H

void MinMaxArray(const double* const Array, const int ArrayLen, double& min, double& max);

double MaxArrayValue(const double* const Array, const int ArrayLen);

double nF0(double const energy);
float  nF0(float  const energy);
double nF(const double beta, const double energy);
float  nF(const float  beta, const float  energy);

double FermiEnerg_cpy(const int num_states, const int filled_states, 
                      double const *const energies, const bool print_bounds=false);
double FermiEnerg(const int num_states, const int filled_states, 
                  double*const energies, const bool print_diagnostics=false);

#endif
