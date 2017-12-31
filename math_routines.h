// math_routines.h
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/
#ifndef MATH_ROUTINES_H
#define MATH_ROUTINES_H

double nF0(double const energy);
float  nF0(float  const energy);
double nF(const double beta, const double energy);
float  nF(const float  beta, const float  energy);

double FermiEnerg(const int num_states, const int filled_states, 
                  double const *const energies, const bool print_bounds=false);

#endif
