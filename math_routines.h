// math_routines.h
/*
Module containing various math functions or algorithms used to be used in other 
applications.
*/
#ifndef MATH_ROUTINES_H
#define MATH_ROUTINES_H

bool nF0(double const energy);
bool nF0(float  const energy);
double nF(const double beta, const double energy);
float  nF(const float  beta, const float  energy);

#endif