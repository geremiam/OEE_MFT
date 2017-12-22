// init_routines.h
/* This module contains routines for initializing arrays to specific values. */
#ifndef INIT_ROUTINES_H
#define INIT_ROUTINES_H

double LinInitArray(const double xMin, const double xMax, const int NumPts, 
                    double*const Array, const bool endpoint=false);
float  LinInitArray(const float  xMin, const float  xMax, const int NumPts, 
                    float*const  Array, const bool endpoint=false);

void ValInitArray(const int NumPts, double*const Array, const double Value=0.);
void ValInitArray(const int NumPts, float*const  Array, const float  Value=0.f);
void ValInitArray(const int NumPts, std::complex<double>*const Array, 
                  const std::complex<double> Value={0.,0.});
void ValInitArray(const int NumPts, std::complex<float>*const  Array, 
                  const std::complex<float>  Value={0.f,0.f});

#endif
