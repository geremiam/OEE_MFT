// init_routines.h
/* This module contains routines for initializing arrays to specific values. */
#ifndef INIT_ROUTINES_H
#define INIT_ROUTINES_H

double LinInitArray(const double xMin, const double xMax, const int NumPts, double* Array, 
                 const bool endpoint=false);
float  LinInitArray(const float  xMin, const float  xMax, const int NumPts, float*  Array, 
                 const bool endpoint=false);

void ZeroInitArray(const int NumPts, double* Array);
void ZeroInitArray(const int NumPts, float*  Array);

#endif