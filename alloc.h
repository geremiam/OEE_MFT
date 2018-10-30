// alloc.h
/* Module containing routines for allocating and deallocating arrays of varying 
   dimensionality. */
#ifndef ALLOC_H
#define ALLOC_H


double*** Alloc4D_d(const int Dim0Len, const int Dim1Len, const int Dim2Len, const int Dim3Len);

float***  Alloc4D_f(const int Dim0Len, const int Dim1Len, const int Dim2Len, const int Dim3Len);

double*** Alloc3D_d(const int Dim0Len, const int Dim1Len, const int Dim2Len);
void Dealloc3D(double const *const *const *const Array, const int Dim0Len);

float***  Alloc3D_f(const int Dim0Len, const int Dim1Len, const int Dim2Len);
void Dealloc3D(float const *const *const *const Array, const int Dim0Len);

double** Alloc2D_d(const int NumRows, const int NumCols);
void Dealloc2D(double const *const *const Array);

float** Alloc2D_f(const int NumRows, const int NumCols);
void Dealloc2D(float const *const *const Array);

std::complex<double>** Alloc2D_z(const int NumRows, const int NumCols);
void Dealloc2D(std::complex<double> const *const *const Array);

std::complex<float>** Alloc2D_c(const int NumRows, const int NumCols);
void Dealloc2D(std::complex<float> const *const *const Array);


#endif
