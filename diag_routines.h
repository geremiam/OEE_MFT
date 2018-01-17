// diag_routines.h
/* Routines for finding eigenvalues and eigenvectors of matrices. */
#ifndef DIAG_ROUTINES_H
#define DIAG_ROUTINES_H

int simple_zheev(const int MatOrder, std::complex<double>*const a, double*const w, 
                 const bool OutputEvecs=false, std::complex<double>*const z=NULL);
int simple_zheevr(const int MatOrder, std::complex<double>*const a, double*const w, 
                  const bool OutputEvecs=false, std::complex<double>*const z=NULL);
                 
#endif
