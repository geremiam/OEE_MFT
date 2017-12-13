// diag_routines.h
/* Description */
#ifndef DIAG_ROUTINES_H
#define DIAG_ROUTINES_H

int simple_zheev(const int MatOrder, std::complex<double>*const a, double*const w, 
                 const bool OutputEvecs=false);

#endif