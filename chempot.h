// chempot.h
/*  */
#ifndef CHEMPOT_H
#define CHEMPOT_H


double ChemPotNewton(const int num_states, const int num_electrons, 
                     const double*const energies, const double T, 
                     const double tol, const int loops_lim, const bool show_output=false);

double ChemPotBisec(const int num_states, const int num_electrons, 
                    const double*const energies, const double T, const double tol=1.e-13, 
                    const bool show_output=false);

#endif
