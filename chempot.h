// chempot.h
/* Module with routines for calculating the chemical potential of a system of free 
fermions. */
#ifndef CHEMPOT_H
#define CHEMPOT_H

#include <limits> // To get machine precision

const double eps_double = std::numeric_limits<double>::epsilon(); // Used as rel. tol.


double ChemPotNewton(const int num_states, const int num_electrons, 
                     const double*const energies, const double T, 
                     const double tol, const int loops_lim, const bool show_output=false);

double ChemPotBisec(const int num_states, const int num_electrons, 
                    const double*const energies, const double T, 
                    const bool show_output=false, const bool usethreads=false);

#endif
