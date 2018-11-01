// chempot.cc
/*
*/

#include <iostream> // Input/output to command line
#include <cmath> // Numerous math functions
#include "math_routines.h"
#include "chempot.h" // Include header file for consistency check

// Internal subroutines

double func(const int num_states, const int num_electrons, 
            const double beta, const double*const energies, const double mu)
{
    /* Function whose root is the physical chemical potential, needed in the numerical 
    root solver.
    func(mu) = - num_electrons + sum (nF(energies-mu))*/
    
    double accumulator = - (double)(num_electrons);
    for (int i=0; i<num_states; ++i)
        accumulator += nF(beta, energies[i] - mu); // nF() is from math_routines
    
    return accumulator;
}

double func_deriv(const int num_states, const int num_electrons, 
                  const double beta, const double*const energies, const double mu)
{
    /* Derivative of func().
    func_deriv(mu) = */
    
    double accumulator = 0.;
    for (int i=0; i<num_states; ++i)
        accumulator += std::pow( std::cosh(0.5*beta*(energies[i]-mu)), -2);
    
    accumulator *= beta/4.; // Multiply sum by beta/4.
    
    return accumulator;
}


// External routines

double ChemPotNewton(const int num_states, const int num_electrons, 
                     const double*const energies, const double beta, 
                     const double tol, const int loops_lim, const bool show_output)
{
    
    // As a first guess, use the Fermi energy
    double mu = FermiEnerg_cpy(num_states, num_electrons, energies);
    double mu_new = mu;
    
    int fail_counter = 0;
    int  win_counter = 0;
    bool converged = false;
    bool fail = false;
    
    do // Update via Newton's method
    {
      mu = mu_new; // Update mu to mu_new
      
      mu_new = mu - func(num_states, num_electrons, beta, energies, mu)/
                               func_deriv(num_states, num_electrons, beta, energies, mu);
      
      if (show_output)
        std::cout << "mu = " << mu  << ", mu_new = " << mu_new 
                  << "\tmu_new-mu = " << mu_new-mu << std::endl;
      
      if (std::abs(mu_new-mu)<tol)
      {
        fail_counter = 0;
        win_counter += 1;
      }
      else
      {
        fail_counter += 1;
        win_counter   = 0;
      }
      
      if (win_counter>=2)
        converged = true;
      if (fail_counter>=loops_lim)
      {
        fail = true;
        std::cerr << " WARNING: chemical potential not converged." << std::endl;
      }
      
    } while (!converged && !fail);
    
    return mu;
}

double ChemPotBisec(const int num_states, const int num_electrons, 
                    const double*const energies, const double beta, const double tol, 
                    const bool show_output)
{
    /* The chemical potential is found using the bisection method. Note that this doesn't 
    work if num_electrons is 0 or num_states, because in that case the function has no 
    root (it tends to zero asymptotically in mu). */
    
    // Find min and max energies
    double minE=0., maxE=0.;
    MinMaxArray(energies, num_states, minE, maxE); // minE and maxE are assigned
    
    /* Starting values for a and b. Choose them slightly outside the energy range to be 
    safe. */
    double a = minE - (maxE-minE)*0.05;
    double b = maxE + (maxE-minE)*0.05;
    
    int counter = 0;
    bool converged = false;
    
    do // Update via Newton's method
    {
      counter++; // Increment 'counter' by one
      
      const double midpoint = (a+b)/2.; // Get the midpoint between a and b
      // Find the image of the function at the midpoint
      const double image = func(num_states, num_electrons, beta, energies, midpoint);
      
      if (show_output)
          std::cout << "a = " << a << ", b = " << b << "\t"
                    << "midpoint = " << midpoint << ", image = " << image << std::endl;
      
      if (image==0.)
      {
        a = midpoint; // If image is 0, the loop terminates and the value is returned.
        b = midpoint;
      }
        
      else if (image<0.)
        a = midpoint; // If 'image' is negative, 'midpoint' is the new 'a'
      else if (image>0.)
        b = midpoint; // If 'image' is positive, 'midpoint' is the new 'b'
      
      if (std::abs((b-a)/2.)<tol)
        converged = true; // Check if the half difference is below tolerance
      
    } while (!converged); // The looping stops if convergence has been reached.
    
    return (b+a)/2.; // Take mu to be the average of a and b.
}
