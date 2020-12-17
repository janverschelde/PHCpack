// The file random_polynomials.cpp defines functions specified
// in random_polynomials.h.

#include <cstdlib>
#include "random_numbers.h"
#include "random_monomials.h"

void make_supports ( int dim, int nbr, int *nvr )
{
   int rnd;

   for(int i=0; i<nbr; i++)
   {
      rnd = rand() % dim;  // in range 0..dim-1
      nvr[i] = 1 + rnd;    // in range 1..dim
   }
}

bool make_real_polynomial
 ( int dim, int nbr, int pwr, int deg, int *nvr, int **idx, int **exp,
   double *cst, double **cff )
{
   bool fail = false;

   for(int i=0; i<=deg; i++) cst[i] = random_double();

   for(int i=0; i<nbr; i++)
   {
      fail = make_real_monomial(dim,nvr[i],pwr,deg,idx[i],exp[i],cff[i]);
      if(fail) return true;
   }
   return fail;
}
