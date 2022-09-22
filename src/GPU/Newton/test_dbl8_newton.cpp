/* Tests operations on Newton for series in octo double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "prompt_newton_setup.h"
#include "dbl8_newton_testers.h"

using namespace std;

int main ( void )
{
   cout << "testing Newton in octo double precision..." << endl;

   int seed,dim,deg,size,posvals,vrblvl,nbritr,nbsteps,szt,nbt,mode,cdata;

   prompt_newton_setup
      (&seed,&szt,&nbt,&dim,&deg,&size,&posvals,&vrblvl,&mode,
       &nbritr,&nbsteps,&cdata);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   if(nbritr == -1)
      cout << "making the lower triangular unit of dimension " << dim
           << " ..." << endl;
   if(nbritr > 0)
      cout << "generating an integer matrix of dimension " << dim
           << " ..." << endl;
/*
 * 1. defining the unimodular matrix of monomials
 */
   int **rowsA = new int*[dim];  // exponents in the rows
   int *nvr = new int[dim];      // number of variables in each monomial
   int **idx = new int*[dim];    // indexes of variables in each monomial
   int **exp = new int*[dim];    // exponents of the variables
   int *nbrfac = new int[dim];   // number of exponents > 1 in each monomial
   int **expfac = new int*[dim]; // exponents of the common factors

   make_monomial_system
      (dim,size,posvals,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);
/*
 * 2. calling the test function
 */
   if(cdata == 0)
      test_dbl8_real_newton
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,nbsteps,mode,vrblvl);
   else
      test_dbl8_complex_newton
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,nbsteps,mode,vrblvl);

   return 0;
}
