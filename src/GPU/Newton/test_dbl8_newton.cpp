/* Tests operations on Newton for series in octo double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl8_newton_testers.h"

using namespace std;

int main ( void )
{
   cout << "testing Newton ..." << endl;

   int dim,deg,size,posvals,vrblvl,nbritr,nbsteps,szt,nbt,mode;
   prompt_dimensions(&dim,&deg,&size,&posvals,&vrblvl,&nbritr,&nbsteps);

   cout << "-> enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> mode;

   if(mode != 1)
   {
      cout << "-> give the number of tiles : "; cin >> nbt;
      cout << "-> give the size of each tile : "; cin >> szt;
      int p = szt*nbt;

      while(p != dim)
      {
          cout << "Dimension = " << dim << " != " << szt << " * " << nbt
               << ", retry." << endl;
          cout << "-> give the size of each tile : "; cin >> szt;
          cout << "-> give the number of tiles : "; cin >> nbt;
          p = szt*nbt;
      }
   }

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
   test_dbl8_real_newton
      (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,nbsteps,mode,vrblvl);

   return 0;
}
