/* Tests a path tracker on series in octo double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "cyclic_columns.h"
#include "prompt_newton_setup.h"
#include "dbl8_path_tracker.h"
#include "cmplx8_path_tracker.h"

using namespace std;

int main ( void )
{
   cout << "test on path tracking in octo double precision ..." << endl;

   int seed,dim,deg,size,vrblvl,nbritr,nbrcol,nbsteps,szt,nbt,mode,cdata;

   prompt_newton_setup
      (&seed,&szt,&nbt,&dim,&deg,&size,&vrblvl,&mode,
       &nbritr,&nbrcol,&nbsteps,&cdata);

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

   if(nbritr == -3)
   {
      for(int i=0; i<dim; i++) rowsA[i] = new int[dim];
   }
   else
   {
      if(nbritr == -4)           // the 1st column is lower triangular
         make_monomial_system
            (dim,size,1,-1,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);
      else
         make_monomial_system
            (dim,size,1,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);

      int *expsol = new int[dim];
      int sing = exponents_check(dim,rowsA,expsol,vrblvl);
      if(sing != 0)
      {
         cout << "The exponent matrix is singular!  Abort." << endl;
         return sing;
      }
   }
/*
 * 2. calling the test function
 */
   int **colnvr = new int*[nbrcol];
   int ***colidx = new int**[nbrcol];
   for(int i=0; i<nbrcol; i++)
   {
      colnvr[i] = new int[dim];
      colidx[i] = new int*[dim];
   }
   if(nbrcol == 1)  // we have a monomial system
   {
      for(int i=0; i<dim; i++)
      {
         colnvr[0][i] = nvr[i];
         colidx[0][i] = new int[nvr[i]];
         for(int j=0; j<nvr[i]; j++) colidx[0][i][j] = idx[i][j];
      }
   }
   else if(nbrcol == 2) // a 2-column system
   {
      for(int i=0; i<dim; i++)  // first column is lower triangular
      {
         colnvr[0][i] = nvr[i];
         colidx[0][i] = new int[nvr[i]];
         for(int j=0; j<nvr[i]; j++) colidx[0][i][j] = idx[i][j];
      }
      for(int i=0; i<dim; i++)  // second column is upper triangular
      {
         int i2 = dim-i-1;
         colnvr[1][i2] = nvr[i];
         colidx[1][i2] = new int[nvr[i]];
         for(int j=0; j<nvr[i]; j++) colidx[1][i2][j] = idx[i][j];
      }
   }
   else
   {
      make_cyclic_variables(nbrcol,dim,colnvr);

      for(int i=0; i<nbrcol; i++)
         for(int j=0; j<dim; j++) colidx[i][j] = new int[colnvr[i][j]];
     
      make_cyclic_columns(nbrcol,dim,colnvr,colidx);
      if(vrblvl > 1)
      {
          cout << "column representation of cyclic "
               << dim << "-roots :" << endl;
          write_cyclic_columns(nbrcol,dim,colnvr,colidx);
      }
   }
   if(cdata == 0)
      test_dbl8_real_track
         (szt,nbt,dim,deg,nbrcol,colnvr,colidx,exp,nbrfac,expfac,rowsA,
          nbsteps,mode,vrblvl);
   else
      test_dbl8_complex_track
         (szt,nbt,dim,deg,nbrcol,colnvr,colidx,exp,nbrfac,expfac,rowsA,
          nbsteps,mode,vrblvl);

   return 0;
}
