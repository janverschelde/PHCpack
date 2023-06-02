/* Tests operations on Newton for series in quad double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "cyclic_columns.h"
#include "cyclic_indices.h"
#include "unimodular_matrices.h"
#include "prompt_newton_setup.h"
#include "dbl4_newton_testers.h"
#include "dbl4_newton_method.h"
#include "cmplx4_newton_method.h"

using namespace std;

int dbl4_test_columns
 ( int dim, int deg, int size, int vrblvl, int nbritr, int nbrcol,
   int nbsteps, int szt, int nbt, int mode, int cdata );
/*
 * DESCRIPTION :
 *   Tests Newton's method using the column representation to
 *   evaluate and differentiate the polynomials in a system.
 *
 * ON ENTRY :
 *   dim      number of equations and variables in the system;
 *   deg      degree at which the series are truncated;
 *   size     size of the numbers in the exponent matrix;
 *   vrblvl   verbose level;
 *   nbritr   number of iterations to generate an exponent matrix,
 *            or the type of input system;
 *   nbrcol   number of column of monomials in the system;
 *   nbsteps  maximum number of steps with Newton's method;
 *   szt      size of each block of threads (tile);
 *   nbt      number of blocks (tiles) for linear algebra;
 *   mode     0 for host, 1 for device, 2 for host + device run;
 *   cdata    0 if on real data, 1 if on complex data. */

int dbl4_test_rows
 ( int dim, int deg, int vrblvl,
   int nbsteps, int szt, int nbt, int mode, int cdata );
/*
 * DESCRIPTION :
 *   Tests Newton's method using the column representation to
 *   evaluate and differentiate the polynomials in a system.
 *
 * ON ENTRY :
 *   dim      number of equations and variables in the system;
 *   deg      degree at which the series are truncated;
 *   vrblvl   verbose level;
 *   nbsteps  maximum number of steps with Newton's method;
 *   szt      size of each block of threads (tile);
 *   nbt      number of blocks (tiles) for linear algebra;
 *   mode     0 for host, 1 for device, 2 for host + device run;
 *   cdata    0 if on real data, 1 if on complex data. */

int main ( void )
{
   cout << "testing Newton in quad double precision..." << endl;

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

   int fail;

   if(nbritr == -5)
      fail = dbl4_test_rows(dim,deg,vrblvl,nbsteps,szt,nbt,mode,cdata);
   else
      fail = dbl4_test_columns(dim,deg,size,vrblvl,nbritr,nbrcol,
                               nbsteps,szt,nbt,mode,cdata);
   
   return fail;
}

int dbl4_test_columns
 ( int dim, int deg, int size, int vrblvl, int nbritr, int nbrcol,
   int nbsteps, int szt, int nbt, int mode, int cdata )
{
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
   // cout << "-> give the damper (1.0 is the default) : "; cin >> dpr;
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
   double dpr = 1.0;

   if(cdata == 0)
      test_dbl4_column_newton
         (szt,nbt,dim,deg,nbrcol,colnvr,colidx,exp,nbrfac,expfac,rowsA,
          dpr,nbsteps,mode,vrblvl);
   else
      test_cmplx4_column_newton
         (szt,nbt,dim,deg,nbrcol,colnvr,colidx,exp,nbrfac,expfac,rowsA,
          dpr,nbsteps,mode,vrblvl);

   return 0;
}

int dbl4_test_rows
 ( int dim, int deg, int vrblvl,
   int nbsteps, int szt, int nbt, int mode, int cdata )
{
   int *nbr = new int[dim];     // number of monomials in each polynomial
   int **nvr = new int*[dim];   // number of variables in dim polynomials
   int ***idx = new int**[dim]; // we have dim polynomials

   make_polynomial_indices(dim,nbr,nvr,idx);
   write_polynomial_indices(dim,nbr,nvr,idx);

   double dpr = 1.0;

   if(cdata == 0)
      test_dbl4_row_newton
         (szt,nbt,dim,deg,nbr,nvr,idx,dpr,nbsteps,mode,vrblvl);
   else
      test_cmplx4_row_newton
         (szt,nbt,dim,deg,nbr,nvr,idx,dpr,nbsteps,mode,vrblvl);

   return 0;
}
