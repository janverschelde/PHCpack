// The file dbl_monomial_systems.cpp defines the functions specified in
// the file dbl_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "random_numbers.h"
#include "random_series.h"
#include "dbl_convolutions_host.h"
#include "dbl_monomial_systems.h"

using namespace std;

void make_real_exponentials ( int dim, int  deg, double **s )
{
   double rnd;

   const double fac = 64.0;      // 1024.0;
   const double inc = 63.0/64.0; // 1023.0/1024.0;

   for(int i=0; i<dim; i++)
   {
      rnd = random_double()/fac;     // rnd in [-1/fac, +1/fac]
      
      if(rnd < 0)         // if -1/fac <= rnd       < 0 
         rnd = rnd - inc; // then   -1 <= rnd - inc < -inc
      else                // if     0  <= rnd       <= 1/fac
         rnd = rnd + inc; // then  inc <= rnd + inc <= 1

      dbl_exponential(deg,rnd,s[i]);
   }
}

void make_complex_exponentials
 ( int dim, int deg, double *angles, double **sre, double **sim )
{
   double xre,xim;

   for(int i=0; i<dim; i++)
   {
      angles[i] = random_angle();
      xre = cos(angles[i]);
      xim = sin(angles[i]);

      cmplx_exponential(deg,xre,xim,sre[i],sim[i]);
   }
}

void evaluate_real_monomials
 ( int dim, int deg, int **rowsA, double **x, double **rhs )
{
   const int degp1 = deg+1;

   double *acc = new double[degp1]; // accumulates product

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_dbl_product(deg,x[j],rhs[i],acc);

               for(int L=0; L<degp1; L++) rhs[i][L] = acc[L];
            }
         }
      }
   }
   free(acc);
}

void evaluate_complex_monomials
 ( int dim, int deg, int **rowsA,
   double **xre, double **xim, double **rhsre, double **rhsim )
{
   const int degp1 = deg+1;

   double *accre = new double[degp1]; // accumulates product
   double *accim = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx_product
                  (deg,xre[j],xim[j],rhsre[i],rhsim[i],accre,accim);

               for(int L=0; L<degp1; L++)
               {
                  rhsre[i][L] = accre[L];
                  rhsim[i][L] = accim[L];
               }
            }
         }
      }
   }
   free(accre); free(accim);
}

void make_real_coefficients ( int nbrcol, int dim, double ***cff )
{
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++) cff[i][j][0] = random_double();
}

void make_complex_coefficients
 ( int nbrcol, int dim, double ***cffre, double ***cffim )
{
   double rnd;

   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++) 
      {
         rnd = random_angle();
         cffre[i][j][0] = cos(rnd);
         cffim[i][j][0] = sin(rnd);
      }
}

void evaluate_real_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cff, double **x, double **rhs, int vrblvl )
{
   const int degp1 = deg+1;

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "evaluating at the series x ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "x[" << i << "][" << j << "] : " << x[i][j] << endl;
   }
   double **prdrhs = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhs[i] = new double[degp1];

      for(int k=0; k<degp1; k++) rhs[i][k] = 0.0; // initialize sum to 0
   }
   rhs[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++)   // initialize product to the coefficient
      {
         prdrhs[i][0] = cff[col][i][0];
         for(int k=1; k<degp1; k++) prdrhs[i][k] = 0.0;
      }
      if(vrblvl > 1)
         cout << "Evaluating at column " << col << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         for(int k=0; k<dim; k++) rowsA[i][k] = 0;
         for(int k=0; k<nvr[col][i]; k++) rowsA[i][idx[col][i][k]] = 1;
         if(vrblvl > 1)
         {
            for(int k=0; k<dim; k++) cout << " " << rowsA[i][k];
            cout << endl;
         }
      }
      evaluate_real_monomials(dim,deg,rowsA,x,prdrhs);
      if(vrblvl > 1)
      {
         cout << scientific << setprecision(16);
         for(int i=0; i<dim; i++)
         {
            cout << "value at dimension " << i << " :" << endl;
            for(int j=0; j<degp1; j++)
               cout << "prdrhs[" << i << "][" << j << "] : "
                    << prdrhs[i][j] << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         // cout << "rowsum[" << i << "] : " << rowsum << endl;
         if(rowsum != 0)
            for(int k=0; k<degp1; k++) rhs[i][k] += prdrhs[i][k];
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "the evaluated series ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : " << rhs[i][j] << endl;
   }
   for(int i=0; i<dim; i++) free(prdrhs[i]);
   free(prdrhs);
}

void evaluate_complex_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffre, double ***cffim, double **xre, double **xim,
   double **rhsre, double **rhsim, int vrblvl )
{
   const int degp1 = deg+1;

   double **prdrhsre = new double*[dim];
   double **prdrhsim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhsre[i] = new double[degp1];
      prdrhsim[i] = new double[degp1];

      for(int k=0; k<degp1; k++)  // initialize sum to zero
      {
         rhsre[i][k] = 0.0; rhsim[i][k] = 0.0;
      }
   }
   rhsre[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++)    // initialize product to coefficient
      {
         prdrhsre[i][0] = cffre[col][i][0]; 
         prdrhsim[i][0] = cffim[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhsre[i][k] = 0.0; prdrhsim[i][k] = 0.0;
         }
      }
      if(vrblvl > 1)
         cout << "Evaluating at column " << col << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         for(int k=0; k<dim; k++) rowsA[i][k] = 0;
         for(int k=0; k<nvr[col][i]; k++) rowsA[i][idx[col][i][k]] = 1;
         if(vrblvl > 1)
         {
            for(int k=0; k<dim; k++) cout << " " << rowsA[i][k];
            cout << endl;
         }
      }
      evaluate_complex_monomials
         (dim,deg,rowsA,xre,xim,prdrhsre,prdrhsim);

      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               rhsre[i][k] += prdrhsre[i][k];
               rhsim[i][k] += prdrhsim[i][k];
            }
      }
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhsre[i]); free(prdrhsim[i]);
   }
   free(prdrhsre); free(prdrhsim);
}
