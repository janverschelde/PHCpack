// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "double_double_functions.h"
#include "random2_vectors.h"
#include "random2_series.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomial_systems.h"

using namespace std;

void make_real2_exponentials
 ( int dim, int  deg, double **shi, double **slo )
{
   double rndhi,rndlo,acchi,acclo;

   const double fac = 64.0;      // 1024.0;
   const double inc = 63.0/64.0; // 1023.0/1024.0;

   for(int i=0; i<dim; i++)
   {
      random_double_double(&rndhi,&rndlo);        // rnd in [-1, +1]
      ddf_div(rndhi,rndlo,fac,0.0,&acchi,&acclo); // acc in [-1/fac, +1/fac]
      
      if(rndhi < 0) // if -1/fac <= rnd       < 0
      {             // then   -1 <= rnd - inc < -inc
         ddf_sub(acchi,acclo,inc,0.0,&rndhi,&rndlo);
      }
      else          // if    0  <= rnd       <= 1/fac
      {             // then inc <= rnd + inc <= 1
         ddf_add(acchi,acclo,inc,0.0,&rndhi,&rndlo);
      }
      dbl2_exponential(deg,rndhi,rndlo,shi[i],slo[i]);
   }
}

void make_complex2_exponentials
 ( int dim, int deg,
   double **srehi, double **srelo, double **simhi, double **simlo )
{
   double xrehi,xrelo,ximhi,ximlo,yhi,ylo;

   for(int i=0; i<dim; i++)
   {
      random_double_double(&xrehi,&xrelo); // cosine of some angle
 
      ddf_sqr(xrehi,xrelo,&yhi,&ylo);        // y = cos^2
      ddf_minus(&yhi,&ylo);                  // y = -cos^2
      ddf_inc_d(&yhi,&ylo,1.0);              // y = 1 - cos^2
      ddf_sqrt(yhi,ylo,&ximhi,&ximlo);       // sin is sqrt(1-cos^2)

      cmplx2_exponential
         (deg,xrehi,xrelo,ximhi,ximlo,srehi[i],srelo[i],simhi[i],simlo[i]);
   }
}

void evaluate_real2_monomials
 ( int dim, int deg, int **rowsA,
   double **shi, double **slo, double **rhshi, double **rhslo )
{
   const int degp1 = deg+1;

   double *acchi = new double[degp1]; // accumulates product
   double *acclo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_dbl2_product
                  (deg,shi[j],slo[j],rhshi[i],rhslo[i],acchi,acclo);

               for(int L=0; L<degp1; L++)
               {
                  rhshi[i][L] = acchi[L];
                  rhslo[i][L] = acclo[L];
               }
            }
         }
      }
   }
   free(acchi); free(acclo);
}

void evaluate_complex2_monomials
 ( int dim, int deg, int **rowsA,
   double **srehi, double **srelo, double **simhi, double **simlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo )
{
   const int degp1 = deg+1;

   double *accrehi = new double[degp1]; // accumulates product
   double *accrelo = new double[degp1];
   double *accimhi = new double[degp1];
   double *accimlo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx2_product
                  (deg,srehi[j],srelo[j],simhi[j],simlo[j],
                   rhsrehi[i],rhsrelo[i],rhsimhi[i],rhsimlo[i],
                   accrehi,accrelo,accimhi,accimlo);

               for(int L=0; L<degp1; L++)
               {
                  rhsrehi[i][L] = accrehi[L];
                  rhsrelo[i][L] = accrelo[L];
                  rhsimhi[i][L] = accimhi[L];
                  rhsimlo[i][L] = accimlo[L];
               }
            }
         }
      }
   }
   free(accrehi); free(accimhi);
   free(accrelo); free(accimlo);
}

void make_real2_coefficients
 ( int nbrcol, int dim, double ***cffhi, double ***cfflo )
{
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         random_double_double(&cffhi[i][j][0],&cfflo[i][j][0]);
}

void make_complex2_coefficients
 ( int nbrcol, int dim,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo )
{
   double xrehi,xrelo,ximhi,ximlo,yhi,ylo;

   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
      {
         random_double_double(&xrehi,&xrelo); // cosine of some angle
 
         ddf_sqr(xrehi,xrelo,&yhi,&ylo);        // y = cos^2
         ddf_minus(&yhi,&ylo);                  // y = -cos^2
         ddf_inc_d(&yhi,&ylo,1.0);              // y = 1 - cos^2
         ddf_sqrt(yhi,ylo,&ximhi,&ximlo);       // sin is sqrt(1-cos^2)

         cffrehi[i][j][0] = xrehi;
         cffrelo[i][j][0] = xrelo;
         cffimhi[i][j][0] = ximhi;
         cffimlo[i][j][0] = ximlo;
      }
}

void evaluate_real2_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhi, double ***cfflo, double **xhi, double **xlo,
   double **rhshi, double **rhslo, int vrblvl )
{
   const int degp1 = deg+1;

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "evaluating at the series x ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "x[" << i << "][" << j << "] : "
                 << xhi[i][j] << "  " << xlo[i][j] << endl;
   }
   double **prdrhshi = new double*[dim];
   double **prdrhslo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhshi[i] = new double[degp1];
      prdrhslo[i] = new double[degp1];

      for(int k=0; k<degp1; k++) // initialize sum to 0
      {
         rhshi[i][k] = 0.0; rhslo[i][k] = 0.0;
      }
   }
   rhshi[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++) // initialize product to coefficient
      {
         prdrhshi[i][0] = cffhi[col][i][0];
         prdrhslo[i][0] = cfflo[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhshi[i][k] = 0.0; prdrhslo[i][k] = 0.0;
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
      evaluate_real2_monomials(dim,deg,rowsA,xhi,xlo,prdrhshi,prdrhslo);
      if(vrblvl > 1)
      {
         cout << scientific << setprecision(16);
         for(int i=0; i<dim; i++)
         {
            cout << "value at dimension " << i << " :" << endl;
            for(int j=0; j<degp1; j++)
               cout << "prdrhs[" << i << "][" << j << "] : "
                    << prdrhshi[i][j] << "  "
                    << prdrhslo[i][j] << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               // rhs[i][k] += prdrhs[i][k];
               ddf_inc(&rhshi[i][k],&rhslo[i][k],
                       prdrhshi[i][k],prdrhslo[i][k]);
            }
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "the evaluated series ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << rhshi[i][j] << "  " << rhslo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhshi[i]);
      free(prdrhslo[i]);
   }
   free(prdrhshi); free(prdrhslo);
}

void evaluate_complex2_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo,
   double **xrehi, double **xrelo, double **ximhi, double **ximlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   int vrblvl )
{
   const int degp1 = deg+1;

   double **prdrhsrehi = new double*[dim];
   double **prdrhsrelo = new double*[dim];
   double **prdrhsimhi = new double*[dim];
   double **prdrhsimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhsrehi[i] = new double[degp1];
      prdrhsrelo[i] = new double[degp1];
      prdrhsimhi[i] = new double[degp1];
      prdrhsimlo[i] = new double[degp1];

      for(int k=0; k<degp1; k++)  // initialize sum to zero
      {
         rhsrehi[i][k] = 0.0; rhsimhi[i][k] = 0.0;
         rhsrelo[i][k] = 0.0; rhsimlo[i][k] = 0.0;
      }
   }
   rhsrehi[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++) // initialize product to coefficient
      {
         prdrhsrehi[i][0] = cffrehi[col][i][0];
         prdrhsrelo[i][0] = cffrelo[col][i][0];
         prdrhsimhi[i][0] = cffimhi[col][i][0];
         prdrhsimlo[i][0] = cffrelo[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhsrehi[i][k] = 0.0; prdrhsimhi[i][k] = 0.0;
            prdrhsrelo[i][k] = 0.0; prdrhsimlo[i][k] = 0.0;
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
      evaluate_complex2_monomials
         (dim,deg,rowsA,xrehi,xrelo,ximhi,ximlo,
          prdrhsrehi,prdrhsrelo,prdrhsimhi,prdrhsimlo);

      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               // rhsre[i][k] += prdrhsre[i][k];
               ddf_inc(&rhsrehi[i][k],&rhsrelo[i][k],
                       prdrhsrehi[i][k],prdrhsrelo[i][k]);
               // rhsim[i][k] += prdrhsim[i][k];
               ddf_inc(&rhsimhi[i][k],&rhsimlo[i][k],
                       prdrhsimhi[i][k],prdrhsimlo[i][k]);
            }
      }
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhsrehi[i]); free(prdrhsimhi[i]);
      free(prdrhsrelo[i]); free(prdrhsimlo[i]);
   }
   free(prdrhsrehi); free(prdrhsimhi);
   free(prdrhsrelo); free(prdrhsimlo);
}
