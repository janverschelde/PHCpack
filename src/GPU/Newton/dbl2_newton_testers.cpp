// The file dbl2_newton_testers.cpp defines the functions with prototypes in
// the file dbl2_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_factorizations.h"
#include "dbl2_bals_host.h"

using namespace std;

void dbl2_unit_series_vector
 ( int dim, int deg, double **cffhi, double **cfflo )
{
   for(int i=0; i<dim; i++)
   {
      cffhi[i][0] = 1.0;
      cfflo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
      }
   }
}

void dbl2_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo, int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_dbl2_product
                     (deg,inputhi[idxvar],inputlo[idxvar],
                      cffhi[i],cfflo[i],acchi,acclo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffhi[i][L] = acchi[L];
                     cfflo[i][L] = acclo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffhi[i][j] << "  " << cfflo[i][j] << endl;
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl2_evaldiff
         (dim,nvr[i],deg,idx[i],cffhi[i],cfflo[i],
          inputhi,inputlo,outputhi[i],outputlo[i]);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputhi[i][dim][j] << "  "
                 << outputlo[i][dim][j] << endl;
      }
   }

   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];
               double tmphi,tmplo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // output[i][idxvar][k] = factor*output[i][idxvar][k];
                  ddf_mul_d_dd(factor,outputhi[i][idxvar][k],
                                      outputlo[i][idxvar][k],&tmphi,&tmplo);
                  outputhi[i][idxvar][k] = tmphi;
                  outputlo[i][idxvar][k] = tmplo;
               }
            }
         }
      }
   }
}

void dbl2_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo, 
   double **rhshi, double **rhslo, double ***jacvalhi, double ***jacvallo,
   int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhi[i][j] = outputhi[i][dim][j];
         funvallo[i][j] = outputlo[i][dim][j];
      }

   if(vrblvl > 0)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalhi[i][0] << "  "
                            << funvallo[i][0] << endl;
   }
   // Linearize the function values in the rhs and swap sign,
   // but keep in mind that the right hand side is 1 - t,
   // so we subtract 1 and add t to the rhs.
   for(int j=0; j<dim; j++) 
   {                                 // rhs[0][j] = -(funval[j][0] - 1.0);
      rhshi[0][j] = -funvalhi[j][0];
      rhslo[0][j] = -funvallo[j][0];
      ddf_inc(&rhshi[0][j],&rhslo[0][j],1.0,0.0);
   }
   if(degp1 > 1)
   {
      for(int j=0; j<dim; j++)       // rhs[1][j] = -(funval[j][1] + 1.0);
      {
         rhshi[1][j] = -funvalhi[j][1];
         rhslo[1][j] = -funvallo[j][1];
         ddf_dec(&rhshi[1][j],&rhslo[1][j],1.0,0.0);
      }
      for(int i=2; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshi[i][j] = -funvalhi[j][i];
            rhslo[i][j] = -funvallo[j][i];
         }
   }
   if(vrblvl > 0)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhshi[i][j] << "  " << rhslo[i][j] << endl;
      }
   }
   // output[i][idx[i][k]] for k from 0 to nvr[i] has the coefficients
   // of the series of the derivative of the i-th monomial with respect
   // to the variable idx[i][k].
   for(int i=0; i<dim; i++)          // the i-th monomial
   {
      int *indexes = idx[i];         // indices of the variables
      for(int k=0; k<nvr[i]; k++)    // derivative w.r.t. idx[i][k]
      {                              // has j-th coefficient
         int idxval = indexes[k];

         for(int j=0; j<degp1; j++) 
         {
            jacvalhi[j][i][idxval] = outputhi[i][idxval][j];
            jacvallo[j][i][idxval] = outputlo[i][idxval][j];
         }
      }
   }
   if(vrblvl > 0)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalhi[0][i][j] << "  " << jacvallo[0][i][j] << endl;
      }
   }
}

void dbl2_update_series
 ( int dim, int degp1, double **xhi, double **xlo,
   double **dxhi, double **dxlo, int vrblvl )
{
   if(vrblvl > 0)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xhi[i][j] << "  " << xlo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         ddf_inc(&xhi[i][j],&xlo[i][j],dxhi[j][i],dxlo[j][i]);
      }

   if(vrblvl > 0)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xhi[i][j] << "  " << xlo[i][j] << endl;
      }
   }
}

void dbl2_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo,
   double ***jacvalhi, double ***jacvallo,
   double **rhshi, double **rhslo, double **solhi, double **sollo,
   double **workmathi, double **workmatlo,
   double *workvechi, double *workveclo,
   double **workrhshi, double **workrhslo,
   double **resvechi, double **resveclo,
   double *resmaxhi, double *resmaxlo, int *ipvt, int vrblvl )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   for(int i=0; i<dim; i++)
   {
      cffhi[i][0] = 1.0;
      cfflo[i][0] = 0.0;

      for(int j=1; j<degp1; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
      }
   }
   dbl2_evaluate_monomials
      (dim,deg,nvr,idx,exp,nbrfac,expfac,
       cffhi,cfflo,acchi,acclo,inputhi,inputlo,outputhi,outputlo,vrblvl);

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhi[i][j][k] = 0.0;
            jacvallo[i][j][k] = 0.0;
         }

   dbl2_linearize_evaldiff_output
      (dim,degp1,nvr,idx,outputhi,outputlo,funvalhi,funvallo,
       rhshi,rhslo,jacvalhi,jacvallo,vrblvl);

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         workrhshi[i][j] = rhshi[i][j];
         workrhslo[i][j] = rhslo[i][j];
      }
   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhi[i][j] = 0.0;
         sollo[i][j] = 0.0;
      }
 
   CPU_dbl2_lusb_solve
      (dim,degp1,jacvalhi,jacvallo,workrhshi,workrhslo,solhi,sollo,
       workmathi,workmatlo,workvechi,workveclo,ipvt,0); // vrblvl);

   CPU_dbl2_linear_residue
      (dim,degp1,jacvalhi,jacvallo,rhshi,rhslo,solhi,sollo,
       resvechi,resveclo,resmaxhi,resmaxlo,vrblvl);

   if(vrblvl > 0)
      cout << "maximum residual : "
           << *resmaxhi << "  " << *resmaxlo << endl;

   dbl2_update_series(dim,degp1,inputhi,inputlo,solhi,sollo,vrblvl);
}

void dbl2_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h,
   double ***outputhi_d, double ***outputlo_d,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d,
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double **workmathi, double **workmatlo,
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   for(int i=0; i<dim; i++)
   {
      cffhi[i][0] = 1.0;
      cfflo[i][0] = 0.0;
      for(int j=1; j<degp1; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
      }
   }
   dbl2_evaluate_monomials
      (dim,deg,nvr,idx,exp,nbrfac,expfac,
       cffhi,cfflo,acchi,acclo,inputhi_h,inputlo_h,outputhi_h,outputlo_h,
       vrblvl);

   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhi_h[i][j][k] = 0.0;
            jacvallo_h[i][j][k] = 0.0;
            jacvalhi_d[i][j][k] = 0.0;
            jacvallo_d[i][j][k] = 0.0;
         }

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   dbl2_linearize_evaldiff_output
      (dim,degp1,nvr,idx,outputhi_h,outputlo_h,funvalhi_h,funvallo_h,
       rhshi_h,rhslo_h,jacvalhi_h,jacvallo_h,vrblvl);

   if(vrblvl > 0) cout << "saving the original rhs ..." << endl;

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhshi_h[i][j] = rhshi_h[i][j];
         urhslo_h[i][j] = rhslo_h[i][j];
      }

   if(vrblvl > 0) cout << "initializing the solution ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhi_h[i][j] = 0.0;
         sollo_h[i][j] = 0.0;
      }
 
   if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_solve ..." << endl;

   CPU_dbl2_qrbs_solve
      (dim,degp1,jacvalhi_h,jacvallo_h,urhshi_h,urhslo_h,solhi_h,sollo_h,
       workmathi,workmatlo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,workvechi,workveclo,
       vrblvl);

   if(vrblvl > 0) cout << "calling CPU_dbl2_linear_residue ..." << endl;

   CPU_dbl2_linear_residue
      (dim,degp1,jacvalhi_h,jacvallo_h,rhshi_h,rhslo_h,solhi_h,sollo_h,
       resvechi,resveclo,resmaxhi,resmaxlo,vrblvl);

   if(vrblvl > 0)
      cout << "maximum residual : "
           << *resmaxhi << "  " << *resmaxlo << endl;

   dbl2_update_series(dim,degp1,inputhi_h,inputlo_h,solhi_h,sollo_h,vrblvl);
}
