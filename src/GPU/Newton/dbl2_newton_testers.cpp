// The file dbl2_newton_testers.cpp defines the functions with prototypes in
// the file dbl2_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_factorizations.h"
#include "dbl2_bals_host.h"
#include "dbl2_bals_kernels.h"
#include "dbl2_systems_host.h"
#include "dbl2_systems_kernels.h"

using namespace std;

void dbl2_unit_series_vector
 ( int dim, int deg, double **cffhi, double **cfflo )
{
   for(int i=0; i<dim; i++)
   {
      cffhi[i][0] = 1.000005;
      cfflo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
      }
   }
}

void cmplx2_unit_series_vector
 ( int dim, int deg,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehi[i][0] = 1.000005; cffrelo[i][0] = 0.0;
      cffimhi[i][0] = 0.000005; cffimlo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffrehi[i][j] = 0.0; cffrelo[i][j] = 0.0;
         cffimhi[i][j] = 0.0; cffimlo[i][j] = 0.0;
      }
   }
}

void dbl2_update_series
 ( int dim, int degp1, double **xhi, double **xlo,
   double **dxhi, double **dxlo, int vrblvl )
{
   if(vrblvl > 1)
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

void cmplx2_update_series
 ( int dim, int degp1,
   double **xrehi, double **xrelo, double **ximhi, double **ximlo,
   double **dxrehi, double **dxrelo, double **dximhi, double **dximlo,
   int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehi[i][j] << "  " << xrelo[i][j] << endl << "  "
                 << ximhi[i][j] << "  " << ximlo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++)
      {
         // xre[i][j] = xre[i][j] + dxre[j][i];
         ddf_inc(&xrehi[i][j],&xrelo[i][j],dxrehi[j][i],dxrelo[j][i]);
         // xim[i][j] = xim[i][j] + dxim[j][i];
         ddf_inc(&ximhi[i][j],&ximlo[i][j],dximhi[j][i],dximlo[j][i]);
      }

   if(vrblvl > 0)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehi[i][j] << "  " << xrelo[i][j] << endl << "  "
                 << ximhi[i][j] << "  " << ximlo[i][j] << endl;
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
   CPU_dbl2_evaluate_monomials
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
      cout << "maximum residual : " << *resmaxhi << endl;

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
   if((mode == 1) || (mode == 2))
   {
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
      CPU_dbl2_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhi,cfflo,acchi,acclo,inputhi_h,inputlo_h,outputhi_h,outputlo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)  // reset the coefficients
      {
         cffhi[i][0] = 1.0;
         cfflo[i][0] = 0.0;
         for(int j=1; j<degp1; j++)
         {
            cffhi[i][j] = 0.0;
            cfflo[i][j] = 0.0;
         }
      }
      GPU_dbl2_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,cffhi,cfflo,acchi,acclo,
          inputhi_d,inputlo_d,outputhi_d,outputlo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;
      for(int k=0; k<dim; k++) // monomial k
         for(int i=0; i<=dim; i++)
            for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "output_h[" << k << "][" << i << "][" << j << "] : "
                    << outputhi_h[k][i][j] << "  "
                    << outputlo_h[k][i][j] << endl;
               cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                    << outputhi_d[k][i][j] << "  "
                    << outputlo_d[k][i][j] << endl;
            }
            errsum += abs(outputhi_h[k][i][j] - outputhi_d[k][i][j])
                    + abs(outputlo_h[k][i][j] - outputlo_d[k][i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }

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

   if((mode == 1) || (mode == 2))
      dbl2_linearize_evaldiff_output
         (dim,degp1,nvr,idx,outputhi_h,outputlo_h,funvalhi_h,funvallo_h,
          rhshi_h,rhslo_h,jacvalhi_h,jacvallo_h,vrblvl);
   if((mode == 0) || (mode == 2))
      dbl2_linearize_evaldiff_output
         (dim,degp1,nvr,idx,outputhi_d,outputlo_d,funvalhi_d,funvallo_d,
          rhshi_d,rhslo_d,jacvalhi_d,jacvallo_d,vrblvl);

   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "funval_h[" << i << "][" << j << "] : "
                    << funvalhi_h[i][j] << "  " << funvallo_h[i][j] << endl;
               cout << "funval_d[" << i << "][" << j << "] : "
                    << funvalhi_d[i][j] << "  " << funvallo_d[i][j] << endl;
            }
            errsum += abs(funvalhi_h[i][j] - funvalhi_d[i][j])
                    + abs(funvallo_h[i][j] - funvallo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            for(int k=0; k<dim; k++)
            {
               if(vrblvl > 1)
               {
                  cout << "jacval_h[" << i << "][" << j << "][" << k << "] : "
                       << jacvalhi_h[i][j][k] << "  "
                       << jacvallo_h[i][j][k] << endl;
                  cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                       << jacvalhi_d[i][j][k] << "  "
                       << jacvallo_d[i][j][k] << endl;
               }
               errsum += abs(jacvalhi_h[i][j][k] - jacvalhi_d[i][j][k])
                       + abs(jacvallo_h[i][j][k] - jacvallo_d[i][j][k]);
            }
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "rhs_h[" << i << "][" << j << "] : "
                    << rhshi_h[i][j] << "  " << rhslo_h[i][j] << endl;
               cout << "rhs_d[" << i << "][" << j << "] : "
                    << rhshi_d[i][j] << "  " << rhslo_d[i][j] << endl;
            }
            errsum += abs(rhshi_h[i][j] - rhshi_d[i][j])
                    + abs(rhslo_h[i][j] - rhslo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }

   if(vrblvl > 0) cout << "saving the original rhs ..." << endl;

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhshi_h[i][j] = rhshi_h[i][j];
         urhslo_h[i][j] = rhslo_h[i][j];
         urhshi_d[i][j] = rhshi_d[i][j];
         urhslo_d[i][j] = rhslo_d[i][j];
      }

   if(vrblvl > 0) cout << "initializing the solution ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhi_h[i][j] = 0.0; sollo_h[i][j] = 0.0;
         solhi_d[i][j] = 0.0; sollo_d[i][j] = 0.0;
      }
 
   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_solve ..." << endl;

      CPU_dbl2_qrbs_solve
         (dim,degp1,jacvalhi_h,jacvallo_h,urhshi_h,urhslo_h,solhi_h,sollo_h,
          workmathi,workmatlo,Qhi_h,Qlo_h,Rhi_h,Rlo_h,workvechi,workveclo,
          vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl2_linear_residue ..." << endl;

         CPU_dbl2_linear_residue
            (dim,degp1,jacvalhi_h,jacvallo_h,rhshi_h,rhslo_h,solhi_h,sollo_h,
             resvechi,resveclo,resmaxhi,resmaxlo,vrblvl);
   
         cout << "maximum residual : " << *resmaxhi << endl;
      }
      dbl2_update_series(dim,degp1,inputhi_h,inputlo_h,solhi_h,sollo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling GPU_dbl2_bals_solve ..." << endl;

      GPU_dbl2_bals_solve
         (dim,degp1,szt,nbt,jacvalhi_d,jacvallo_d,Qhi_d,Qlo_d,Rhi_d,Rlo_d,
          urhshi_d,urhslo_d,solhi_d,sollo_d,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl2_linear_residue ..." << endl;

         CPU_dbl2_linear_residue
            (dim,degp1,jacvalhi_d,jacvallo_d,rhshi_d,rhslo_d,solhi_d,sollo_d,
             resvechi,resveclo,resmaxhi,resmaxlo,vrblvl);
   
         cout << "maximum residual : " << *resmaxhi << endl;
      }
      dbl2_update_series(dim,degp1,inputhi_d,inputlo_d,solhi_d,sollo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "Q_h[" << i << "][" << j << "] : "
                    << Qhi_h[i][j] << "  " << Qlo_h[i][j] << endl;
               cout << "Q_d[" << i << "][" << j << "] : "
                    << Qhi_d[i][j] << "  " << Qlo_d[i][j] << endl;
            }
            errsum += abs(Qhi_h[i][j] - Qhi_d[i][j])
                    + abs(Qlo_h[i][j] - Qlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "R_h[" << i << "][" << j << "] : "
                    << Rhi_h[i][j] << "  " << Rlo_h[i][j] << endl;
               cout << "R_d[" << i << "][" << j << "] : "
                    << Rhi_d[i][j] << "  " << Rlo_d[i][j] << endl;
            }
            errsum += abs(Rhi_h[i][j] - Rhi_d[i][j])
                    + abs(Rlo_h[i][j] - Rlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "urhs_h[" << i << "][" << j << "] : "
                    << urhshi_h[i][j] << "  " << urhslo_h[i][j] << endl;
               cout << "urhs_d[" << i << "][" << j << "] : "
                    << urhshi_d[i][j] << "  " << urhslo_d[i][j] << endl;
            }
            errsum += abs(urhshi_h[i][j] - urhshi_d[i][j])
                    + abs(urhslo_h[i][j] - urhslo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "sol_h[" << i << "][" << j << "] : "
                    << solhi_h[i][j] << "  " << sollo_h[i][j] << endl;
               cout << "sol_d[" << i << "][" << j << "] : "
                    << solhi_d[i][j] << "  " << sollo_d[i][j] << endl;
            }
            errsum += abs(solhi_h[i][j] - solhi_d[i][j])
                    + abs(sollo_h[i][j] - sollo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "input_h[" << i << "][" << j << "] : "
                    << inputhi_h[i][j] << "  "
                    << inputlo_h[i][j] << endl;
               cout << "input_d[" << i << "][" << j << "] : "
                    << inputhi_d[i][j] << "  "
                    << inputlo_d[i][j] << endl;
            }
            errsum += abs(inputhi_h[i][j] - inputhi_d[i][j])
                    + abs(inputlo_h[i][j] - inputlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}

void cmplx2_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double *accrehi, double *accrelo, double *accimhi, double *accimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d, 
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h,
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d,
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h, 
   double **solrehi_d, double **solrelo_d, 
   double **solimhi_d, double **solimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h,
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d, 
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double **workmatrehi, double **workmatrelo,
   double **workmatimhi, double **workmatimlo,
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo, int vrblvl, int mode )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         cffrehi[i][0] = 1.0; cffrelo[i][0] = 0.0;
         cffimhi[i][0] = 0.0; cffimlo[i][0] = 0.0;

         for(int j=1; j<degp1; j++)
         {
            cffrehi[i][j] = 0.0; cffrelo[i][j] = 0.0;
            cffimhi[i][j] = 0.0; cffimlo[i][j] = 0.0;
         }
      }
      CPU_cmplx2_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffrehi,cffrelo,cffimhi,cffimlo,
          accrehi,accrelo,accimhi,accimlo,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)  // reset the coefficients
      {
         cffrehi[i][0] = 1.0; cffrelo[i][0] = 0.0;
         cffimhi[i][0] = 0.0; cffimlo[i][0] = 0.0;

         for(int j=1; j<degp1; j++)
         {
            cffrehi[i][j] = 0.0; cffrelo[i][j] = 0.0;
            cffimhi[i][j] = 0.0; cffimlo[i][j] = 0.0;
         }
      }
      GPU_cmplx2_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffrehi,cffrelo,cffimhi,cffimlo,
          accrehi,accrelo,accimhi,accimlo,
          inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;
      for(int k=0; k<dim; k++) // monomial k
         for(int i=0; i<=dim; i++)
            for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "output_h[" << k << "][" << i << "][" << j << "] : "
                    << outputrehi_h[k][i][j] << "  "
                    << outputrelo_h[k][i][j] << endl << "  "
                    << outputimhi_h[k][i][j] << "  "
                    << outputimlo_h[k][i][j] << endl;
               cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                    << outputrehi_d[k][i][j] << "  "
                    << outputrelo_d[k][i][j] << endl << "  "
                    << outputimhi_d[k][i][j] << "  "
                    << outputimlo_d[k][i][j] << endl;
            }
            errsum += abs(outputrehi_h[k][i][j] - outputrehi_d[k][i][j])
                    + abs(outputrelo_h[k][i][j] - outputrelo_d[k][i][j])
                    + abs(outputimhi_h[k][i][j] - outputimhi_d[k][i][j])
                    + abs(outputimlo_h[k][i][j] - outputimlo_d[k][i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalrehi_h[i][j][k] = 0.0; jacvalimhi_h[i][j][k] = 0.0;
            jacvalrelo_h[i][j][k] = 0.0; jacvalimlo_h[i][j][k] = 0.0;
            jacvalrehi_d[i][j][k] = 0.0; jacvalimhi_d[i][j][k] = 0.0;
            jacvalrelo_d[i][j][k] = 0.0; jacvalimlo_d[i][j][k] = 0.0;
         }

   if((mode == 1) || (mode == 2))
      cmplx2_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
          funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
          rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,vrblvl);

   if((mode == 0) || (mode == 2))
      cmplx2_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
          funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
          rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,vrblvl);

   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "funval_h[" << i << "][" << j << "] : "
                    << funvalrehi_h[i][j] << "  "
                    << funvalrelo_h[i][j] << endl << "  "
                    << funvalimhi_h[i][j] << "  "
                    << funvalimlo_h[i][j] << endl;
               cout << "funval_d[" << i << "][" << j << "] : "
                    << funvalrehi_d[i][j] << "  "
                    << funvalrelo_d[i][j] << endl << "  "
                    << funvalimhi_d[i][j] << "  "
                    << funvalimlo_d[i][j] << endl;
            }
            errsum += abs(funvalrehi_h[i][j] - funvalrehi_d[i][j])
                    + abs(funvalrelo_h[i][j] - funvalrelo_d[i][j])
                    + abs(funvalimhi_h[i][j] - funvalimhi_d[i][j])
                    + abs(funvalimlo_h[i][j] - funvalimlo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            for(int k=0; k<dim; k++)
            {
               if(vrblvl > 1)
               {
                  cout << "jacval_h[" << i << "][" << j << "][" << k << "] : "
                       << jacvalrehi_h[i][j][k] << "  "
                       << jacvalrelo_h[i][j][k] << endl << "  "
                       << jacvalimhi_h[i][j][k] << "  "
                       << jacvalimlo_h[i][j][k] << endl;
                  cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                       << jacvalrehi_d[i][j][k] << "  "
                       << jacvalrelo_d[i][j][k] << endl << "  "
                       << jacvalimhi_d[i][j][k] << "  "
                       << jacvalimlo_d[i][j][k] << endl;
               }
               errsum += abs(jacvalrehi_h[i][j][k] - jacvalrehi_d[i][j][k])
                       + abs(jacvalrelo_h[i][j][k] - jacvalrelo_d[i][j][k])
                       + abs(jacvalimhi_h[i][j][k] - jacvalimhi_d[i][j][k])
                       + abs(jacvalimlo_h[i][j][k] - jacvalimlo_d[i][j][k]);
            }
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "rhs_h[" << i << "][" << j << "] : "
                    << rhsrehi_h[i][j] << "  "
                    << rhsrelo_h[i][j] << endl << "  "
                    << rhsimhi_h[i][j] << "  "
                    << rhsimlo_h[i][j] << endl;
               cout << "rhs_d[" << i << "][" << j << "] : "
                    << rhsrehi_d[i][j] << "  "
                    << rhsrelo_d[i][j] << endl << "  "
                    << rhsimhi_d[i][j] << "  "
                    << rhsimlo_d[i][j] << endl;
            }
            errsum += abs(rhsrehi_h[i][j] - rhsrehi_d[i][j])
                    + abs(rhsrelo_h[i][j] - rhsrelo_d[i][j])
                    + abs(rhsimhi_h[i][j] - rhsimhi_d[i][j])
                    + abs(rhsimlo_h[i][j] - rhsimlo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhsrehi_h[i][j] = rhsrehi_h[i][j];
         urhsimhi_h[i][j] = rhsimhi_h[i][j];
         urhsrelo_h[i][j] = rhsrelo_h[i][j];
         urhsimlo_h[i][j] = rhsimlo_h[i][j];
         urhsrehi_d[i][j] = rhsrehi_d[i][j];
         urhsimhi_d[i][j] = rhsimhi_d[i][j];
         urhsrelo_d[i][j] = rhsrelo_d[i][j];
         urhsimlo_d[i][j] = rhsimlo_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solrehi_h[i][j] = 0.0; solimhi_h[i][j] = 0.0;
         solrelo_h[i][j] = 0.0; solimlo_h[i][j] = 0.0;
         solrehi_d[i][j] = 0.0; solimhi_d[i][j] = 0.0;
         solrelo_d[i][j] = 0.0; solimlo_d[i][j] = 0.0;
      }

   if((mode == 1) || (mode == 2))
   {
      CPU_cmplx2_qrbs_solve
         (dim,degp1,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
          urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,
          workmatrehi,workmatrelo,workmatimhi,workmatimlo,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,vrblvl);
 
      if(vrblvl > 0)
      {
         CPU_cmplx2_linear_residue
            (dim,degp1,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,solrehi_h,
             solrelo_h,solimhi_h,solimlo_h,
             resvecrehi,resvecrelo,resvecimhi,resvecimlo,
             resmaxhi,resmaxlo,vrblvl);
         cout << "maximum residual : " << *resmaxhi << endl;
      }
      cmplx2_update_series
         (dim,degp1,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      GPU_cmplx2_bals_solve
         (dim,degp1,szt,nbt,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
          Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,vrblvl);

      if(vrblvl > 0)
      {
         CPU_cmplx2_linear_residue
            (dim,degp1,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
             solrehi_d,solrelo_d,solimhi_d,solimlo_d,
             resvecrehi,resvecrelo,resvecimhi,resvecimlo,
             resmaxhi,resmaxlo,vrblvl);
         cout << "maximum residual : " << *resmaxhi << endl;
      }
      cmplx2_update_series
         (dim,degp1,inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "Q_h[" << i << "][" << j << "] : "
                    << Qrehi_h[i][j] << "  " << Qrelo_h[i][j] << endl << "  "
                    << Qimhi_h[i][j] << "  " << Qimlo_h[i][j] << endl;
               cout << "Q_d[" << i << "][" << j << "] : "
                    << Qrehi_d[i][j] << "  " << Qrelo_d[i][j] << endl << "  "
                    << Qimhi_d[i][j] << "  " << Qimlo_d[i][j] << endl;
            }
            errsum += abs(Qrehi_h[i][j] - Qrehi_d[i][j])
                    + abs(Qrelo_h[i][j] - Qrelo_d[i][j])
                    + abs(Qimhi_h[i][j] - Qimhi_d[i][j])
                    + abs(Qimlo_h[i][j] - Qimlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "R_h[" << i << "][" << j << "] : "
                    << Rrehi_h[i][j] << "  " << Rrelo_h[i][j] << endl << "  "
                    << Rimhi_h[i][j] << "  " << Rimlo_h[i][j] << endl;
               cout << "R_d[" << i << "][" << j << "] : "
                    << Rrehi_d[i][j] << "  " << Rrelo_d[i][j] << endl << "  "
                    << Rimhi_d[i][j] << "  " << Rimlo_d[i][j] << endl;
            }
            errsum += abs(Rrehi_h[i][j] - Rrehi_d[i][j])
                    + abs(Rrelo_h[i][j] - Rrelo_d[i][j])
                    + abs(Rimhi_h[i][j] - Rimhi_d[i][j])
                    + abs(Rimlo_h[i][j] - Rimlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "urhs_h[" << i << "][" << j << "] : "
                    << urhsrehi_h[i][j] << "  "
                    << urhsrelo_h[i][j] << endl << "  "
                    << urhsimhi_h[i][j] << "  "
                    << urhsimlo_h[i][j] << endl;
               cout << "urhs_d[" << i << "][" << j << "] : "
                    << urhsrehi_d[i][j] << "  "
                    << urhsrelo_d[i][j] << endl << "  "
                    << urhsimhi_d[i][j] << "  "
                    << urhsimlo_d[i][j] << endl;
            }
            errsum += abs(urhsrehi_h[i][j] - urhsrehi_d[i][j])
                    + abs(urhsrelo_h[i][j] - urhsrelo_d[i][j])
                    + abs(urhsimhi_h[i][j] - urhsimhi_d[i][j])
                    + abs(urhsimlo_h[i][j] - urhsimlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "sol_h[" << i << "][" << j << "] : "
                    << solrehi_h[i][j] << "  "
                    << solrelo_h[i][j] << endl << "  "
                    << solimhi_h[i][j] << "  "
                    << solimlo_h[i][j] << endl;
               cout << "sol_d[" << i << "][" << j << "] : "
                    << solrehi_d[i][j] << "  "
                    << solrelo_d[i][j] << endl << "  "
                    << solimhi_d[i][j] << "  "
                    << solimlo_d[i][j] << endl;
            }
            errsum += abs(solrehi_h[i][j] - solrehi_d[i][j])
                    + abs(solrelo_h[i][j] - solrelo_d[i][j])
                    + abs(solimhi_h[i][j] - solimhi_d[i][j])
                    + abs(solimlo_h[i][j] - solimlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "input_h[" << i << "][" << j << "] : "
                    << inputrehi_h[i][j] << "  "
                    << inputrelo_h[i][j] << endl << "  "
                    << inputimhi_h[i][j] << "  "
                    << inputimlo_h[i][j] << endl;
               cout << "input_d[" << i << "][" << j << "] : "
                    << inputrehi_d[i][j] << "  "
                    << inputrelo_d[i][j] << endl << "  "
                    << inputimhi_d[i][j] << "  "
                    << inputimlo_d[i][j] << endl;
            }
            errsum += abs(inputrehi_h[i][j] - inputrehi_d[i][j])
                    + abs(inputrelo_h[i][j] - inputrelo_d[i][j])
                    + abs(inputimhi_h[i][j] - inputimhi_d[i][j])
                    + abs(inputimlo_h[i][j] - inputimlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl2_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhi_h = new double*[dim];
   double **inputlo_h = new double*[dim];
   double **inputhi_d = new double*[dim];
   double **inputlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhi_h[i] = new double[degp1];
      inputlo_h[i] = new double[degp1];
      inputhi_d[i] = new double[degp1];
      inputlo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acchi = new double[degp1]; // accumulated power series
   double *acclo = new double[degp1];
   double **cffhi = new double*[dim]; // the coefficients of monomials
   double **cfflo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhi[i] = new double[degp1];
      cfflo[i] = new double[degp1];
   }
   double ***outputhi_h = new double**[dim];
   double ***outputlo_h = new double**[dim];
   double ***outputhi_d = new double**[dim];
   double ***outputlo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhi_h[i] = new double*[dim+1];
      outputlo_h[i] = new double*[dim+1];
      outputhi_d[i] = new double*[dim+1];
      outputlo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhi_h[i][j] = new double[degp1];
         outputlo_h[i][j] = new double[degp1];
         outputhi_d[i][j] = new double[degp1];
         outputlo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhi_h = new double*[dim];
   double **funvallo_h = new double*[dim];
   double **funvalhi_d = new double*[dim];
   double **funvallo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhi_h[i] = new double[degp1];
      funvallo_h[i] = new double[degp1];
      funvalhi_d[i] = new double[degp1];
      funvallo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h = new double**[degp1];
   double ***jacvallo_h = new double**[degp1];
   double ***jacvalhi_d = new double**[degp1];
   double ***jacvallo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhi_h[i] = new double*[dim];
      jacvallo_h[i] = new double*[dim];
      jacvalhi_d[i] = new double*[dim];
      jacvallo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhi_h[i][j] = new double[dim];
         jacvallo_h[i][j] = new double[dim];
         jacvalhi_d[i][j] = new double[dim];
         jacvallo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **solhi_h = new double*[degp1];
   double **sollo_h = new double*[degp1];
   double **solhi_d = new double*[degp1];
   double **sollo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhi_h[i] = new double[dim];
      sollo_h[i] = new double[dim];
      solhi_d[i] = new double[dim];
      sollo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshi_h = new double*[degp1];
   double **rhslo_h = new double*[degp1];
   double **rhshi_d = new double*[degp1];
   double **rhslo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshi_h[i] = new double[dim];
      rhslo_h[i] = new double[dim];
      rhshi_d[i] = new double[dim];
      rhslo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshi_h = new double*[degp1];
   double **urhslo_h = new double*[degp1];
   double **urhshi_d = new double*[degp1];
   double **urhslo_d = new double*[degp1];
   for(int i=0; i<degp1; i++)
   {
      urhshi_h[i] = new double[dim];
      urhslo_h[i] = new double[dim];
      urhshi_d[i] = new double[dim];
      urhslo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathi = new double*[dim];
   double **workmatlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathi[i] = new double[dim];
      workmatlo[i] = new double[dim];
   }
   // int *ipvt = new int[dim];
   double *workvechi = new double[dim];
   double *workveclo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshi = new double*[degp1];
   double **workrhslo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshi[i] = new double[dim];
      workrhslo[i] = new double[dim];
   }
   double **resvechi = new double*[degp1];
   double **resveclo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechi[i] = new double[dim];
      resveclo[i] = new double[dim];
   }
   double resmaxhi,resmaxlo;
   double **Qhi_h = new double*[dim];
   double **Qlo_h = new double*[dim];
   double **Qhi_d = new double*[dim];
   double **Qlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qhi_h[i] = new double[dim];
      Qlo_h[i] = new double[dim];
      Qhi_d[i] = new double[dim];
      Qlo_d[i] = new double[dim];
   }
   double **Rhi_h = new double*[dim];
   double **Rlo_h = new double*[dim];
   double **Rhi_d = new double*[dim];
   double **Rlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhi_h[i] = new double[dim];
      Rlo_h[i] = new double[dim];
      Rhi_d[i] = new double[dim];
      Rlo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl2_unit_series_vector(dim,deg,inputhi_h,inputlo_h);
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhi_d[i][j] = inputhi_h[i][j];
         inputlo_d[i][j] = inputlo_h[i][j];
      }

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputhi_h[i][0] << "  "
                            << inputlo_h[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;
/*
      dbl2_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,cffhi,cfflo,acchi,acclo,
          inputhi,inputlo,outputhi,outputlo,funvalhi,funvallo,
          jacvalhi,jacvallo,rhshi,rhslo,solhi,sollo,workmathi,workmatlo,
          workvechi,workveclo,workrhshi,workrhslo,resvechi,resveclo,
          &resmaxhi,&resmaxlo,ipvt,vrblvl);
 */
      dbl2_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,cffhi,cfflo,acchi,acclo,
          inputhi_h,inputlo_h,inputhi_d,inputlo_d,
          outputhi_h,outputlo_h,outputhi_d,outputlo_d,
          funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,urhshi_h,urhslo_h,urhshi_d,urhslo_d,
          solhi_h,sollo_h,solhi_d,sollo_d,
          Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
          workmathi,workmatlo,workvechi,workveclo,
          resvechi,resveclo,&resmaxhi,&resmaxlo,vrblvl,mode);
   }
   return 0;
}

int test_dbl2_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputrehi_h = new double*[dim];
   double **inputrelo_h = new double*[dim];
   double **inputimhi_h = new double*[dim];
   double **inputimlo_h = new double*[dim];
   double **inputrehi_d = new double*[dim];
   double **inputrelo_d = new double*[dim];
   double **inputimhi_d = new double*[dim];
   double **inputimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputrehi_h[i] = new double[degp1];
       inputrelo_h[i] = new double[degp1];
       inputimhi_h[i] = new double[degp1];
       inputimlo_h[i] = new double[degp1];
       inputrehi_d[i] = new double[degp1];
       inputrelo_d[i] = new double[degp1];
       inputimhi_d[i] = new double[degp1];
       inputimlo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *accrehi = new double[degp1]; // accumulated power series
   double *accrelo = new double[degp1];
   double *accimhi = new double[degp1];
   double *accimlo = new double[degp1];
   double **cffrehi = new double*[dim]; // the coefficients of monomials
   double **cffrelo = new double*[dim];
   double **cffimhi = new double*[dim]; 
   double **cffimlo = new double*[dim]; 
   for(int i=0; i<dim; i++)
   {
      cffrehi[i] = new double[degp1];
      cffrelo[i] = new double[degp1];
      cffimhi[i] = new double[degp1];
      cffimlo[i] = new double[degp1];
   }
   double ***outputrehi_h = new double**[dim];
   double ***outputrelo_h = new double**[dim];
   double ***outputimhi_h = new double**[dim];
   double ***outputimlo_h = new double**[dim];
   double ***outputrehi_d = new double**[dim];
   double ***outputrelo_d = new double**[dim];
   double ***outputimhi_d = new double**[dim];
   double ***outputimlo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehi_h[i] = new double*[dim+1];
      outputrelo_h[i] = new double*[dim+1];
      outputimhi_h[i] = new double*[dim+1];
      outputimlo_h[i] = new double*[dim+1];
      outputrehi_d[i] = new double*[dim+1];
      outputrelo_d[i] = new double*[dim+1];
      outputimhi_d[i] = new double*[dim+1];
      outputimlo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehi_h[i][j] = new double[degp1];
         outputrelo_h[i][j] = new double[degp1];
         outputimhi_h[i][j] = new double[degp1];
         outputimlo_h[i][j] = new double[degp1];
         outputrehi_d[i][j] = new double[degp1];
         outputrelo_d[i][j] = new double[degp1];
         outputimhi_d[i][j] = new double[degp1];
         outputimlo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalrehi_h = new double*[dim];
   double **funvalrelo_h = new double*[dim];
   double **funvalimhi_h = new double*[dim];
   double **funvalimlo_h = new double*[dim];
   double **funvalrehi_d = new double*[dim];
   double **funvalrelo_d = new double*[dim];
   double **funvalimhi_d = new double*[dim];
   double **funvalimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalrehi_h[i] = new double[degp1];
      funvalrelo_h[i] = new double[degp1];
      funvalimhi_h[i] = new double[degp1];
      funvalimlo_h[i] = new double[degp1];
      funvalrehi_d[i] = new double[degp1];
      funvalrelo_d[i] = new double[degp1];
      funvalimhi_d[i] = new double[degp1];
      funvalimlo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalrehi_h = new double**[degp1];
   double ***jacvalrelo_h = new double**[degp1];
   double ***jacvalimhi_h = new double**[degp1];
   double ***jacvalimlo_h = new double**[degp1];
   double ***jacvalrehi_d = new double**[degp1];
   double ***jacvalrelo_d = new double**[degp1];
   double ***jacvalimhi_d = new double**[degp1];
   double ***jacvalimlo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalrehi_h[i] = new double*[dim];
      jacvalrelo_h[i] = new double*[dim];
      jacvalimhi_h[i] = new double*[dim];
      jacvalimlo_h[i] = new double*[dim];
      jacvalrehi_d[i] = new double*[dim];
      jacvalrelo_d[i] = new double*[dim];
      jacvalimhi_d[i] = new double*[dim];
      jacvalimlo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalrehi_h[i][j] = new double[dim];
         jacvalrelo_h[i][j] = new double[dim];
         jacvalimhi_h[i][j] = new double[dim];
         jacvalimlo_h[i][j] = new double[dim];
         jacvalrehi_d[i][j] = new double[dim];
         jacvalrelo_d[i][j] = new double[dim];
         jacvalimhi_d[i][j] = new double[dim];
         jacvalimlo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solrehi_h = new double*[degp1];
   double **solrelo_h = new double*[degp1];
   double **solimhi_h = new double*[degp1];
   double **solimlo_h = new double*[degp1];
   double **solrehi_d = new double*[degp1];
   double **solrelo_d = new double*[degp1];
   double **solimhi_d = new double*[degp1];
   double **solimlo_d = new double*[degp1];

   for(int i=0; i<degp1; i++) 
   {
      solrehi_h[i] = new double[dim];
      solrelo_h[i] = new double[dim];
      solimhi_h[i] = new double[dim];
      solimlo_h[i] = new double[dim];
      solrehi_d[i] = new double[dim];
      solrelo_d[i] = new double[dim];
      solimhi_d[i] = new double[dim];
      solimlo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsrehi_h = new double*[degp1];
   double **rhsrelo_h = new double*[degp1];
   double **rhsimhi_h = new double*[degp1];
   double **rhsimlo_h = new double*[degp1];
   double **rhsrehi_d = new double*[degp1];
   double **rhsrelo_d = new double*[degp1];
   double **rhsimhi_d = new double*[degp1];
   double **rhsimlo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsrehi_h[i] = new double[dim];
      rhsrelo_h[i] = new double[dim];
      rhsimhi_h[i] = new double[dim];
      rhsimlo_h[i] = new double[dim];
      rhsrehi_d[i] = new double[dim];
      rhsrelo_d[i] = new double[dim];
      rhsimhi_d[i] = new double[dim];
      rhsimlo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmatrehi = new double*[dim];
   double **workmatrelo = new double*[dim];
   double **workmatimhi = new double*[dim];
   double **workmatimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmatrehi[i] = new double[dim];
      workmatrelo[i] = new double[dim];
      workmatimhi[i] = new double[dim];
      workmatimlo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvecrehi = new double[dim];
   double *workvecrelo = new double[dim];
   double *workvecimhi = new double[dim];
   double *workvecimlo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **urhsrehi_h = new double*[degp1];
   double **urhsrelo_h = new double*[degp1];
   double **urhsimhi_h = new double*[degp1];
   double **urhsimlo_h = new double*[degp1];
   double **urhsrehi_d = new double*[degp1];
   double **urhsrelo_d = new double*[degp1];
   double **urhsimhi_d = new double*[degp1];
   double **urhsimlo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhsrehi_h[i] = new double[dim];
      urhsrelo_h[i] = new double[dim];
      urhsimhi_h[i] = new double[dim];
      urhsimlo_h[i] = new double[dim];
      urhsrehi_d[i] = new double[dim];
      urhsrelo_d[i] = new double[dim];
      urhsimhi_d[i] = new double[dim];
      urhsimlo_d[i] = new double[dim];
   }
   double **resvecrehi = new double*[degp1];
   double **resvecrelo = new double*[degp1];
   double **resvecimhi = new double*[degp1];
   double **resvecimlo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecrehi[i] = new double[dim];
      resvecrelo[i] = new double[dim];
      resvecimhi[i] = new double[dim];
      resvecimlo[i] = new double[dim];
   }
   double resmaxhi,resmaxlo;
   double **Qrehi_h = new double*[dim];
   double **Qrelo_h = new double*[dim];
   double **Qimhi_h = new double*[dim];
   double **Qimlo_h = new double*[dim];
   double **Qrehi_d = new double*[dim];
   double **Qrelo_d = new double*[dim];
   double **Qimhi_d = new double*[dim];
   double **Qimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qrehi_h[i] = new double[dim];
      Qrelo_h[i] = new double[dim];
      Qimhi_h[i] = new double[dim];
      Qimlo_h[i] = new double[dim];
      Qrehi_d[i] = new double[dim];
      Qrelo_d[i] = new double[dim];
      Qimhi_d[i] = new double[dim];
      Qimlo_d[i] = new double[dim];
   }
   double **Rrehi_h = new double*[dim];
   double **Rrelo_h = new double*[dim];
   double **Rimhi_h = new double*[dim];
   double **Rimlo_h = new double*[dim];
   double **Rrehi_d = new double*[dim];
   double **Rrelo_d = new double*[dim];
   double **Rimhi_d = new double*[dim];
   double **Rimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rrehi_h[i] = new double[dim];
      Rrelo_h[i] = new double[dim];
      Rimhi_h[i] = new double[dim];
      Rimlo_h[i] = new double[dim];
      Rrehi_d[i] = new double[dim];
      Rrelo_d[i] = new double[dim];
      Rimhi_d[i] = new double[dim];
      Rimlo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   cmplx2_unit_series_vector
      (dim,deg,inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehi_d[i][j] = inputrehi_h[i][j];
         inputrelo_d[i][j] = inputrelo_h[i][j];
         inputimhi_d[i][j] = inputimhi_h[i][j];
         inputimlo_d[i][j] = inputimlo_h[i][j];
      }

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << inputrehi_h[i][0] << "  "
              << inputrelo_h[i][0] << endl << "  "
              << inputimhi_h[i][0] << "  "
              << inputimlo_h[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;

      cmplx2_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffrehi,cffrelo,cffimhi,cffimlo,accrehi,accrelo,accimhi,accimlo,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
          inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
          funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
          funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
          rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
          rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
          urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
          urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
          Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          workmatrehi,workmatrelo,workmatimhi,workmatimlo,
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,
          resvecrehi,resvecrelo,resvecimhi,resvecimlo,
          &resmaxhi,&resmaxlo,vrblvl,mode);
   }
   return 0;
}
