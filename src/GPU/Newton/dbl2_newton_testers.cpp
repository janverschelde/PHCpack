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
      cffhi[i][0] = 1.0;
      cfflo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
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
             cout << "output_h[" << k << "][" << i << "][" << j << "] : "
                  << outputhi_h[k][i][j] << "  "
                  << outputlo_h[k][i][j] << endl;
             cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                  << outputhi_d[k][i][j] << "  "
                  << outputlo_d[k][i][j] << endl;
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
            cout << "funval_h[" << i << "][" << j << "] : "
                 << funvalhi_h[i][j] << "  " << funvallo_h[i][j] << endl;
            cout << "funval_d[" << i << "][" << j << "] : "
                 << funvalhi_d[i][j] << "  " << funvallo_d[i][j] << endl;
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
               cout << "jacval_h[" << i << "][" << j << "][" << k << "] : "
                    << jacvalhi_h[i][j][k] << "  "
                    << jacvallo_h[i][j][k] << endl;
               cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                    << jacvalhi_d[i][j][k] << "  "
                    << jacvallo_d[i][j][k] << endl;
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
            cout << "rhs_h[" << i << "][" << j << "] : "
                 << rhshi_h[i][j] << "  " << rhslo_h[i][j] << endl;
            cout << "rhs_d[" << i << "][" << j << "] : "
                 << rhshi_d[i][j] << "  " << rhslo_d[i][j] << endl;
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
   
         cout << "maximum residual : "
              << *resmaxhi << "  " << *resmaxlo << endl;
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
   
         cout << "maximum residual : "
              << *resmaxhi << "  " << *resmaxlo << endl;
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
             cout << "Q_h[" << i << "][" << j << "] : "
                  << Qhi_h[i][j] << "  " << Qlo_h[i][j] << endl;
             cout << "Q_d[" << i << "][" << j << "] : "
                  << Qhi_d[i][j] << "  " << Qlo_d[i][j] << endl;
             errsum += abs(Qhi_h[i][j] - Qhi_d[i][j])
                     + abs(Qlo_h[i][j] - Qlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "R_h[" << i << "][" << j << "] : "
                  << Rhi_h[i][j] << "  " << Rlo_h[i][j] << endl;
             cout << "R_d[" << i << "][" << j << "] : "
                  << Rhi_d[i][j] << "  " << Rlo_d[i][j] << endl;
             errsum += abs(Rhi_h[i][j] - Rhi_d[i][j])
                     + abs(Rlo_h[i][j] - Rlo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "urhs_h[" << i << "][" << j << "] : "
                  << urhshi_h[i][j] << "  " << urhslo_h[i][j] << endl;
             cout << "urhs_d[" << i << "][" << j << "] : "
                  << urhshi_d[i][j] << "  " << urhslo_d[i][j] << endl;
             errsum += abs(urhshi_h[i][j] - urhshi_d[i][j])
                     + abs(urhslo_h[i][j] - urhslo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "sol_h[" << i << "][" << j << "] : "
                  << solhi_h[i][j] << "  " << sollo_h[i][j] << endl;
             cout << "sol_d[" << i << "][" << j << "] : "
                  << solhi_d[i][j] << "  " << sollo_d[i][j] << endl;
             errsum += abs(solhi_h[i][j] - solhi_d[i][j])
                     + abs(sollo_h[i][j] - sollo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
             cout << "input_h[" << i << "][" << j << "] : "
                  << inputhi_h[i][j] << "  "
                  << inputlo_h[i][j] << endl;
             cout << "input_d[" << i << "][" << j << "] : "
                  << inputhi_d[i][j] << "  "
                  << inputlo_d[i][j] << endl;
             errsum += abs(inputhi_h[i][j] - inputhi_d[i][j])
                     + abs(inputlo_h[i][j] - inputlo_d[i][j]);
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
      cout << "step " << step << " ..." << endl;
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
