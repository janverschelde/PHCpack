// The file dbl_newton_testers.cpp defines the functions with prototypes in
// the file dbl_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "dbl_factorizations.h"
#include "dbl_systems_host.h"
#include "dbl_systems_kernels.h"
#include "dbl_bals_host.h"
#include "dbl_bals_kernels.h"

using namespace std;

void prompt_newton_setup
 ( int *seed, int *szt, int*nbt, int *dim, int *deg, int *size, int *posvals,
   int *vrblvl, int *mode, int *nbritr, int *nbsteps )
{
   cout << "-> give the seed (0 for time) : "; cin >> *seed;

   prompt_dimensions(dim,deg,size,posvals,vrblvl,nbritr,nbsteps);

   cout << "-> enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> *mode;

   if(*mode != 1)
   {
      cout << "-> give the number of tiles : "; cin >> *nbt;
      cout << "-> give the size of each tile : "; cin >> *szt;
      int p = (*szt)*(*nbt);

      while(p != *dim)
      {
          cout << "Dimension = " << *dim << " != " << *szt << " * " << *nbt
               << ", retry." << endl;
          cout << "-> give the size of each tile : "; cin >> *szt;
          cout << "-> give the number of tiles : "; cin >> *nbt;
          p = (*szt)*(*nbt);
      }
   }
}

void dbl_unit_series_vector ( int dim, int deg, double **cff )
{
   for(int i=0; i<dim; i++)
   {
      cff[i][0] = 1.0;
      for(int j=1; j<=deg; j++) cff[i][j] = 0.0;
   }
}

void dbl_update_series
 ( int dim, int degp1, double **x, double **dx, int vrblvl )
{
   if(vrblvl > 0)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++) cout << x[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++) x[i][j] = x[i][j] + dx[j][i];

   if(vrblvl > 0)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++) cout << x[i][j] << endl;
      }
   }
}

void dbl_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cff, double *acc, double **input, double ***output,
   double **funval, double ***jacval, double **rhs, double **sol,
   double **workmat, double *workvec, double **workrhs, double **resvec,
   double *resmax, int *ipvt, int vrblvl )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   for(int i=0; i<dim; i++)
   {
      cff[i][0] = 1.0;
      for(int j=1; j<degp1; j++) cff[i][j] = 0.0;
   }
   CPU_dbl_evaluate_monomials
      (dim,deg,nvr,idx,exp,nbrfac,expfac,cff,acc,input,output,vrblvl);

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++) jacval[i][j][k] = 0.0;

   dbl_linearize_evaldiff_output
      (dim,degp1,nvr,idx,output,funval,rhs,jacval,vrblvl);

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++) workrhs[i][j] = rhs[i][j];
   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++) sol[i][j] = 0.0;

   CPU_dbl_lusb_solve
      (dim,degp1,jacval,workrhs,sol,workmat,workvec,ipvt,0); // vrblvl);

   CPU_dbl_linear_residue(dim,degp1,jacval,rhs,sol,resvec,resmax,vrblvl);

   if(vrblvl > 0) cout << "maximum residual : " << *resmax << endl;

   dbl_update_series(dim,degp1,input,sol,vrblvl);
}

void dbl_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cff, double *acc,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **workmat, double *workvec, double **resvec, double *resmax,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   for(int i=0; i<dim; i++)
   {
      cff[i][0] = 1.0;
      for(int j=1; j<degp1; j++) cff[i][j] = 0.0;
      // for(int j=0; j<degp1; j++) input_d[i][j] = input_h[i][j];
   }
   if((mode == 1) || (mode == 2))
      CPU_dbl_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,cff,acc,input_h,output_h,vrblvl);
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)  // reset the coefficients
      {
         cff[i][0] = 1.0;
         for(int j=1; j<degp1; j++) cff[i][j] = 0.0;
      }
      GPU_dbl_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,cff,acc,
          input_d,output_d,vrblvl);
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
                  << output_h[k][i][j] << endl;
             cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                  << output_d[k][i][j] << endl;
             errsum += abs(output_h[k][i][j] - output_d[k][i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacval_h[i][j][k] = 0.0;
            jacval_d[i][j][k] = 0.0;
         }

   if((mode == 1) || (mode == 2))
      dbl_linearize_evaldiff_output
         (dim,degp1,nvr,idx,output_h,funval_h,rhs_h,jacval_h,vrblvl);
   if((mode == 0) || (mode == 2))
      dbl_linearize_evaldiff_output
         (dim,degp1,nvr,idx,output_d,funval_d,rhs_d,jacval_d,vrblvl);

   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
         {
            cout << "funval_h[" << i << "][" << j << "] : "
                 << funval_h[i][j] << endl;
            cout << "funval_d[" << i << "][" << j << "] : "
                 << funval_d[i][j] << endl;
            errsum += abs(funval_h[i][j] - funval_d[i][j]);
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
                    << jacval_h[i][j][k] << endl;
               cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                    << jacval_d[i][j][k] << endl;
               errsum += abs(jacval_h[i][j][k] - jacval_d[i][j][k]);
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
                 << rhs_h[i][j] << endl;
            cout << "rhs_d[" << i << "][" << j << "] : "
                 << rhs_d[i][j] << endl;
            errsum += abs(rhs_h[i][j] - rhs_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhs_h[i][j] = rhs_h[i][j];
         urhs_d[i][j] = rhs_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         sol_h[i][j] = 0.0;
         sol_d[i][j] = 0.0;
      }

   if((mode == 1) || (mode == 2))
   {
      CPU_dbl_qrbs_solve
         (dim,degp1,jacval_h,urhs_h,sol_h,workmat,Q_h,R_h,workvec,vrblvl);
 
      if(vrblvl > 0)
      {
         CPU_dbl_linear_residue
            (dim,degp1,jacval_h,rhs_h,sol_h,resvec,resmax,vrblvl);
         cout << "maximum residual : " << *resmax << endl;
      }
      dbl_update_series(dim,degp1,input_h,sol_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      GPU_dbl_bals_solve
         (dim,degp1,szt,nbt,jacval_d,Q_d,R_d,urhs_d,sol_d,vrblvl);

      if(vrblvl > 0)
      {
         CPU_dbl_linear_residue
            (dim,degp1,jacval_d,rhs_d,sol_d,resvec,resmax,vrblvl);
         cout << "maximum residual : " << *resmax << endl;
      }
      dbl_update_series(dim,degp1,input_d,sol_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "Q_h[" << i << "][" << j << "] : "
                  << Q_h[i][j] << endl;
             cout << "Q_d[" << i << "][" << j << "] : "
                  << Q_d[i][j] << endl;
             errsum += abs(Q_h[i][j] - Q_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "R_h[" << i << "][" << j << "] : "
                  << R_h[i][j] << endl;
             cout << "R_d[" << i << "][" << j << "] : "
                  << R_d[i][j] << endl;
             errsum += abs(R_h[i][j] - R_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "urhs_h[" << i << "][" << j << "] : "
                  << urhs_h[i][j] << endl;
             cout << "urhs_d[" << i << "][" << j << "] : "
                  << urhs_d[i][j] << endl;
             errsum += abs(urhs_h[i][j] - urhs_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "sol_h[" << i << "][" << j << "] : "
                  << sol_h[i][j] << endl;
             cout << "sol_d[" << i << "][" << j << "] : "
                  << sol_d[i][j] << endl;
             errsum += abs(sol_h[i][j] - sol_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
             cout << "input_h[" << i << "][" << j << "] : "
                  << input_h[i][j] << endl;
             cout << "input_d[" << i << "][" << j << "] : "
                  << input_d[i][j] << endl;
             errsum += abs(input_h[i][j] - input_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}
