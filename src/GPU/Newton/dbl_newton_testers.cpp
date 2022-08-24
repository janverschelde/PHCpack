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
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"
#include "dbl_factorizations.h"
#include "dbl_bals_host.h"
#include "dbl_bals_kernels.h"

using namespace std;

void prompt_newton_setup
 ( int *seed, int *szt, int*nbt, int *dim, int *deg, int *size, int *posvals,
   int *vrblvl, int *mode, int *nbritr, int *nbsteps )
{
   cout << "-> give the seed (0 for time) : "; cin >> *seed;

   prompt_dimensions(dim,deg,size,posvals,vrblvl,nbritr,nbsteps);

   cout << "-> give the size of each tile : "; cin >> *szt;
   cout << "-> give the number of tiles : "; cin >> *nbt;

   cout << "-> enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> *mode;

   if(*mode != 1)
   {
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

void dbl_evaluate_monomials
 ( int dim, int deg, int *nvr, int **idx, int **exp, int *nbrfac,
   int **expfac, double **cff, double *acc, double **input,
   double ***output, int vrblvl )
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
                  CPU_dbl_product(deg,input[idxvar],cff[i],acc);
                  for(int L=0; L<=deg; L++) cff[i][L] = acc[L];
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
         for(int j=0; j<=deg; j++) cout << cff[i][j] << endl;
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
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl_evaldiff(dim,nvr[i],deg,idx[i],cff[i],input,output[i]);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << output[i][dim][j] << endl;
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

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
                  output[i][idxvar][k] = factor*output[i][idxvar][k];
            }
         }
      }
   }
}

void dbl_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx, double ***output,
   double **funval, double **rhs, double ***jacval, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++) funval[i][j] = output[i][dim][j];

   if(vrblvl > 0)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++) cout << i << " : " << funval[i][0] << endl;
   }
   // Linearize the function values in the rhs and swap sign,
   // but keep in mind that the right hand side is 1 - t,
   // so we subtract 1 and add t to the rhs.
   for(int j=0; j<dim; j++) rhs[0][j] = -(funval[j][0] - 1.0);
   if(degp1 > 1)
   {
      for(int j=0; j<dim; j++) rhs[1][j] = -(funval[j][1] + 1.0);
      for(int i=2; i<degp1; i++)
         for(int j=0; j<dim; j++) rhs[i][j] = -funval[j][i];
   }
   if(vrblvl > 0)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << rhs[i][j] << endl;
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
            jacval[j][i][idxval] = output[i][idxval][j];
      }
   }
   if(vrblvl > 0)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacval[0][i][j] << endl;
      }
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
   dbl_evaluate_monomials
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
   double **input_h, double **input_d, double ***output,
   double **funval, double ***jacval, double **rhs,
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
      for(int j=0; j<degp1; j++) input_d[i][j] = input_h[i][j];
   }
   if((mode == 1) || (mode == 2))
      dbl_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,cff,acc,input_h,output,vrblvl);
   if(mode == 0)
      dbl_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,cff,acc,input_d,output,vrblvl);

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++) jacval[i][j][k] = 0.0;

   dbl_linearize_evaldiff_output
      (dim,degp1,nvr,idx,output,funval,rhs,jacval,vrblvl);

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhs_h[i][j] = rhs[i][j];
         urhs_d[i][j] = rhs[i][j];
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
         (dim,degp1,jacval,urhs_h,sol_h,workmat,Q_h,R_h,workvec,vrblvl);

      CPU_dbl_linear_residue(dim,degp1,jacval,rhs,sol_h,resvec,resmax,vrblvl);
      if(vrblvl > 0) cout << "maximum residual : " << *resmax << endl;
      dbl_update_series(dim,degp1,input_h,sol_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      GPU_dbl_bals_solve
         (dim,degp1,szt,nbt,jacval,Q_d,R_d,urhs_d,sol_d,vrblvl);

      CPU_dbl_linear_residue(dim,degp1,jacval,rhs,sol_d,resvec,resmax,vrblvl);
      if(vrblvl > 0) cout << "maximum residual : " << *resmax << endl;
      dbl_update_series(dim,degp1,input_d,sol_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i< degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "urhs_h[" << i << "][" << j << "] : "
                  << urhs_h[i][j] << endl;
             cout << "urhs_d[" << i << "][" << j << "] : "
                  << urhs_d[i][j] << endl;
         }

      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i< degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "sol_h[" << i << "][" << j << "] : "
                  << sol_h[i][j] << endl;
             cout << "sol_d[" << i << "][" << j << "] : "
                  << sol_d[i][j] << endl;
         }

      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
             cout << "input_h[" << i << "][" << j << "] : "
                  << input_h[i][j] << endl;
             cout << "input_d[" << i << "][" << j << "] : "
                  << input_d[i][j] << endl;
         }
   }
}
