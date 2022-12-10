// The file dbl_newton_method.cpp defines the functions with prototypes in
// the file dbl_newton_method.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "random_numbers.h"
#include "random_monomials.h"
#include "dbl_factorizations.h"
#include "dbl_monomial_systems.h"
#include "unimodular_matrices.h"
#include "dbl_systems_host.h"
#include "dbl_systems_kernels.h"
#include "dbl_bals_host.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_kernels.h"
#include "dbl_newton_testers.h"

using namespace std;

void dbl_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mb, double dpr, double ***cff, double **acc,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   bool *noqr_h, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_dbl_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl_unit_series_vector(dim,deg,cff[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.
         CPU_dbl_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,cff[0],acc[0],
             input_h,output_h,vrblvl);
      }
      else
         CPU_dbl_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,cff,acc,input_h,funval_h,jacval_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl_unit_series_vector(dim,deg,cff[0]); // reset coefficients

         GPU_dbl_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,cff[0],acc[0],
             input_d,output_d,vrblvl);
      }
      else
         GPU_dbl_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cff,input_d,output_d,
             funval_d,jacval_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;

      errsum = dbl_error3sum
                  (dim,dim+1,degp1,output_h,output_d,"output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(16);
      cout << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl_define_rhs(dim,degp1,mb,funval_h,rhs_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++) jacval_h[i][j][k] = 0.0;

         dbl_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mb,dpr,output_h,funval_h,rhs_h,
             jacval_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl_define_rhs(dim,degp1,mb,funval_d,rhs_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++) jacval_d[i][j][k] = 0.0;

         dbl_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mb,dpr,output_d,funval_d,rhs_d,
             jacval_d,vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      errsum = dbl_error2sum(dim,degp1,funval_h,funval_d,"funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = dbl_error3sum
                  (degp1,dim,dim,jacval_h,jacval_d,"jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = dbl_error2sum(degp1,dim,rhs_h,rhs_d,"rhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++) urhs_h[i][j] = rhs_h[i][j];

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++) sol_h[i][j] = 0.0;

      if(vrblvl > 0)
         cout << "calling CPU_dbl_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_dbl_qrbs_solve
         (dim,degp1,oldtail,jacval_h,urhs_h,sol_h,Q_h,R_h,workvec,
          noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl_linear_residue ..." << endl;

         CPU_dbl_linear_residue
            (dim,degp1,*tailidx_h-1,jacval_h,rhs_h,sol_h,
             resvec,resmax,vrblvl);
         cout << "maximum residual : " << *resmax << endl;
      }
      dbl_update_series(dim,degp1,*tailidx_h-1,input_h,sol_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++) urhs_d[i][j] = rhs_d[i][j];

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++) sol_d[i][j] = 0.0;

      if(vrblvl > 0)
         cout << "calling GPU_dbl_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_dbl_bals_solve
         (dim,degp1,szt,nbt,oldtail,jacval_d,Q_d,R_d,urhs_d,sol_d,
          noqr_d,upidx_d,bsidx_d,&newtail,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl_linear_residue ..." << endl;

         double elapsedms;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_dbl_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,jacval_d,rhs_d,sol_d,
             resvec,resmax,&elapsedms,&addcnt,&mulcnt,vrblvl);
         cout << "maximum residual : " << *resmax << endl;
      }
      dbl_update_series(dim,degp1,*tailidx_d-1,input_d,sol_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      double errsum = dbl_error2sum(dim,dim,Q_h,Q_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = dbl_error2sum(dim,dim,R_h,R_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = dbl_error2sum(degp1,dim,urhs_h,urhs_d,"urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = dbl_error2sum(degp1,dim,sol_h,sol_d,"sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = dbl_error2sum(dim,degp1,input_h,input_d,"input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl_real_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **input_h = new double*[dim];
   double **input_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
       input_h[i] = new double[degp1];
       input_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **acc = new double*[dim+1]; // accumulate series in one column
   for(int i=0; i<=dim; i++) acc[i] = new double[degp1];
   double ***cff = new double**[nbrcol];
   for(int i=0; i<nbrcol; i++)
   {
      cff[i] = new double*[dim]; // the coefficients of monomials
      for(int j=0; j<dim; j++) cff[i][j] = new double[degp1];
   }
   double ***output_h;
   double ***output_d;

   if((mode == 1) || (mode == 2))
   {
      output_h = new double**[dim];
      for(int i=0; i<dim; i++)
      {
         output_h[i] = new double*[dim+1];
         for(int j=0; j<=dim; j++) output_h[i][j] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      output_d = new double**[dim];
      for(int i=0; i<dim; i++)
      {
         output_d[i] = new double*[dim+1];
         for(int j=0; j<=dim; j++) output_d[i][j] = new double[degp1];
      }
   }
   double **funval_h;  // function values on host
   double **funval_d;  // function values on device

   if((mode == 1) || (mode == 2))
   {
      funval_h = new double*[dim];
      for(int i=0; i<dim; i++) funval_h[i] = new double[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      funval_d = new double*[dim];
      for(int i=0; i<dim; i++) funval_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacval_h;
   double ***jacval_d;

   if((mode == 1) || (mode == 2))
   {
      jacval_h = new double**[degp1];
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacval_h[i] = new double*[dim];
         for(int j=0; j<dim; j++) jacval_h[i][j] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      jacval_d = new double**[degp1];
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacval_d[i] = new double*[dim];
         for(int j=0; j<dim; j++) jacval_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **sol_h;
   double **sol_d;

   if((mode == 1) || (mode == 2))
   {
      sol_h = new double*[degp1];
      for(int i=0; i<degp1; i++) sol_h[i] = new double[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      sol_d = new double*[degp1];
      for(int i=0; i<degp1; i++) sol_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhs_h = new double*[degp1];
   double **rhs_d = new double*[degp1];

   if((mode == 1) || (mode == 2))
   {
      rhs_h = new double*[degp1];
      for(int i=0; i<degp1; i++) rhs_h[i] = new double[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhs_d = new double*[degp1];
      for(int i=0; i<degp1; i++) rhs_d[i] = new double[dim];
   }
   double *workvec = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **urhs_h;
   double **urhs_d;

   if((mode == 1) || (mode == 2))
   {
      urhs_h = new double*[degp1];
      for(int i=0; i<degp1; i++) urhs_h[i] = new double[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      urhs_d = new double*[degp1];
      for(int i=0; i<degp1; i++) urhs_d[i] = new double[dim];
   }
   double **resvec = new double*[degp1];
   for(int i=0; i<degp1; i++) resvec[i] = new double[dim];
   double resmax;

   double **Q_h;
   double **Q_d;

   if((mode == 1) || (mode == 2))
   {
      Q_h = new double*[dim];
      for(int i=0; i<dim; i++) Q_h[i] = new double[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      Q_d = new double*[dim];
      for(int i=0; i<dim; i++) Q_d[i] = new double[dim];
   }
   double **R_h;
   double **R_d;

   if((mode == 1) || (mode == 2))
   {
      R_h = new double*[dim];
      for(int i=0; i<dim; i++) R_h[i] = new double[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      R_d = new double*[dim];
      for(int i=0; i<dim; i++) R_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   double **sol = new double*[dim];

   for(int i=0; i<dim; i++) sol[i] = new double[degp1];

   make_real_exponentials(dim,deg,sol);
   if(nbrcol != 1) // randomize the leading term
      for(int i=0; i<dim; i++)
         // sol[i][0] = sol[i][0] + random_double()/2.0;
         sol[i][0] = random_double();

   // compute the right hand sides via evaluation

   double **mbrhs = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhs[i] = new double[degp1];

      mbrhs[i][0] = 1.0;     // initialize product to one

      for(int k=1; k<degp1; k++) mbrhs[i][k] = 0.0;
   }
   if(nbrcol == 1)
      evaluate_real_monomials(dim,deg,rowsA,sol,mbrhs);
   else
   {
      evaluate_real_columns(dim,deg,nbrcol,nvr,idx,rowsA,sol,mbrhs,vrblvl);
      dbl_unit_series_vectors(nbrcol,dim,deg,cff);
   }
   
   double *start0 = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0[i] = sol[i][0];
   }
   real_start_series_vector(dim,deg,start0,input_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         // input_h[i][j] = sol[i][j]; // check if evaluation is done right
         input_d[i][j] = input_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << input_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool noqr_h = false;
   bool noqr_d = false;
   int tailidx_h = 1;
   int tailidx_d = 1;

   int wrkdeg = 0; // working degree of precision

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      dbl_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,mbrhs,dpr,cff,acc,
          input_h,input_d,output_h,output_d,funval_h,funval_d,
          jacval_h,jacval_d,rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
          Q_h,Q_d,R_h,R_d,workvec,resvec,&resmax,
          &noqr_h,&noqr_d,&upidx_h,&bsidx_h,&upidx_d,&bsidx_d,vrblvl,mode);

      if(vrblvl > 0)
         cout << "up_h : " << upidx_h << "  bs_h : " << bsidx_h
              << "  tail_h : " << tailidx_h
              << "  up_d : " << upidx_d << "  bs_d : " << bsidx_d
              << "  tail_d : " << tailidx_d
              << "  wdeg : " << wrkdeg << endl;

      if((mode == 1) || (mode == 2)) if(tailidx_h >= deg) break;
      if((mode == 0) || (mode == 2)) if(tailidx_d >= deg) break;

      wrkdeg = wrkdeg + 1 + wrkdeg/2;
      if(wrkdeg > deg) wrkdeg = deg;
   }
   if(vrblvl < 2)
   {
      double errsum = 0.0;

      cout << scientific << setprecision(16); // just in case vrblvl == 0
      cout << "The solution series : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << "sol[" << i << "][" << j << "] : " << sol[i][j] << endl;
            if((mode == 0) || (mode == 2))
            {
               cout << "x_d[" << i << "][" << j << "] : "
                              << input_d[i][j] << endl;
               errsum += abs(sol[i][j] - input_d[i][j]);
            }
            if((mode == 1) || (mode == 2))
            {
               cout << "x_h[" << i << "][" << j << "] : "
                              << input_h[i][j] << endl;
               errsum += abs(sol[i][j] - input_h[i][j]);
            }
         }
      }
      cout << "error : " << errsum << endl;
   }
   return 0;
}
