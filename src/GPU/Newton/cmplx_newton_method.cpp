// The file cmplx_newton_method.cpp defines the functions with prototypes in
// the file cmplx_newton_method.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "random_numbers.h"
#include "random_monomials.h"
#include "dbl_factorizations.h"
#include "unimodular_matrices.h"
#include "dbl_monomial_systems.h"
#include "dbl_systems_host.h"
#include "dbl_systems_kernels.h"
#include "dbl_bals_host.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_kernels.h"
#include "dbl_newton_testers.h"

using namespace std;

void cmplx_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **mbre, double **mbim, double dpr,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d,
   double ***outputre_h, double ***outputim_h,
   double ***outputre_d, double ***outputim_d,
   double **funvalre_h, double **funvalim_h,
   double **funvalre_d, double **funvalim_d,
   double ***jacvalre_h, double ***jacvalim_h,
   double ***jacvalre_d, double ***jacvalim_d,
   double **rhsre_h, double **rhsim_h, double **rhsre_d, double **rhsim_d,
   double **urhsre_h, double **urhsim_h, double **urhsre_d, double **urhsim_d,
   double **solre_h, double **solim_h, double **solre_d, double **solim_d,
   double **Qre_h, double **Qim_h, double **Qre_d, double **Qim_d, 
   double **Rre_h, double **Rim_h, double **Rre_d, double **Rim_d,
   double **workmatre, double **workmatim,
   double *workvecre, double *workvecim,
   double **resvecre, double **resvecim, double *resmax,
   bool *noqr_h, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      // The series coefficients accumulate common factors,
      // initially the coefficients are set to one.
      cmplx_unit_series_vector(dim,deg,cffre,cffim);

      if(vrblvl > 0)
         cout << "calling CPU_cmplx_evaluate_monomials ..." << endl;

      CPU_cmplx_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffre,cffim,accre,accim,inputre_h,inputim_h,
          outputre_h,outputim_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx_unit_series_vector(dim,deg,cffre,cffim); // reset coefficients

      if(vrblvl > 0)
         cout << "calling GPU_cmplx_evaluate_monomials ..." << endl;

      GPU_cmplx_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffre,cffim,accre,accim,inputre_d,inputim_d,
          outputre_d,outputim_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;

      double errsum = cmplx_error3sum
                         (dim,dim+1,degp1,
                          outputre_h,outputim_h,
                          outputre_d,outputim_d,"output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(16);
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalre_h[i][j][k] = 0.0; jacvalim_h[i][j][k] = 0.0;
            jacvalre_d[i][j][k] = 0.0; jacvalim_d[i][j][k] = 0.0;
         }

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      cmplx_linearize_evaldiff_output
         (dim,degp1,nvr,idx,mbre,mbim,dpr,
          outputre_h,outputim_h,funvalre_h,funvalim_h,
          rhsre_h,rhsim_h,jacvalre_h,jacvalim_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx_linearize_evaldiff_output
         (dim,degp1,nvr,idx,mbre,mbim,dpr,
          outputre_d,outputim_d,funvalre_d,funvalim_d,
          rhsre_d,rhsim_d,jacvalre_d,jacvalim_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = cmplx_error2sum
                         (dim,degp1,funvalre_h,funvalim_h,
                                    funvalre_d,funvalim_d,"funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = cmplx_error3sum
                  (degp1,dim,dim,jacvalre_h,jacvalim_h,
                                 jacvalre_d,jacvalim_d,"jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = cmplx_error2sum(degp1,dim,rhsre_h,rhsim_h,
                                         rhsre_d,rhsim_d,"rhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhsre_h[i][j] = rhsre_h[i][j]; urhsim_h[i][j] = rhsim_h[i][j];
         urhsre_d[i][j] = rhsre_d[i][j]; urhsim_d[i][j] = rhsim_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solre_h[i][j] = 0.0; solim_h[i][j] = 0.0;
         solre_d[i][j] = 0.0; solim_d[i][j] = 0.0;
      }

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_cmplx_qrbs_solve ..." << endl;

      CPU_cmplx_qrbs_solve
         (dim,degp1,jacvalre_h,jacvalim_h,urhsre_h,urhsim_h,
          solre_h,solim_h,workmatre,workmatim,Qre_h,Qim_h,Rre_h,Rim_h,
          workvecre,workvecim,noqr_h,upidx_h,bsidx_h,vrblvl);
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_cmplx_linear_residue ..." << endl;

         CPU_cmplx_linear_residue
            (dim,degp1,jacvalre_h,jacvalim_h,rhsre_h,rhsim_h,
             solre_h,solim_h,resvecre,resvecim,resmax,vrblvl);
         cout << "maximum residual : " << *resmax << endl;
      }
      cmplx_update_series
         (dim,degp1,inputre_h,inputim_h,solre_h,solim_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_cmplx_bals_solve ..." << endl;

      GPU_cmplx_bals_solve
         (dim,degp1,szt,nbt,jacvalre_d,jacvalim_d,Qre_d,Qim_d,Rre_d,Rim_d,
          urhsre_d,urhsim_d,solre_d,solim_d,noqr_d,upidx_d,bsidx_d,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling GPU_cmplx_linear_residue ..." << endl;

         double elapsedms;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_cmplx_linear_residue
            (dim,degp1,szt,nbt,jacvalre_d,jacvalim_d,rhsre_d,rhsim_d,
             solre_d,solim_d,resvecre,resvecim,resmax,
             &elapsedms,&addcnt,&mulcnt,vrblvl);
         cout << "maximum residual : " << *resmax << endl;
      }
      cmplx_update_series
         (dim,degp1,inputre_d,inputim_d,solre_d,solim_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      errsum = cmplx_error2sum(dim,dim,Qre_h,Qim_h,Qre_d,Qim_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = cmplx_error2sum(dim,dim,Rre_h,Rim_h,Rre_d,Rim_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = cmplx_error2sum(degp1,dim,urhsre_h,urhsim_h,
                                         urhsre_d,urhsim_d,"urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = cmplx_error2sum(degp1,dim,solre_h,solim_h,
                                         solre_d,solim_d,"sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = cmplx_error2sum(dim,degp1,inputre_h,inputim_h,
                                         inputre_d,inputim_d,"input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputre_h = new double*[dim];
   double **inputim_h = new double*[dim];
   double **inputre_d = new double*[dim];
   double **inputim_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputre_h[i] = new double[degp1];
       inputim_h[i] = new double[degp1];
       inputre_d[i] = new double[degp1];
       inputim_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *accre = new double[degp1]; // accumulated power series
   double *accim = new double[degp1];
   double **cffre = new double*[dim]; // the coefficients of monomials
   double **cffim = new double*[dim]; 
   for(int i=0; i<dim; i++)
   {
      cffre[i] = new double[degp1];
      cffim[i] = new double[degp1];
   }
   double ***outputre_h = new double**[dim];
   double ***outputim_h = new double**[dim];
   double ***outputre_d = new double**[dim];
   double ***outputim_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputre_h[i] = new double*[dim+1];
      outputim_h[i] = new double*[dim+1];
      outputre_d[i] = new double*[dim+1];
      outputim_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputre_h[i][j] = new double[degp1];
         outputim_h[i][j] = new double[degp1];
         outputre_d[i][j] = new double[degp1];
         outputim_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalre_h = new double*[dim];
   double **funvalim_h = new double*[dim];
   double **funvalre_d = new double*[dim];
   double **funvalim_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalre_h[i] = new double[degp1];
      funvalim_h[i] = new double[degp1];
      funvalre_d[i] = new double[degp1];
      funvalim_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalre_h = new double**[degp1];
   double ***jacvalim_h = new double**[degp1];
   double ***jacvalre_d = new double**[degp1];
   double ***jacvalim_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalre_h[i] = new double*[dim];
      jacvalim_h[i] = new double*[dim];
      jacvalre_d[i] = new double*[dim];
      jacvalim_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalre_h[i][j] = new double[dim];
         jacvalim_h[i][j] = new double[dim];
         jacvalre_d[i][j] = new double[dim];
         jacvalim_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solre_h = new double*[degp1];
   double **solim_h = new double*[degp1];
   double **solre_d = new double*[degp1];
   double **solim_d = new double*[degp1];

   for(int i=0; i<degp1; i++) 
   {
      solre_h[i] = new double[dim];
      solim_h[i] = new double[dim];
      solre_d[i] = new double[dim];
      solim_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsre_h = new double*[degp1];
   double **rhsim_h = new double*[degp1];
   double **rhsre_d = new double*[degp1];
   double **rhsim_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsre_h[i] = new double[dim];
      rhsim_h[i] = new double[dim];
      rhsre_d[i] = new double[dim];
      rhsim_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmatre = new double*[dim];
   double **workmatim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmatre[i] = new double[dim];
      workmatim[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvecre = new double[dim];
   double *workvecim = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **urhsre_h = new double*[degp1];
   double **urhsim_h = new double*[degp1];
   double **urhsre_d = new double*[degp1];
   double **urhsim_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhsre_h[i] = new double[dim];
      urhsim_h[i] = new double[dim];
      urhsre_d[i] = new double[dim];
      urhsim_d[i] = new double[dim];
   }
   double **resvecre = new double*[degp1];
   double **resvecim = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecre[i] = new double[dim];
      resvecim[i] = new double[dim];
   }
   double resmax;
   double **Qre_h = new double*[dim];
   double **Qim_h = new double*[dim];
   double **Qre_d = new double*[dim];
   double **Qim_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qre_h[i] = new double[dim];
      Qim_h[i] = new double[dim];
      Qre_d[i] = new double[dim];
      Qim_d[i] = new double[dim];
   }
   double **Rre_h = new double*[dim];
   double **Rim_h = new double*[dim];
   double **Rre_d = new double*[dim];
   double **Rim_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rre_h[i] = new double[dim];
      Rim_h[i] = new double[dim];
      Rre_d[i] = new double[dim];
      Rim_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   double *angles = new double[dim];
   double **solre = new double*[dim];
   double **solim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solre[i] = new double[degp1];
      solim[i] = new double[degp1];
   }
   make_complex_exponentials(dim,deg,angles,solre,solim);

   // compute the right hand sides via evaluation

   double **mbrhsre = new double*[dim];
   double **mbrhsim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhsre[i] = new double[degp1];
      mbrhsim[i] = new double[degp1];

      mbrhsre[i][0] = 1.0;     // initialize product to one
      mbrhsim[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhsre[i][k] = 0.0;
         mbrhsim[i][k] = 0.0;
      }
   }
   evaluate_complex_monomials(dim,deg,rowsA,solre,solim,mbrhsre,mbrhsim);
   
   double *start0re = new double[dim];
   double *start0im = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0re[i] = solre[i][0];
      start0im[i] = solim[i][0]; 
   }
   cmplx_start_series_vector(dim,deg,start0re,start0im,inputre_h,inputim_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputre_d[i][j] = inputre_h[i][j];
         inputim_d[i][j] = inputim_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << inputre_h[i][0] << "  " << inputim_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool noqr_h = false;
   bool noqr_d = false;

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;

      cmplx_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          mbrhsre,mbrhsim,dpr,cffre,cffim,accre,accim,
          inputre_h,inputim_h,inputre_d,inputim_d,
          outputre_h,outputim_h,outputre_d,outputim_d,
          funvalre_h,funvalim_h,funvalre_d,funvalim_d,
          jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
          rhsre_h,rhsim_h,rhsre_d,rhsim_d,
          urhsre_h,urhsim_h,urhsre_d,urhsim_d,
          solre_h,solim_h,solre_d,solim_d,
          Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
          workmatre,workmatim,workvecre,workvecim,resvecre,resvecim,&resmax,
          &noqr_h,&noqr_d,&upidx_h,&bsidx_h,&upidx_d,&bsidx_d,vrblvl,mode);

      if(vrblvl > 0)
         cout << "upidx_h : " << upidx_h << "  bsidx_h : " << bsidx_h
              << "  upidx_d : " << upidx_d << "  bsidx_d : " << bsidx_d
              << "  deg : " << deg << endl;

      if((mode == 1) || (mode == 2)) if(bsidx_h >= deg) break;
      if((mode == 0) || (mode == 2)) if(bsidx_d >= deg) break;
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
            cout << "sol[" << i << "][" << j << "] : "
                           << solre[i][j] << "  "
                           << solim[i][j] << endl;
            if((mode == 0) || (mode == 2))
            {
               cout << "x_d[" << i << "][" << j << "] : "
                              << inputre_d[i][j] << "  "
                              << inputim_d[i][j] << endl;
               errsum += abs(solre[i][j] - inputre_d[i][j])
                       + abs(solim[i][j] - inputim_d[i][j]);
            }
            if((mode == 1) || (mode == 2))
            {
               cout << "x_h[" << i << "][" << j << "] : "
                              << inputre_h[i][j] << "  "
                              << inputim_h[i][j] << endl;
               errsum += abs(solre[i][j] - inputre_h[i][j])
                       + abs(solim[i][j] - inputim_h[i][j]);
            }
         }
      }
      cout << "error : " << errsum << endl;
   }
   return 0;
}
