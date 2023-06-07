// The file dbl_newton_method.cpp defines the functions with prototypes in
// the file dbl_newton_method.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "random_numbers.h"
#include "random_monomials.h"
#include "dbl_indexed_coefficients.h"
#include "dbl_polynomials_host.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"
#include "job_makers.h"
#include "dbl_polynomials_kernels.h"
#include "dbl_factorizations.h"
#include "dbl_monomial_systems.h"
#include "unimodular_matrices.h"
#include "dbl_systems_host.h"
#include "dbl_systems_kernels.h"
#include "dbl_bals_host.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_kernels.h"
#include "dbl_newton_testers.h"
#include "write_newton_times.h"

using namespace std;

int dbl_errors_funjacrhs
 ( int dim, int deg, double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d,
   double **rhs_h, double **rhs_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-8;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = dbl_error2sum(dim,degp1,funval_h,funval_d,"funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = dbl_error3sum
               (degp1,dim,dim,jacval_h,jacval_d,"jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   errsum = dbl_error2sum(degp1,dim,rhs_h,rhs_d,"rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl_errors_inurhsQRsol
 ( int dim, int deg, double **input_h, double **input_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-8;
   int fail = 0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   double errsum = dbl_error2sum(dim,dim,Q_h,Q_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = dbl_error2sum(dim,dim,R_h,R_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = dbl_error2sum(degp1,dim,urhs_h,urhs_d,"urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = dbl_error2sum(degp1,dim,sol_h,sol_d,"sol",vrblvl);
   cout << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = dbl_error2sum(dim,degp1,input_h,input_d,"input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d, double **input_h, double **input_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

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
          zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl_linear_residue ..." << endl;

         CPU_dbl_linear_residue
            (dim,degp1,*tailidx_h-1,jacval_h,rhs_h,sol_h,
             resvec,resmax,vrblvl);
         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmax << endl;
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
          zeroQ_d,noqr_d,upidx_d,bsidx_d,&newtail,
          totqrlapsedms,totqtblapsedms,totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_dbl_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,jacval_d,rhs_d,sol_d,
             resvec,resmax,&elapsedms,&addcnt,&mulcnt,vrblvl);
         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmax;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;
         *totreslapsedms += elapsedms;
      }
      dbl_update_series(dim,degp1,*tailidx_d-1,input_d,sol_d,vrblvl);
   }
   return 0;
}

int dbl_column_newton_qrstep
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
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
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
             input_d,output_d,totcnvlapsedms,vrblvl);
      }
      else
         GPU_dbl_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cff,input_d,output_d,
             funval_d,jacval_d,totcnvlapsedms,vrblvl);
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

      cout << scientific << setprecision(3);
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
      int fail = dbl_errors_funjacrhs
         (dim,deg,funval_h,funval_d,jacval_h,jacval_d,rhs_h,rhs_d,vrblvl);
   }
   dbl_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,input_h,input_d,
       funval_h,funval_d,jacval_h,jacval_d,
       rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
       Q_h,Q_d,R_h,R_d,workvec,resvec,resmax,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
      return dbl_errors_inurhsQRsol(dim,deg,
         input_h,input_d,Q_h,Q_d,R_h,R_d,urhs_h,urhs_d,sol_h,sol_d,vrblvl);
   else
      return 0;
}

int dbl_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cst, double ***cff, double dpr,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_dbl_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],input_h,
            output_h[i],&lapsed,0); // no output: vrblvl is zero

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << lapsed << " seconds." << endl;

         timelapsed_h += lapsed;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_h << " seconds." << endl;
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ConvolutionJobs cnvjobs(dim);
         AdditionJobs addjobs(dim,nbr[i]);

         make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,false);

         GPU_dbl_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],input_d,
             output_d[i],cnvjobs,addjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0); // vrblvl);

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << walltimes_d << " seconds." << endl;

         timelapsed_d += walltimes_d;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_d << " seconds." << endl;
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << scientific << setprecision(16)
           << "comparing CPU with GPU evaluations ... " << endl;
   
      double errsum = 0.0;

      errsum = dbl_error3sum
                  (dim,dim+1,degp1,output_h,output_d,"output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(2);
      cout << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0)
      cout << "mapping the output to values of function and matrix series ..."
           << endl;

   if((mode == 1) || (mode == 2))
   {
      dbl_map_evaldiff_output(dim,deg,output_h,funval_h,jacval_h,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++) rhs_h[i][j] = -funval_h[j][i];
   }
   if((mode == 0) || (mode == 2))
   {
      dbl_map_evaldiff_output(dim,deg,output_d,funval_d,jacval_d,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++) rhs_d[i][j] = -funval_d[j][i];
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl_errors_funjacrhs
         (dim,deg,funval_h,funval_d,jacval_h,jacval_d,rhs_h,rhs_d,vrblvl);
   }
   dbl_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,input_h,input_d,
       funval_h,funval_d,jacval_h,jacval_d,
       rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
       Q_h,Q_d,R_h,R_d,workvec,resvec,resmax,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
      return dbl_errors_inurhsQRsol(dim,deg,
         input_h,input_d,Q_h,Q_d,R_h,R_d,urhs_h,urhs_d,sol_h,sol_d,vrblvl);
   else
      return 0;
}

int dbl_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **input_h, double **input_d,
   double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         input_h[i] = new double[degp1];
         output_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
            output_h[i][j] = new double[degp1];

         funval_h[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacval_h[i] = new double*[dim];
         for(int j=0; j<dim; j++)
            jacval_h[i][j] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         input_d[i] = new double[degp1];
         output_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
            output_d[i][j] = new double[degp1];

         funval_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacval_d[i] = new double*[dim];
         for(int j=0; j<dim; j++)
            jacval_d[i][j] = new double[dim];
      }
   }
   return 0;
}

int dbl_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhs_h, double **rhs_d, double **urhs_h, double **urhs_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double **sol_h, double **sol_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) rhs_h[i] = new double[dim];
      for(int i=0; i<degp1; i++) urhs_h[i] = new double[dim];
      for(int i=0; i<dim; i++) Q_h[i] = new double[dim];
      for(int i=0; i<dim; i++) R_h[i] = new double[dim];
      for(int i=0; i<degp1; i++) sol_h[i] = new double[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) rhs_d[i] = new double[dim];
      for(int i=0; i<degp1; i++) urhs_d[i] = new double[dim];
      for(int i=0; i<dim; i++) Q_d[i] = new double[dim];
      for(int i=0; i<dim; i++) R_d[i] = new double[dim];
      for(int i=0; i<degp1; i++) sol_d[i] = new double[dim];
   }
   return 0;
}

void dbl_start_setup
 ( int dim, int deg, double **testsol,
   double **input_h, double **input_d, int mode, int vrblvl )
{
   double *start0 = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0[i] = testsol[i][0];
   }
   if((mode == 1) || (mode == 2))
      real_start_series_vector(dim,deg,start0,input_h);
   else
      real_start_series_vector(dim,deg,start0,input_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
            input_d[i][j] = input_h[i][j];
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;

      if((mode == 1) || (mode == 2))
         for(int i=0; i<dim; i++)
            cout << i << " : " << input_h[i][0] << endl;
      else
         for(int i=0; i<dim; i++)
            cout << i << " : " << input_d[i][0] << endl;
   }
}

void dbl_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA, double ***cff,
   double **testsol, double **mbrhs,
   double **input_h, double **input_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++) testsol[i] = new double[degp1];

   make_real_exponentials(dim,deg,testsol);

   // compute the right hand sides via evaluation

   for(int i=0; i<dim; i++)
   {
      mbrhs[i] = new double[degp1];

      mbrhs[i][0] = 1.0;     // initialize product to one

      for(int k=1; k<degp1; k++) mbrhs[i][k] = 0.0;
   }
   if(nbrcol == 1)
      evaluate_real_monomials(dim,deg,rowsA,testsol,mbrhs);
   else
      evaluate_real_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cff,testsol,mbrhs,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhs[i][j] << endl;
   }
   dbl_start_setup(dim,deg,testsol,input_h,input_d,mode,vrblvl);
}

void dbl_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cst, double ***cff, double **testsol, 
   double **input_h, double **input_d,
   double ***output_h, double ***output_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++) testsol[i] = new double[degp1];

   make_real_exponentials(dim,deg,testsol);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],testsol,
            output_h[i],&lapsed,0); // no output: vrblvl is zero

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << lapsed << " seconds." << endl;

         timelapsed_h += lapsed;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_h << " seconds." << endl;

      for(int i=0; i<dim; i++) // adjust constant coefficients
      {
         for(int j=0; j<=deg; j++)
            cst[i][j] = cst[i][j] - output_h[i][dim][j];
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_dbl_poly_evaldiff
              (dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],testsol,
               output_h[i],&lapsed,0); // no output: vrblvl is zero

            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << lapsed << " seconds." << endl;

            timelapsed_h += lapsed;
         }
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_h << " seconds." << endl;

         double errsum = 0.0;

         for(int i=0; i<dim; i++)
            for(int j=0; j<=deg; j++)
               errsum = errsum + output_h[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
      }
   }
   else // GPU is faster
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the device ..." << endl;

      const bool vrb = false; // no output (vrblvl > 1);
      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ConvolutionJobs cnvjobs(dim);
         AdditionJobs addjobs(dim,nbr[i]);

         make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrb);

         GPU_dbl_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],testsol,
             output_d[i],cnvjobs,addjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0); // vrblvl);

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << walltimes_d << " seconds." << endl;

         timelapsed_d += walltimes_d;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_d << " seconds." << endl;

      for(int i=0; i<dim; i++) // adjust constant coefficients
      {
         for(int j=0; j<=deg; j++)
            cst[i][j] = cst[i][j] - output_d[i][dim][j];
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double cnvlapms,addlapms,timelapms_d,walltimes_d;
   
            ConvolutionJobs cnvjobs(dim);
            AdditionJobs addjobs(dim,nbr[i]);

            make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrb);

            GPU_dbl_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],cst[i],cff[i],testsol,
                output_d[i],cnvjobs,addjobs,
                &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0); // vrblvl);

            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << walltimes_d << " seconds." << endl;
   
            timelapsed_d += walltimes_d;
         }
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_d << " seconds." << endl;

         double errsum = 0.0;

         for(int i=0; i<dim; i++)
            for(int j=0; j<=deg; j++)
               errsum = errsum + output_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
      }
   }
   dbl_start_setup(dim,deg,testsol,input_h,input_d,mode,vrblvl);
}

int dbl_error_testsol
 ( int dim, int deg, int mode,
   double **testsol, double **input_h, double **input_d )
{
   double errsum = 0.0;

   cout << scientific << setprecision(16); // just in case vrblvl == 0
   cout << "The solution series : " << endl;

   for(int j=0; j<=deg; j++)
   {
      cout << "coefficient of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "sol[" << i << "][" << j << "] : "
              << testsol[i][j] << endl;

         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << input_d[i][j] << endl;
            errsum += abs(testsol[i][j] - input_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << input_h[i][j] << endl;
            errsum += abs(testsol[i][j] - input_h[i][j]);
         }
      }
   }
   cout << scientific << setprecision(2)
        << "error : " << errsum << endl;

   return (errsum > 1.0e-8);
}

int test_dbl_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **acc = new double*[dim+1]; // accumulate series in one column
   for(int i=0; i<=dim; i++) acc[i] = new double[degp1];
   double ***cff = new double**[nbrcol];
   for(int i=0; i<nbrcol; i++)
   {
      cff[i] = new double*[dim]; // the coefficients of monomials
      for(int j=0; j<dim; j++) cff[i][j] = new double[degp1];
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_real_coefficients(nbrcol,dim,cff);

   double **input_h;
   double **input_d;
   double ***output_h;
   double ***output_d;
   double **funval_h;  // function values on host
   double **funval_d;  // function values on device
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacval_h;
   double ***jacval_d;

   if((mode == 1) || (mode == 2))
   {
      input_h = new double*[dim];
      output_h = new double**[dim];
      funval_h = new double*[dim];
      jacval_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      input_d = new double*[dim];
      output_d = new double**[dim];
      funval_d = new double*[dim];
      jacval_d = new double**[degp1];
   }
   dbl_allocate_inoutfunjac
      (dim,deg,mode,input_h,input_d,output_h,output_d,
       funval_h,funval_d,jacval_h,jacval_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhs_h;
   double **rhs_d;
   double **urhs_h;
   double **urhs_d;
   double **Q_h;
   double **Q_d;
   double **R_h;
   double **R_d;
   double **sol_h;
   double **sol_d;

   if((mode == 1) || (mode == 2))
   {
      rhs_h = new double*[degp1];
      urhs_h = new double*[degp1];
      Q_h = new double*[dim];
      R_h = new double*[dim];
      sol_h = new double*[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      rhs_d = new double*[degp1];
      urhs_d = new double*[degp1];
      Q_d = new double*[dim];
      R_d = new double*[dim];
      sol_d = new double*[degp1];
   }
   dbl_allocate_rhsqrsol
      (dim,deg,mode,rhs_h,rhs_d,urhs_h,urhs_d,Q_h,Q_d,R_h,R_d,sol_h,sol_d);

   double *workvec = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **resvec = new double*[degp1];
   for(int i=0; i<degp1; i++) resvec[i] = new double[dim];
   double resmax;
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsol = new double*[dim];
   double **mbrhs = new double*[dim];

   dbl_column_setup(dim,deg,nbrcol,nvr,idx,rowsA,cff,testsol,mbrhs,
                    input_h,input_d,mode,vrblvl);

   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = true;
   bool zeroQ_d = true;
   bool noqr_h = false;
   bool noqr_d = false;
   int tailidx_h = 1;
   int tailidx_d = 1;
   int wrkdeg = 0; // working degree of precision
   int stepcnt = 0;

   double totcnvlapsedms = 0.0;
   double totqrlapsedms = 0.0;
   double totqtblapsedms = 0.0;
   double totbslapsedms = 0.0;
   double totupdlapsedms = 0.0;
   double totreslapsedms = 0.0;

   struct timeval begintime,endtime; // wall clock time of computations
   gettimeofday(&begintime,0);

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      dbl_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,mbrhs,dpr,cff,acc,
          input_h,input_d,output_h,output_d,funval_h,funval_d,
          jacval_h,jacval_d,rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
          Q_h,Q_d,R_h,R_d,workvec,resvec,&resmax,
          &zeroQ_h,&noqr_h,&zeroQ_d,&noqr_d,
          &upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
          &totcnvlapsedms,&totqrlapsedms,&totqtblapsedms,&totbslapsedms,
          &totupdlapsedms,&totreslapsedms,vrblvl,mode);

      stepcnt = stepcnt + 1;

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
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   double walltimesec = seconds + microseconds*1.0e-6;

   dbl_error_testsol(dim,deg,mode,testsol,input_h,input_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_dbl_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;
   double **cst = new double*[dim];
   double ***cff = new double**[dim];

   dbl_make_coefficients(dim,deg,nbr,nvr,idx,cst,cff,vrblvl);
   // must randomize the leading coefficients,
   // because for exponenials all leading coefficients are one
   for(int i=0; i<dim; i++)
      for(int j=0; j<nbr[i]; j++) cff[i][j][0] = random_double();

   double **input_h;
   double **input_d;
   double ***output_h;
   double ***output_d;
   double **funval_h;  // function values on host
   double **funval_d;  // function values on device
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacval_h;
   double ***jacval_d;

   if((mode == 1) || (mode == 2))
   {
      input_h = new double*[dim];
      output_h = new double**[dim];
      funval_h = new double*[dim];
      jacval_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      input_d = new double*[dim];
      output_d = new double**[dim];
      funval_d = new double*[dim];
      jacval_d = new double**[degp1];
   }
   dbl_allocate_inoutfunjac
      (dim,deg,mode,input_h,input_d,output_h,output_d,
       funval_h,funval_d,jacval_h,jacval_d);

   double **rhs_h;
   double **rhs_d;
   double **urhs_h;
   double **urhs_d;
   double **Q_h;
   double **Q_d;
   double **R_h;
   double **R_d;
   double **sol_h;
   double **sol_d;

   if((mode == 1) || (mode == 2))
   {
      rhs_h = new double*[degp1];
      urhs_h = new double*[degp1];
      Q_h = new double*[dim];
      R_h = new double*[dim];
      sol_h = new double*[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      rhs_d = new double*[degp1];
      urhs_d = new double*[degp1];
      Q_d = new double*[dim];
      R_d = new double*[dim];
      sol_d = new double*[degp1];
   }
   dbl_allocate_rhsqrsol
      (dim,deg,mode,rhs_h,rhs_d,urhs_h,urhs_d,Q_h,Q_d,R_h,R_d,sol_h,sol_d);

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsol = new double*[dim];

   dbl_row_setup(dim,deg,nbr,nvr,idx,cst,cff,testsol,
                 input_h,input_d,output_h,output_d,mode,vrblvl);

// alocating some extra work space

   double *workvec = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **resvec = new double*[degp1];
   for(int i=0; i<degp1; i++) resvec[i] = new double[dim];
   double resmax;

   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = true;
   bool zeroQ_d = true;
   bool noqr_h = false;
   bool noqr_d = false;
   int tailidx_h = 1;
   int tailidx_d = 1;
   int wrkdeg = 0; // working degree of precision
   int stepcnt = 0;

   double totcnvlapsedms = 0.0;
   double totqrlapsedms = 0.0;
   double totqtblapsedms = 0.0;
   double totbslapsedms = 0.0;
   double totupdlapsedms = 0.0;
   double totreslapsedms = 0.0;

   struct timeval begintime,endtime; // wall clock time of computations
   gettimeofday(&begintime,0);

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      dbl_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,
          nbr,nvr,idx,cst,cff,dpr,
          input_h,input_d,output_h,output_d,funval_h,funval_d,
          jacval_h,jacval_d,rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
          Q_h,Q_d,R_h,R_d,workvec,resvec,&resmax,
          &zeroQ_h,&noqr_h,&zeroQ_d,&noqr_d,
          &upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
          &totcnvlapsedms,&totqrlapsedms,&totqtblapsedms,&totbslapsedms,
          &totupdlapsedms,&totreslapsedms,vrblvl,mode);

      stepcnt = stepcnt + 1;

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
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   double walltimesec = seconds + microseconds*1.0e-6;

   dbl_error_testsol(dim,deg,mode,testsol,input_h,input_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
