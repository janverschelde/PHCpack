// The file cmplx_newton_method.cpp defines the functions with prototypes in
// the file cmplx_newton_method.h.

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
#include "unimodular_matrices.h"
#include "dbl_monomial_systems.h"
#include "dbl_systems_host.h"
#include "dbl_systems_kernels.h"
#include "dbl_bals_host.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_kernels.h"
#include "dbl_newton_testers.h"
#include "write_newton_times.h"

using namespace std;

int cmplx_errors_funjacrhs
 ( int dim, int deg,
   double **funvalre_h, double **funvalim_h,
   double **funvalre_d, double **funvalim_d,
   double ***jacvalre_h, double ***jacvalim_h,
   double ***jacvalre_d, double ***jacvalim_d,
   double **rhsre_h, double **rhsim_h,
   double **rhsre_d, double **rhsim_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-8;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = cmplx_error2sum
              (dim,degp1,funvalre_h,funvalim_h,
                         funvalre_d,funvalim_d,"funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = cmplx_error3sum
              (degp1,dim,dim,jacvalre_h,jacvalim_h,
                             jacvalre_d,jacvalim_d,"jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   errsum = cmplx_error2sum
               (degp1,dim,rhsre_h,rhsim_h,
                          rhsre_d,rhsim_d,"rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int cmplx_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d,
   double **Qre_h, double **Qim_h, double **Qre_d, double **Qim_d,
   double **Rre_h, double **Rim_h, double **Rre_d, double **Rim_d,
   double **urhsre_h, double **urhsim_h,
   double **urhsre_d, double **urhsim_d,
   double **solre_h, double **solim_h, double **solre_d, double **solim_d,
   int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-8;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   errsum = cmplx_error2sum(dim,dim,Qre_h,Qim_h,Qre_d,Qim_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = cmplx_error2sum(dim,dim,Rre_h,Rim_h,Rre_d,Rim_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = cmplx_error2sum
               (degp1,dim,urhsre_h,urhsim_h,
                          urhsre_d,urhsim_d,"urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = cmplx_error2sum
               (degp1,dim,solre_h,solim_h,
                          solre_d,solim_d,"sol",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = cmplx_error2sum
               (dim,degp1,inputre_h,inputim_h,
                          inputre_d,inputim_d,"input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int cmplx_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d,
   double **funvalre_h, double **funvalim_h,
   double **funvalre_d, double **funvalim_d,
   double ***jacvalre_h, double ***jacvalim_h,
   double ***jacvalre_d, double ***jacvalim_d,
   double **rhsre_h, double **rhsim_h, double **rhsre_d, double **rhsim_d,
   double **urhsre_h, double **urhsim_h, double **urhsre_d, double **urhsim_d,
   double **solre_h, double **solim_h, double **solre_d, double **solim_d,
   double **Qre_h, double **Qim_h, double **Qre_d, double **Qim_d, 
   double **Rre_h, double **Rim_h, double **Rre_d, double **Rim_d,
   double *workvecre, double *workvecim,
   double **resvecre, double **resvecim, double *resmax,
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
         for(int j=0; j<dim; j++)
         {
            urhsre_h[i][j] = rhsre_h[i][j]; urhsim_h[i][j] = rhsim_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solre_h[i][j] = 0.0; solim_h[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling CPU_cmplx_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_cmplx_qrbs_solve
         (dim,degp1,oldtail,jacvalre_h,jacvalim_h,urhsre_h,urhsim_h,
          solre_h,solim_h,Qre_h,Qim_h,Rre_h,Rim_h,
          workvecre,workvecim,zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_cmplx_linear_residue ..." << endl;

         CPU_cmplx_linear_residue
            (dim,degp1,*tailidx_h-1,jacvalre_h,jacvalim_h,rhsre_h,rhsim_h,
             solre_h,solim_h,resvecre,resvecim,resmax,vrblvl);
         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmax << endl;
      }
      cmplx_update_series
         (dim,degp1,*tailidx_h-1,inputre_h,inputim_h,solre_h,solim_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhsre_d[i][j] = rhsre_d[i][j]; urhsim_d[i][j] = rhsim_d[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solre_d[i][j] = 0.0; solim_d[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling GPU_cmplx_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_cmplx_bals_solve
         (dim,degp1,szt,nbt,oldtail,jacvalre_d,jacvalim_d,
          Qre_d,Qim_d,Rre_d,Rim_d,urhsre_d,urhsim_d,solre_d,solim_d,
          zeroQ_d,noqr_d,upidx_d,bsidx_d,&newtail,
          totqrlapsedms,totqtblapsedms,totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_cmplx_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_cmplx_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,jacvalre_d,jacvalim_d,
             rhsre_d,rhsim_d,solre_d,solim_d,resvecre,resvecim,resmax,
             &elapsedms,&addcnt,&mulcnt,vrblvl);
         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmax;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;
         *totreslapsedms += elapsedms;
      }
      cmplx_update_series
         (dim,degp1,*tailidx_d-1,inputre_d,inputim_d,solre_d,solim_d,vrblvl);
   }
   return 0;
}

int cmplx_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbre, double **mbim, double dpr,
   double ***cffre, double ***cffim, double **accre, double **accim,
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
   double *workvecre, double *workvecim,
   double **resvecre, double **resvecim, double *resmax,
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
         cout << "calling CPU_cmplx_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx_unit_series_vector(dim,deg,cffre[0],cffim[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.
         CPU_cmplx_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffre[0],cffim[0],accre[0],accim[0],inputre_h,inputim_h,
             outputre_h,outputim_h,vrblvl);
      }
      else
         CPU_cmplx_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,cffre,cffim,accre,accim,
             inputre_h,inputim_h,funvalre_h,funvalim_h,
             jacvalre_h,jacvalim_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_cmplx_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx_unit_series_vector(dim,deg,cffre[0],cffim[0]);
         // reset coefficients
/*
         GPU_cmplx_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffre[0],cffim[0],accre[0],accim[0],inputre_d,inputim_d,
             outputre_d,outputim_d,totcnvlapsedms,vrblvl);
 */
         GPU_cmplxvectorized_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffre[0],cffim[0],accre[0],accim[0],inputre_d,inputim_d,
             outputre_d,outputim_d,totcnvlapsedms,vrblvl);
      }
      else
         GPU_cmplx_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cffre,cffim,
             inputre_d,inputim_d,outputre_d,outputim_d,
             funvalre_d,funvalim_d,jacvalre_d,jacvalim_d,
             totcnvlapsedms,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;

      double errsum = cmplx_error3sum
                         (dim,dim+1,degp1,
                          outputre_h,outputim_h,
                          outputre_d,outputim_d,"output",vrblvl);
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
         cmplx_define_rhs
            (dim,degp1,mbre,mbim,funvalre_h,funvalim_h,rhsre_h,rhsim_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalre_h[i][j][k] = 0.0; jacvalim_h[i][j][k] = 0.0;
               }

         cmplx_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbre,mbim,dpr,
             outputre_h,outputim_h,funvalre_h,funvalim_h,
             rhsre_h,rhsim_h,jacvalre_h,jacvalim_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         cmplx_define_rhs
            (dim,degp1,mbre,mbim,funvalre_d,funvalim_d,rhsre_d,rhsim_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalre_d[i][j][k] = 0.0; jacvalim_d[i][j][k] = 0.0;
               }
         cmplx_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbre,mbim,dpr,
             outputre_d,outputim_d,funvalre_d,funvalim_d,
             rhsre_d,rhsim_d,jacvalre_d,jacvalim_d,vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = cmplx_errors_funjacrhs
         (dim,deg,funvalre_h,funvalim_h,funvalre_d,funvalim_d,
          jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
          rhsre_h,rhsim_h,rhsre_d,rhsim_d,vrblvl);
   }
   cmplx_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputre_h,inputim_h,inputre_d,inputim_d,
       funvalre_h,funvalim_h,funvalre_d,funvalim_d,
       jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
       rhsre_h,rhsim_h,rhsre_d,rhsim_d,urhsre_h,
       urhsim_h,urhsre_d,urhsim_d,solre_h,solim_h,solre_d,solim_d,
       Qre_h, Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
       workvecre,workvecim,resvecre,resvecim,resmax,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return cmplx_errors_inurhsQRsol(dim,deg,
         inputre_h,inputim_h,inputre_d,inputim_d,Qre_h,Qim_h,Qre_d,Qim_d,
         Rre_h,Rim_h,Rre_d,Rim_d,urhsre_h,urhsim_h,urhsre_d,urhsim_d,
         solre_h,solim_h,solre_d,solim_d,vrblvl);
   }
   else
      return 0;
}

int cmplx_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cstre, double **cstim,
   double ***cffre, double ***cffim, double dpr,
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
   double *workvecre, double *workvecim,
   double **resvecre, double **resvecim, double *resmax,
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
         cout << "calling CPU_cmplx_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_cmplx_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
            cffre[i],cffim[i],inputre_h,inputim_h,
            outputre_h[i],outputim_h[i],&lapsed,0); // no output

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
         cout << "calling GPU_cmplx_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ComplexConvolutionJobs cnvjobs(dim);
         ComplexIncrementJobs incjobs(cnvjobs,false);
         ComplexAdditionJobs addjobs(dim,nbr[i]);

         make_all_complex_jobs
            (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,false);

         GPU_cmplxvectorized_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
             cffre[i],cffim[i],inputre_d,inputim_d,
             outputre_d[i],outputim_d[i],cnvjobs,incjobs,addjobs,
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

      errsum = cmplx_error3sum
                  (dim,dim+1,degp1,outputre_h,outputim_h,
                                   outputre_d,outputim_d,"output",vrblvl);
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
      cmplx_map_evaldiff_output
         (dim,deg,outputre_h,outputim_h,funvalre_h,funvalim_h,
          jacvalre_h,jacvalim_h,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsre_h[i][j] = -funvalre_h[j][i];
            rhsim_h[i][j] = -funvalim_h[j][i];
         }
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx_map_evaldiff_output
         (dim,deg,outputre_d,outputim_d,funvalre_d,funvalim_d,
          jacvalre_d,jacvalim_d,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsre_d[i][j] = -funvalre_d[j][i];
            rhsim_d[i][j] = -funvalim_d[j][i];
         }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = cmplx_errors_funjacrhs
         (dim,deg,funvalre_h,funvalim_h,funvalre_d,funvalim_d,
                  jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
          rhsre_h,rhsim_h,rhsre_d,rhsim_d,vrblvl);
   }
   cmplx_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputre_h,inputim_h,inputre_d,inputim_d,
       funvalre_h,funvalim_h,funvalre_d,funvalim_d,
       jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
       rhsre_h,rhsim_h,rhsre_d,rhsim_d,
       urhsre_h,urhsim_h,urhsre_d,urhsim_d,
       solre_h,solim_h,solre_d,solim_d,
       Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
       workvecre,workvecim,resvecre,resvecim,resmax,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
      return cmplx_errors_inurhsQRsol(dim,deg,
                inputre_h,inputim_h,inputre_d,inputim_d,
                Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
                urhsre_h,urhsim_h,urhsre_d,urhsim_d,
                solre_h,solim_h,solre_d,solim_d,vrblvl);
   else
      return 0;
}

int cmplx_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d,
   double ***outputre_h, double ***outputim_h,
   double ***outputre_d, double ***outputim_d,
   double **funvalre_h, double **funvalim_h,
   double **funvalre_d, double **funvalim_d,
   double ***jacvalre_h, double ***jacvalim_h,
   double ***jacvalre_d, double ***jacvalim_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputre_h[i] = new double[degp1];
         inputim_h[i] = new double[degp1];

         outputre_h[i] = new double*[dim+1];
         outputim_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputre_h[i][j] = new double[degp1];
            outputim_h[i][j] = new double[degp1];
         }
         funvalre_h[i] = new double[degp1];
         funvalim_h[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalre_h[i] = new double*[dim];
         jacvalim_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalre_h[i][j] = new double[dim];
            jacvalim_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputre_d[i] = new double[degp1];
         inputim_d[i] = new double[degp1];

         outputre_d[i] = new double*[dim+1];
         outputim_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputre_d[i][j] = new double[degp1];
            outputim_d[i][j] = new double[degp1];
         }
         funvalre_d[i] = new double[degp1];
         funvalim_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalre_d[i] = new double*[dim];
         jacvalim_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalre_d[i][j] = new double[dim];
            jacvalim_d[i][j] = new double[dim];
         }
      }
   }
   return 0;
}

int cmplx_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhsre_h, double **rhsim_h, double **rhsre_d, double **rhsim_d,
   double **urhsre_h, double **urhsim_h, double **urhsre_d, double **urhsim_d,
   double **Qre_h, double **Qim_h, double **Qre_d, double **Qim_d,
   double **Rre_h, double **Rim_h, double **Rre_d, double **Rim_d,
   double **solre_h, double **solim_h, double **solre_d, double **solim_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) 
      {
         rhsre_h[i] = new double[dim];
         rhsim_h[i] = new double[dim];
         urhsre_h[i] = new double[dim];
         urhsim_h[i] = new double[dim];
         solre_h[i] = new double[dim];
         solim_h[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qre_h[i] = new double[dim];
         Qim_h[i] = new double[dim];
         Rre_h[i] = new double[dim];
         Rim_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) 
      {
         rhsre_d[i] = new double[dim];
         rhsim_d[i] = new double[dim];
         urhsre_d[i] = new double[dim];
         urhsim_d[i] = new double[dim];
         solre_d[i] = new double[dim];
         solim_d[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qre_d[i] = new double[dim];
         Qim_d[i] = new double[dim];
         Rre_d[i] = new double[dim];
         Rim_d[i] = new double[dim];
      }
   }
   return 0;
}

void cmplx_start_setup
 ( int dim, int deg, double **testsolre, double **testsolim,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d, int mode, int vrblvl )
{
   double *start0re = new double[dim];
   double *start0im = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0re[i] = testsolre[i][0];
      start0im[i] = testsolim[i][0]; 
   }
   if((mode == 1) || (mode == 2))
      cmplx_start_series_vector(dim,deg,start0re,start0im,inputre_h,inputim_h);
   else
      cmplx_start_series_vector(dim,deg,start0re,start0im,inputre_d,inputim_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
         {
            inputre_d[i][j] = inputre_h[i][j];
            inputim_d[i][j] = inputim_h[i][j];
         }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;

      if((mode == 1) || (mode == 2))
      {
         for(int i=0; i<dim; i++)
            cout << i << " : "
                 << inputre_h[i][0] << "  " << inputim_h[i][0] << endl;
      }
      else
      {
         for(int i=0; i<dim; i++)
            cout << i << " : "
                 << inputre_d[i][0] << "  " << inputim_d[i][0] << endl;
      }
   }
}

void cmplx_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA, double ***cffre, double ***cffim,
   double **testsolre, double **testsolim,
   double **mbrhsre, double **mbrhsim,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   double *angles = new double[dim];

   for(int i=0; i<dim; i++)
   {
      testsolre[i] = new double[degp1];
      testsolim[i] = new double[degp1];
   }
   make_complex_exponentials(dim,deg,angles,testsolre,testsolim);

   // compute the right hand sides via evaluation

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
   if(nbrcol == 1)
      evaluate_complex_monomials
         (dim,deg,rowsA,testsolre,testsolim,mbrhsre,mbrhsim);
   else
      evaluate_complex_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffre,cffim,
          testsolre,testsolim,mbrhsre,mbrhsim,vrblvl);
 
   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsre[i][j] << "  " << mbrhsim[i][j] << endl;
   }
   cmplx_start_setup
      (dim,deg,testsolre,testsolim,inputre_h,inputim_h,inputre_d,inputim_d,
       mode,vrblvl);
}

void cmplx_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstre, double **cstim, double ***cffre, double ***cffim,
   double **testsolre, double **testsolim,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d,
   double ***outputre_h, double ***outputim_h,
   double ***outputre_d, double ***outputim_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   double *angles = new double[dim];

   for(int i=0; i<dim; i++)
   {
      testsolre[i] = new double[degp1];
      testsolim[i] = new double[degp1];
   }
   make_complex_exponentials(dim,deg,angles,testsolre,testsolim);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_cmplx_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
            cffre[i],cffim[i],testsolre,testsolim,
            outputre_h[i],outputim_h[i],&lapsed,0); // no output

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
         {
            cstre[i][j] = cstre[i][j] - outputre_h[i][dim][j];
            cstim[i][j] = cstim[i][j] - outputim_h[i][dim][j];
         }
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_cmplx_poly_evaldiff
              (dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
               cffre[i],cffim[i],testsolre,testsolim,
               outputre_h[i],outputim_h[i],&lapsed,0); // no output

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
               errsum = errsum + outputre_h[i][dim][j]
                               + outputim_h[i][dim][j]; 

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

         ComplexConvolutionJobs cnvjobs(dim);
         ComplexIncrementJobs incjobs(cnvjobs,false);
         ComplexAdditionJobs addjobs(dim,nbr[i]);

         make_all_complex_jobs
            (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrb);

         GPU_cmplxvectorized_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
             cffre[i],cffim[i],testsolre,testsolim,
             outputre_d[i],outputim_d[i],cnvjobs,incjobs,addjobs,
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
         {
            cstre[i][j] = cstre[i][j] - outputre_d[i][dim][j];
            cstim[i][j] = cstim[i][j] - outputim_d[i][dim][j];
         }
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double cnvlapms,addlapms,timelapms_d,walltimes_d;
   
            ComplexConvolutionJobs cnvjobs(dim);
            ComplexIncrementJobs incjobs(cnvjobs,false);
            ComplexAdditionJobs addjobs(dim,nbr[i]);

            make_all_complex_jobs
               (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrb);

            GPU_cmplxvectorized_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],cstre[i],cstim[i],
                cffre[i],cffim[i],testsolre,testsolim,
                outputre_d[i],outputim_d[i],cnvjobs,incjobs,addjobs,
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
               errsum = errsum + outputre_d[i][dim][j]
                               + outputim_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
      }
   }
   cmplx_start_setup
      (dim,deg,testsolre,testsolim,inputre_h,inputim_h,inputre_d,inputim_d,
       mode,vrblvl);
}

int cmplx_error_testsol
 ( int dim, int deg, int mode,
   double **testsolre, double **testsolim,
   double **inputre_h, double **inputim_h,
   double **inputre_d, double **inputim_d )
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
                        << testsolre[i][j] << "  "
                        << testsolim[i][j] << endl;
         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputre_d[i][j] << "  "
                           << inputim_d[i][j] << endl;
            errsum += abs(testsolre[i][j] - inputre_d[i][j])
                    + abs(testsolim[i][j] - inputim_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputre_h[i][j] << "  "
                           << inputim_h[i][j] << endl;
            errsum += abs(testsolre[i][j] - inputre_h[i][j])
                    + abs(testsolim[i][j] - inputim_h[i][j]);
         }
      }
   }
   cout << scientific << setprecision(2)
        << "error : " << errsum << endl;

   return (errsum > 1.0e-8);
}

int test_cmplx_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **accre = new double*[dim+1]; // accumulate series in one column
   double **accim = new double*[dim+1];
   for(int i=0; i<=dim; i++)
   {
      accre[i] = new double[degp1];
      accim[i] = new double[degp1];
   }
   double ***cffre = new double**[nbrcol]; // the coefficients of monomials
   double ***cffim = new double**[nbrcol]; 
   for(int i=0; i<nbrcol; i++)
   {
      cffre[i] = new double*[dim];
      cffim[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffre[i][j] = new double[degp1];
         cffim[i][j] = new double[degp1];
      }
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_complex_coefficients(nbrcol,dim,cffre,cffim);

   double **inputre_h;
   double **inputim_h;
   double **inputre_d;
   double **inputim_d;
   double ***outputre_h;
   double ***outputim_h;
   double ***outputre_d;
   double ***outputim_d;
   double **funvalre_h;
   double **funvalim_h;
   double **funvalre_d;
   double **funvalim_d;
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalre_h;
   double ***jacvalim_h;
   double ***jacvalre_d;
   double ***jacvalim_d;

   if((mode == 1) || (mode == 2))
   {
      inputre_h = new double*[dim];
      inputim_h = new double*[dim];
      outputre_h = new double**[dim];
      outputim_h = new double**[dim];
      funvalre_h = new double*[dim];
      funvalim_h = new double*[dim];
      jacvalre_h = new double**[degp1];
      jacvalim_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputre_d = new double*[dim];
      inputim_d = new double*[dim];
      outputre_d = new double**[dim];
      outputim_d = new double**[dim];
      funvalre_d = new double*[dim];
      funvalim_d = new double*[dim];
      jacvalre_d = new double**[degp1];
      jacvalim_d = new double**[degp1];
   }
   cmplx_allocate_inoutfunjac
      (dim,deg,mode,inputre_h,inputim_h,inputre_d,inputim_d,
       outputre_h,outputim_h,outputre_d,outputim_d,
       funvalre_h,funvalim_h,funvalre_d,funvalim_d,
       jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d);

/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhsre_h;
   double **rhsim_h;
   double **rhsre_d;
   double **rhsim_d;
   double **urhsre_h;
   double **urhsim_h;
   double **urhsre_d;
   double **urhsim_d;
   double **solre_h;
   double **solim_h;
   double **solre_d;
   double **solim_d;
   double **Qre_h;
   double **Qim_h;
   double **Qre_d;
   double **Qim_d;
   double **Rre_h;
   double **Rim_h;
   double **Rre_d;
   double **Rim_d;

   if((mode == 1) || (mode == 2))
   {
      rhsre_h = new double*[degp1];
      rhsim_h = new double*[degp1];
      urhsre_h = new double*[degp1];
      urhsim_h = new double*[degp1];
      Qre_h = new double*[dim];
      Qim_h = new double*[dim];
      Rre_h = new double*[dim];
      Rim_h = new double*[dim];
      solre_h = new double*[degp1];
      solim_h = new double*[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      rhsre_d = new double*[degp1];
      rhsim_d = new double*[degp1];
      urhsre_d = new double*[degp1];
      urhsim_d = new double*[degp1];
      Qre_d = new double*[dim];
      Qim_d = new double*[dim];
      Rre_d = new double*[dim];
      Rim_d = new double*[dim];
      solre_d = new double*[degp1];
      solim_d = new double*[degp1];
   }
   cmplx_allocate_rhsqrsol
      (dim,deg,mode,rhsre_h,rhsim_h,rhsre_d,rhsim_d,
       urhsre_h,urhsim_h,urhsre_d,urhsim_d,
       Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
       solre_h,solim_h,solre_d,solim_d);
  
   double *workvecre = new double[dim];
   double *workvecim = new double[dim];
   double **resvecre = new double*[degp1];
   double **resvecim = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecre[i] = new double[dim];
      resvecim[i] = new double[dim];
   }
   double resmax;
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolre = new double*[dim];
   double **testsolim = new double*[dim];
   double **mbrhsre = new double*[dim];
   double **mbrhsim = new double*[dim];

   cmplx_column_setup
      (dim,deg,nbrcol,nvr,idx,rowsA,cffre,cffim,testsolre,testsolim,
       mbrhsre,mbrhsim,inputre_h,inputim_h,inputre_d,inputim_d,mode,vrblvl);

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

      cmplx_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,
          mbrhsre,mbrhsim,dpr,cffre,cffim,accre,accim,
          inputre_h,inputim_h,inputre_d,inputim_d,
          outputre_h,outputim_h,outputre_d,outputim_d,
          funvalre_h,funvalim_h,funvalre_d,funvalim_d,
          jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
          rhsre_h,rhsim_h,rhsre_d,rhsim_d,
          urhsre_h,urhsim_h,urhsre_d,urhsim_d,
          solre_h,solim_h,solre_d,solim_d,
          Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
          workvecre,workvecim,resvecre,resvecim,&resmax,
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

   cmplx_error_testsol
      (dim,deg,mode,testsolre,testsolim,
       inputre_h,inputim_h,inputre_d,inputim_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_cmplx_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;
   double **cstre = new double*[dim];
   double **cstim = new double*[dim];
   double ***cffre = new double**[dim];
   double ***cffim = new double**[dim];

   cmplx_make_coefficients
      (dim,deg,nbr,nvr,idx,cstre,cstim,cffre,cffim,vrblvl);
   // no randomization of the leading coefficients needed,
   // because no complex exponenials used for complex data

   double **inputre_h;
   double **inputim_h;
   double **inputre_d;
   double **inputim_d;
   double ***outputre_h;
   double ***outputim_h;
   double ***outputre_d;
   double ***outputim_d;
   double **funvalre_h;  // function values on host
   double **funvalim_h; 
   double **funvalre_d;  // function values on device
   double **funvalim_d;
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalre_h;
   double ***jacvalim_h;
   double ***jacvalre_d;
   double ***jacvalim_d;

   if((mode == 1) || (mode == 2))
   {
      inputre_h = new double*[dim];
      inputim_h = new double*[dim];
      outputre_h = new double**[dim];
      outputim_h = new double**[dim];
      funvalre_h = new double*[dim];
      funvalim_h = new double*[dim];
      jacvalre_h = new double**[degp1];
      jacvalim_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputre_d = new double*[dim];
      inputim_d = new double*[dim];
      outputre_d = new double**[dim];
      outputim_d = new double**[dim];
      funvalre_d = new double*[dim];
      funvalim_d = new double*[dim];
      jacvalre_d = new double**[degp1];
      jacvalim_d = new double**[degp1];
   }
   cmplx_allocate_inoutfunjac
      (dim,deg,mode,inputre_h,inputim_h,inputre_d,inputim_d,
       outputre_h,outputim_h,outputre_d,outputim_d,
       funvalre_h,funvalim_h,funvalre_d,funvalim_d,
       jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d);

   double **rhsre_h;
   double **rhsim_h;
   double **rhsre_d;
   double **rhsim_d;
   double **urhsre_h;
   double **urhsim_h;
   double **urhsre_d;
   double **urhsim_d;
   double **Qre_h;
   double **Qim_h;
   double **Qre_d;
   double **Qim_d;
   double **Rre_h;
   double **Rim_h;
   double **Rre_d;
   double **Rim_d;
   double **solre_h;
   double **solim_h;
   double **solre_d;
   double **solim_d;

   if((mode == 1) || (mode == 2))
   {
      rhsre_h = new double*[degp1];
      rhsim_h = new double*[degp1];
      urhsre_h = new double*[degp1];
      urhsim_h = new double*[degp1];
      Qre_h = new double*[dim];
      Qim_h = new double*[dim];
      Rre_h = new double*[dim];
      Rim_h = new double*[dim];
      solre_h = new double*[degp1];
      solim_h = new double*[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      rhsre_d = new double*[degp1];
      rhsim_d = new double*[degp1];
      urhsre_d = new double*[degp1];
      urhsim_d = new double*[degp1];
      Qre_d = new double*[dim];
      Qim_d = new double*[dim];
      Rre_d = new double*[dim];
      Rim_d = new double*[dim];
      solre_d = new double*[degp1];
      solim_d = new double*[degp1];
   }
   cmplx_allocate_rhsqrsol
      (dim,deg,mode,rhsre_h,rhsim_h,rhsre_d,rhsim_d,
       urhsre_h,urhsim_h,urhsre_d,urhsim_d,
       Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
       solre_h,solim_h,solre_d,solim_d);

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolre = new double*[dim];
   double **testsolim = new double*[dim];

   cmplx_row_setup
      (dim,deg,nbr,nvr,idx,cstre,cstim,cffre,cffim,testsolre,testsolim,
       inputre_h,inputim_h,inputre_d,inputim_d,
       outputre_h,outputim_h,outputre_d,outputim_d,mode,vrblvl);

   // alocating some extra work space

   double *workvecre = new double[dim];
   double *workvecim = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **resvecre = new double*[degp1];
   double **resvecim = new double*[degp1];
   for(int i=0; i<degp1; i++)
   {
      resvecre[i] = new double[dim];
      resvecim[i] = new double[dim];
   }
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

      cmplx_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,
          nbr,nvr,idx,cstre,cstim,cffre,cffim,dpr,
          inputre_h,inputim_h,inputre_d,inputim_d,
          outputre_h,outputim_h,outputre_d,outputim_d,
          funvalre_h,funvalim_h,funvalre_d,funvalim_d,
          jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
          rhsre_h,rhsim_h,rhsre_d,rhsim_d,urhsre_h,
          urhsim_h,urhsre_d,urhsim_d,solre_h,solim_h,solre_d,solim_d,
          Qre_h, Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
          workvecre,workvecim,resvecre,resvecim,&resmax,
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

   cmplx_error_testsol
      (dim,deg,mode,testsolre,testsolim,
       inputre_h,inputim_h,inputre_d,inputim_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
