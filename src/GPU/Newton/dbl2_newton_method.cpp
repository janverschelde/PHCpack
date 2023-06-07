// The file dbl2_newton_method.cpp defines the functions with prototypes in
// the file dbl2_newton_method.h.

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
#include "unimodular_matrices.h"
#include "random_numbers.h"
#include "random2_vectors.h"
#include "dbl2_indexed_coefficients.h"
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_polynomials_host.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"
#include "job_makers.h"
#include "dbl2_polynomials_kernels.h"
#include "dbl2_factorizations.h"
#include "dbl2_monomial_systems.h"
#include "dbl2_bals_host.h"
#include "dbl2_bals_kernels.h"
#include "dbl2_tail_kernels.h"
#include "dbl2_systems_host.h"
#include "dbl2_systems_kernels.h"
#include "dbl_bals_flopcounts.h"
#include "dbl2_newton_testers.h"
#include "write_newton_times.h"

using namespace std;

int dbl2_errors_funjacrhs
 ( int dim, int deg,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-20;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = dbl2_error2sum
               (dim,degp1,funvalhi_h,funvallo_h,
                          funvalhi_d,funvallo_d,"funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = dbl2_error3sum
               (degp1,dim,dim,jacvalhi_h,jacvallo_h,
                              jacvalhi_d,jacvallo_d,"jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   errsum = dbl2_error2sum
               (degp1,dim,rhshi_h,rhslo_h,rhshi_d,rhslo_d,"rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl2_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d,
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d,
   int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-20;
   double errsum = 0.0;
   int fail = 0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   errsum = dbl2_error2sum(dim,dim,Qhi_h,Qlo_h,Qhi_d,Qlo_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = dbl2_error2sum(dim,dim,Rhi_h,Rlo_h,Rhi_d,Rlo_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = dbl2_error2sum(degp1,dim,urhshi_h,urhslo_h,
                                     urhshi_d,urhslo_d,"urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = dbl2_error2sum(degp1,dim,solhi_h,sollo_h,
                                     solhi_d,sollo_d,"sol",vrblvl);
   cout << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = dbl2_error2sum(dim,degp1,inputhi_h,inputlo_h,
                                     inputhi_d,inputlo_d,"input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl2_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d,
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
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
            urhshi_h[i][j] = rhshi_h[i][j];
            urhslo_h[i][j] = rhslo_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhi_h[i][j] = 0.0;
            sollo_h[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling CPU_dbl2_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_dbl2_qrbs_solve
         (dim,degp1,oldtail,jacvalhi_h,jacvallo_h,
          urhshi_h,urhslo_h,solhi_h,sollo_h,
          Qhi_h,Qlo_h,Rhi_h,Rlo_h,workvechi,workveclo,
          zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl2_linear_residue ..." << endl;

         CPU_dbl2_linear_residue
            (dim,degp1,*tailidx_h-1,jacvalhi_h,jacvallo_h,rhshi_h,rhslo_h,
             solhi_h,sollo_h,resvechi,resveclo,resmaxhi,resmaxlo,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi << endl;
      }
      dbl2_update_series
         (dim,degp1,*tailidx_h-1,inputhi_h,inputlo_h,
          solhi_h,sollo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhshi_d[i][j] = rhshi_d[i][j];
            urhslo_d[i][j] = rhslo_d[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhi_d[i][j] = 0.0;
            sollo_d[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling GPU_dbl2_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_dbl2_bals_solve
         (dim,degp1,szt,nbt,oldtail,jacvalhi_d,jacvallo_d,
          Qhi_d,Qlo_d,Rhi_d,Rlo_d,urhshi_d,urhslo_d,solhi_d,sollo_d,
          zeroQ_d,noqr_d,upidx_d,bsidx_d,&newtail,
          totqrlapsedms,totqtblapsedms,totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl2_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_dbl2_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,jacvalhi_d,jacvallo_d,
             rhshi_d,rhslo_d,solhi_d,sollo_d,resvechi,resveclo,
             resmaxhi,resmaxlo,&elapsedms,&addcnt,&mulcnt,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;
         *totreslapsedms += elapsedms;
      }
      dbl2_update_series
         (dim,degp1,*tailidx_d-1,inputhi_d,inputlo_d,solhi_d,sollo_d,vrblvl);
   }
   return 0;
}

int dbl2_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbhi, double **mblo, double dpr,
   double ***cffhi, double ***cfflo, double **acchi, double **acclo,
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
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   bool* zeroQ_h, bool *noqr_h, bool* zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_dbl2_evaluate_monomials ..." << endl;
 
      if(nbrcol == 1)
      {
         dbl2_unit_series_vector(dim,deg,cffhi[0],cfflo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_dbl2_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhi[0],cfflo[0],acchi[0],acclo[0],
             inputhi_h,inputlo_h,outputhi_h,outputlo_h,vrblvl);
      }
      else
         CPU_dbl2_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,cffhi,cfflo,acchi,acclo,
             inputhi_h,inputlo_h,funvalhi_h,funvallo_h,
             jacvalhi_h,jacvallo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl2_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl2_unit_series_vector(dim,deg,cffhi[0],cfflo[0]);
         // reset coefficients

         GPU_dbl2_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhi[0],cfflo[0],acchi[0],acclo[0],
             inputhi_d,inputlo_d,outputhi_d,outputlo_d,totcnvlapsedms,vrblvl);
      }
      else
         GPU_dbl2_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cffhi,cfflo,
             inputhi_d,inputlo_d,outputhi_d,outputlo_d,
             funvalhi_d,funvallo_d,jacvalhi_d,jacvallo_d,
             totcnvlapsedms,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;

      errsum = dbl2_error3sum
                  (dim,dim+1,degp1,outputhi_h,outputlo_h,
                                   outputhi_d,outputlo_d,"output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;
   cout << flush;

   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl2_define_rhs
            (dim,degp1,mbhi,mblo,funvalhi_h,funvallo_h,
             rhshi_h,rhslo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  cout << "  i = " << i
                       << "  j = " << j
                       << "  k = " << k << endl << flush;

                  jacvalhi_h[i][j][k] = 0.0; jacvallo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbhi,mblo,dpr,outputhi_h,outputlo_h,
             funvalhi_h,funvallo_h,rhshi_h,rhslo_h,
             jacvalhi_h,jacvallo_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl2_define_rhs
            (dim,degp1,mbhi,mblo,funvalhi_d,funvallo_d,
             rhshi_d,rhslo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalhi_d[i][j][k] = 0.0; jacvallo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbhi,mblo,dpr,outputhi_d,outputlo_d,
             funvalhi_d,funvallo_d,rhshi_d,rhslo_d,jacvalhi_d,jacvallo_d,
             vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl2_errors_funjacrhs
         (dim,deg,funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,vrblvl);
   }
   dbl2_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
       jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
       rhshi_h,rhslo_h,rhshi_d,rhslo_d,
       urhshi_h,urhslo_h,urhshi_d,urhslo_d,
       solhi_h,sollo_h,solhi_d,sollo_d,
       Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
       workvechi,workveclo,resvechi,resveclo,resmaxhi,resmaxlo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
      return dbl2_errors_inurhsQRsol(dim,deg,
                inputhi_h,inputlo_h,inputhi_d,inputlo_d,
                Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
                urhshi_h,urhslo_h,urhshi_d,urhslo_d,
                solhi_h,sollo_h,solhi_d,sollo_d,vrblvl);
   else
      return 0;
}

int dbl2_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **csthi, double **cstlo,
   double ***cffhi, double ***cfflo, double dpr,
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
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
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
         cout << "calling CPU_dbl2_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl2_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],csthi[i],cstlo[i],cffhi[i],cfflo[i],
            inputhi_h,inputlo_h,outputhi_h[i],outputlo_h[i],&lapsed,0);

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
         cout << "calling GPU_dbl2_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ConvolutionJobs cnvjobs(dim);
         AdditionJobs addjobs(dim,nbr[i]);

         make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,false);

         GPU_dbl2_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],csthi[i],cstlo[i],
             cffhi[i],cfflo[i],inputhi_d,inputlo_d,
             outputhi_d[i],outputlo_d[i],cnvjobs,addjobs,
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

      errsum = dbl2_error3sum
                  (dim,dim+1,degp1,outputhi_h,outputlo_h,
                                   outputhi_d,outputlo_d,"output",vrblvl);
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
      dbl2_map_evaldiff_output
         (dim,deg,outputhi_h,outputlo_h,funvalhi_h,funvallo_h,
          jacvalhi_h,jacvallo_h,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshi_h[i][j] = -funvalhi_h[j][i];
            rhslo_h[i][j] = -funvallo_h[j][i];
         }
   }
   if((mode == 0) || (mode == 2))
   {
      dbl2_map_evaldiff_output
         (dim,deg,outputhi_d,outputlo_d,funvalhi_d,funvallo_d,
          jacvalhi_d,jacvallo_d,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshi_d[i][j] = -funvalhi_d[j][i];
            rhslo_d[i][j] = -funvallo_d[j][i];
         }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl2_errors_funjacrhs(dim,deg,
                     funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
                     jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
                     rhshi_h,rhslo_h,rhshi_d,rhslo_d,vrblvl);
   }
   dbl2_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
       jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
       rhshi_h,rhslo_h,rhshi_d,rhslo_d,
       urhshi_h,urhslo_h,urhshi_d,urhslo_d,
       solhi_h,sollo_h,solhi_d,sollo_d,
       Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
       workvechi,workveclo,resvechi,resveclo,resmaxhi,resmaxlo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
      return dbl2_errors_inurhsQRsol(dim,deg,
                 inputhi_h,inputlo_h,inputhi_d,inputlo_d,
                 Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
                 urhshi_h,urhslo_h,urhshi_d,urhslo_d,
                 solhi_h,sollo_h,solhi_d,sollo_d,vrblvl);
   else
      return 0;
}

int dbl2_allocate_inoutfunjac
 ( int dim, int deg, int mode,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h, 
   double ***outputhi_d, double ***outputlo_d,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputhi_h[i] = new double[degp1];
         inputlo_h[i] = new double[degp1];
         outputhi_h[i] = new double*[dim+1];
         outputlo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhi_h[i][j] = new double[degp1];
            outputlo_h[i][j] = new double[degp1];
         }
         funvalhi_h[i] = new double[degp1];
         funvallo_h[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhi_h[i] = new double*[dim];
         jacvallo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhi_h[i][j] = new double[dim];
            jacvallo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputhi_d[i] = new double[degp1];
         inputlo_d[i] = new double[degp1];
         outputhi_d[i] = new double*[dim+1];
         outputlo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhi_d[i][j] = new double[degp1];
            outputlo_d[i][j] = new double[degp1];
         }
         funvalhi_d[i] = new double[degp1];
         funvallo_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhi_d[i] = new double*[dim];
         jacvallo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhi_d[i][j] = new double[dim];
            jacvallo_d[i][j] = new double[dim];
         }
      }
   }
   return 0;
}

int dbl2_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d, 
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) 
      {
         rhshi_h[i] = new double[dim];
         rhslo_h[i] = new double[dim];
         urhshi_h[i] = new double[dim];
         urhslo_h[i] = new double[dim];
         solhi_h[i] = new double[dim];
         sollo_h[i] = new double[dim];
      }
      for(int i=0; i<dim; i++) 
      {
         Qhi_h[i] = new double[dim];
         Qlo_h[i] = new double[dim];
         Rhi_h[i] = new double[dim];
         Rlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++)
      {
         rhshi_d[i] = new double[dim];
         rhslo_d[i] = new double[dim];
         urhshi_d[i] = new double[dim];
         urhslo_d[i] = new double[dim];
         solhi_d[i] = new double[dim];
         sollo_d[i] = new double[dim];
      }
      for(int i=0; i<dim; i++) 
      {
         Qhi_d[i] = new double[dim];
         Qlo_d[i] = new double[dim];
         Rhi_d[i] = new double[dim];
         Rlo_d[i] = new double[dim];
      }
   }
   return 0;
}

void dbl2_start_setup
 ( int dim, int deg, double **testsolhi, double **testsollo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d, int mode, int vrblvl )
{
   double *start0hi = new double[dim];
   double *start0lo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0hi[i] = testsolhi[i][0];
      start0lo[i] = testsollo[i][0];
   }
   if((mode == 1) || (mode == 2))
      real2_start_series_vector
         (dim,deg,start0hi,start0lo,inputhi_h,inputlo_h);
   else
      real2_start_series_vector
         (dim,deg,start0hi,start0lo,inputhi_d,inputlo_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
         {
            inputhi_d[i][j] = inputhi_h[i][j];
            inputlo_d[i][j] = inputlo_h[i][j];
         }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;

      if((mode == 1) || (mode == 2))
         for(int i=0; i<dim; i++)
            cout << i << " : " << inputhi_h[i][0] << "  "
                               << inputlo_h[i][0] << endl;
      else
         for(int i=0; i<dim; i++)
            cout << i << " : " << inputhi_d[i][0] << "  "
                               << inputlo_d[i][0] << endl;
   }
}

void dbl2_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA, double ***cffhi, double ***cfflo,
   double **testsolhi, double **testsollo,
   double **mbrhshi, double **mbrhslo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolhi[i] = new double[degp1];
      testsollo[i] = new double[degp1];
   }
   make_real2_exponentials(dim,deg,testsolhi,testsollo);

   // compute the right hand sides via evaluation

   for(int i=0; i<dim; i++)
   {
      mbrhshi[i] = new double[degp1];
      mbrhslo[i] = new double[degp1];

      mbrhshi[i][0] = 1.0;     // initialize product to one
      mbrhslo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhshi[i][k] = 0.0;
         mbrhslo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_real2_monomials
         (dim,deg,rowsA,testsolhi,testsollo,mbrhshi,mbrhslo);
   else
      evaluate_real2_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffhi,cfflo,
          testsolhi,testsollo,mbrhshi,mbrhslo,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshi[i][j] << "  "
                 << mbrhslo[i][j] << endl;
   }
   dbl2_start_setup
      (dim,deg,testsolhi,testsollo,inputhi_h,inputlo_h,
       inputhi_d,inputlo_d,mode,vrblvl);
}

void dbl2_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthi, double **cstlo, double ***cffhi, double ***cfflo, 
   double **testsolhi, double **testsollo, 
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h,
   double ***outputhi_d, double ***outputlo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++) 
   {
      testsolhi[i] = new double[degp1];
      testsollo[i] = new double[degp1];
   }
   make_real2_exponentials(dim,deg,testsolhi,testsollo);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl2_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
            csthi[i],cstlo[i],cffhi[i],cfflo[i],
            testsolhi,testsollo,outputhi_h[i],outputlo_h[i],&lapsed,0);

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
         for(int j=0; j<=deg; j++) // cst[i][j] -= output_h[i][dim][j];
            ddf_dec(&csthi[i][j],&cstlo[i][j],
                    outputhi_h[i][dim][j],outputlo_h[i][dim][j]);
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_dbl2_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
               csthi[i],cstlo[i],cffhi[i],cfflo[i],
               testsolhi,testsollo,outputhi_h[i],outputlo_h[i],&lapsed,0);

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
               errsum = errsum + outputhi_h[i][dim][j]
                               + outputlo_h[i][dim][j];

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

         GPU_dbl2_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],csthi[i],cstlo[i],
             cffhi[i],cfflo[i],testsolhi,testsollo,
             outputhi_d[i],outputlo_d[i],cnvjobs,addjobs,
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
         for(int j=0; j<=deg; j++) // cst[i][j] -= output_d[i][dim][j];
            ddf_dec(&csthi[i][j],&cstlo[i][j],
                    outputhi_d[i][dim][j],outputlo_d[i][dim][j]);
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

            GPU_dbl2_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],csthi[i],cstlo[i],
                cffhi[i],cfflo[i],testsolhi,testsollo,
                outputhi_d[i],outputlo_d[i],cnvjobs,addjobs,
                &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0);

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
               errsum = errsum + outputhi_d[i][dim][j]
                               + outputlo_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "residual of test solution : " << errsum << endl;
      }
   }
   dbl2_start_setup(dim,deg,
       testsolhi,testsollo,inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       mode,vrblvl);
}

int dbl2_error_testsol
 ( int dim, int deg, int mode, double **testsolhi, double **testsollo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d )
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
              << testsolhi[i][j] << "  "
              << testsollo[i][j] << endl;

         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhi_d[i][j] << "  "
                           << inputlo_d[i][j] << endl;
            errsum += abs(testsolhi[i][j] - inputhi_d[i][j])
                    + abs(testsollo[i][j] - inputlo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhi_h[i][j] << "  "
                           << inputlo_h[i][j] << endl;
            errsum += abs(testsolhi[i][j] - inputhi_h[i][j])
                    + abs(testsollo[i][j] - inputlo_h[i][j]);
         }
      }
   }
   cout << "error : " << errsum << endl;

   return (errsum > 1.0e-20);
}

int test_dbl2_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **acchi = new double*[dim+1]; // accumulated power series
   double **acclo = new double*[dim+1]; // in one column

   for(int i=0; i<=dim; i++)
   {
      acchi[i] = new double[degp1];
      acclo[i] = new double[degp1];
   }
   double ***cffhi = new double**[nbrcol]; // coefficients of monomials
   double ***cfflo = new double**[nbrcol];
   for(int i=0; i<nbrcol; i++)
   {
      cffhi[i] = new double*[dim];
      cfflo[i] = new double*[dim];
      for(int j=0; j<dim; j++)
      {
         cffhi[i][j] = new double[degp1];
         cfflo[i][j] = new double[degp1];
      }
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_real2_coefficients(nbrcol,dim,cffhi,cfflo);

   double **inputhi_h;
   double **inputlo_h;
   double **inputhi_d;
   double **inputlo_d;
   double ***outputhi_h;
   double ***outputlo_h;
   double ***outputhi_d;
   double ***outputlo_d;
   double **funvalhi_h;
   double **funvallo_h;
   double **funvalhi_d;
   double **funvallo_d;
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h;
   double ***jacvallo_h;
   double ***jacvalhi_d;
   double ***jacvallo_d;

   if((mode == 1) || (mode == 2))
   {
      inputhi_h = new double*[dim];
      inputlo_h = new double*[dim];
      outputhi_h = new double**[dim];
      outputlo_h = new double**[dim];
      funvalhi_h = new double*[dim];
      funvallo_h = new double*[dim];
      jacvalhi_h = new double**[degp1];
      jacvallo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputhi_d = new double*[dim];
      inputlo_d = new double*[dim];
      outputhi_d = new double**[dim];
      outputlo_d = new double**[dim];
      funvalhi_d = new double*[dim];
      funvallo_d = new double*[dim];
      jacvalhi_d = new double**[degp1];
      jacvallo_d = new double**[degp1];
   }
   dbl2_allocate_inoutfunjac(dim,deg,mode,
       inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       outputhi_h,outputlo_h,outputhi_d,outputlo_d,
       funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
       jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhshi_h;
   double **rhslo_h;
   double **rhshi_d;
   double **rhslo_d;
   double **urhshi_h;
   double **urhslo_h;
   double **urhshi_d;
   double **urhslo_d;
   double **solhi_h;
   double **sollo_h;
   double **solhi_d;
   double **sollo_d;
   double **Qhi_h;
   double **Qlo_h;
   double **Qhi_d;
   double **Qlo_d;
   double **Rhi_h;
   double **Rlo_h;
   double **Rhi_d;
   double **Rlo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshi_h = new double*[degp1];
      rhslo_h = new double*[degp1];
      urhshi_h = new double*[degp1];
      urhslo_h = new double*[degp1];
      solhi_h = new double*[degp1];
      sollo_h = new double*[degp1];
      Qhi_h = new double*[dim];
      Qlo_h = new double*[dim];
      Rhi_h = new double*[dim];
      Rlo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhshi_d = new double*[degp1];
      rhslo_d = new double*[degp1];
      urhshi_d = new double*[degp1];
      urhslo_d = new double*[degp1];
      solhi_d = new double*[degp1];
      sollo_d = new double*[degp1];
      Qhi_d = new double*[dim];
      Qlo_d = new double*[dim];
      Rhi_d = new double*[dim];
      Rlo_d = new double*[dim];
   }
   dbl2_allocate_rhsqrsol(dim,deg,mode,
       rhshi_h,rhslo_h,rhshi_d,rhslo_d,
       urhshi_h,urhslo_h,urhshi_d,urhslo_d,
       Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
       solhi_h,sollo_h,solhi_d,sollo_d);

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
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolhi = new double*[dim];
   double **testsollo = new double*[dim];
   double **mbrhshi = new double*[dim];
   double **mbrhslo = new double*[dim];

   dbl2_column_setup
      (dim,deg,nbrcol,nvr,idx,rowsA,cffhi,cfflo,testsolhi,testsollo,
       mbrhshi,mbrhslo,inputhi_h,inputlo_h,inputhi_d,inputlo_d,mode,vrblvl);

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

      dbl2_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,
          mbrhshi,mbrhslo,dpr,cffhi,cfflo,acchi,acclo,
          inputhi_h,inputlo_h,inputhi_d,inputlo_d,
          outputhi_h,outputlo_h,outputhi_d,outputlo_d,
          funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,urhshi_h,urhslo_h,urhshi_d,urhslo_d,
          solhi_h,sollo_h,solhi_d,sollo_d,
          Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
          workvechi,workveclo,resvechi,resveclo,&resmaxhi,&resmaxlo,
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

   dbl2_error_testsol
      (dim,deg,mode,testsolhi,testsollo,
       inputhi_h,inputlo_h,inputhi_d,inputlo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_dbl2_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;
   double **csthi = new double*[dim];
   double **cstlo = new double*[dim];
   double ***cffhi = new double**[dim];
   double ***cfflo = new double**[dim];

   dbl2_make_coefficients(dim,deg,nbr,nvr,idx,csthi,cstlo,cffhi,cfflo,vrblvl);
   // must randomize the leading coefficients,
   // because for exponenials all leading coefficients are one
   for(int i=0; i<dim; i++)
      for(int j=0; j<nbr[i]; j++)
         random_double_double(&cffhi[i][j][0],&cfflo[i][j][0]);

   double **inputhi_h;
   double **inputlo_h;
   double **inputhi_d;
   double **inputlo_d;
   double ***outputhi_h;
   double ***outputlo_h;
   double ***outputhi_d;
   double ***outputlo_d;
   double **funvalhi_h;  // function values on host
   double **funvallo_h; 
   double **funvalhi_d;  // function values on device
   double **funvallo_d;
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h;
   double ***jacvallo_h;
   double ***jacvalhi_d;
   double ***jacvallo_d;

   if((mode == 1) || (mode == 2))
   {
      inputhi_h = new double*[dim];
      inputlo_h = new double*[dim];
      outputhi_h = new double**[dim];
      outputlo_h = new double**[dim];
      funvalhi_h = new double*[dim];
      funvallo_h = new double*[dim];
      jacvalhi_h = new double**[degp1];
      jacvallo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputhi_d = new double*[dim];
      inputlo_d = new double*[dim];
      outputhi_d = new double**[dim];
      outputlo_d = new double**[dim];
      funvalhi_d = new double*[dim];
      funvallo_d = new double*[dim];
      jacvalhi_d = new double**[degp1];
      jacvallo_d = new double**[degp1];
   }
   dbl2_allocate_inoutfunjac(dim,deg,mode,
       inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       outputhi_h,outputlo_h,outputhi_d,outputlo_d,
       funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
       jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d);

   double **rhshi_h;
   double **rhslo_h;
   double **rhshi_d;
   double **rhslo_d;
   double **urhshi_h;
   double **urhslo_h;
   double **urhshi_d;
   double **urhslo_d;
   double **Qhi_h;
   double **Qlo_h;
   double **Qhi_d;
   double **Qlo_d;
   double **Rhi_h;
   double **Rlo_h;
   double **Rhi_d;
   double **Rlo_d;
   double **solhi_h;
   double **sollo_h;
   double **solhi_d;
   double **sollo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshi_h = new double*[degp1];
      rhslo_h = new double*[degp1];
      urhshi_h = new double*[degp1];
      urhslo_h = new double*[degp1];
      solhi_h = new double*[degp1];
      sollo_h = new double*[degp1];
      Qhi_h = new double*[dim];
      Qlo_h = new double*[dim];
      Rhi_h = new double*[dim];
      Rlo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhshi_d = new double*[degp1];
      rhslo_d = new double*[degp1];
      urhshi_d = new double*[degp1];
      urhslo_d = new double*[degp1];
      solhi_d = new double*[degp1];
      sollo_d = new double*[degp1];
      Qhi_d = new double*[dim];
      Qlo_d = new double*[dim];
      Rhi_d = new double*[dim];
      Rlo_d = new double*[dim];
   }
   dbl2_allocate_rhsqrsol(dim,deg,mode,
       rhshi_h,rhslo_h,rhshi_d,rhslo_d,
       urhshi_h,urhslo_h,urhshi_d,urhslo_d,
       Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
       solhi_h,sollo_h,solhi_d,sollo_d);

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolhi = new double*[dim];
   double **testsollo = new double*[dim];

   dbl2_row_setup
      (dim,deg,nbr,nvr,idx,csthi,cstlo,cffhi,cfflo,
       testsolhi,testsollo,inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       outputhi_h,outputlo_h,outputhi_d,outputlo_d,mode,vrblvl);

   // alocating some extra work space

   double *workvechi = new double[dim];
   double *workveclo = new double[dim];
   double **resvechi = new double*[degp1];
   double **resveclo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechi[i] = new double[dim];
      resveclo[i] = new double[dim];
   }
   double resmaxhi,resmaxlo;

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

      dbl2_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,
          nbr,nvr,idx,csthi,cstlo,cffhi,cfflo,dpr,
          inputhi_h,inputlo_h,inputhi_d,inputlo_d,
          outputhi_h,outputlo_h,outputhi_d,outputlo_d,
          funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,
          urhshi_h,urhslo_h,urhshi_d,urhslo_d,
          solhi_h,sollo_h,solhi_d,sollo_d,
          Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
          workvechi,workveclo,resvechi,resveclo,&resmaxhi,&resmaxlo,
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

   dbl2_error_testsol(dim,deg,mode,
       testsolhi,testsollo,inputhi_h,inputlo_h,inputhi_d,inputlo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
