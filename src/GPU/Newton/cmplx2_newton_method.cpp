// The file cmplx2_newton_method.cpp defines the functions with prototypes in
// the file cmplx2_newton_method.h.

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
#include "random_monomials.h"
#include "double_double_functions.h"
#include "dbl2_indexed_coefficients.h"
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

int cmplx2_errors_funjacrhs
 ( int dim, int deg,
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
   double **rhsimhi_d, double **rhsimlo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-20;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = cmplx2_error2sum(dim,degp1,
               funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
               funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
               "funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = cmplx2_error3sum(degp1,dim,dim,
               jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
               jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
               "jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   errsum = cmplx2_error2sum(degp1,dim,
               rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
               rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
               "rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int cmplx2_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h, 
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d,
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double **urhsrehi_h, double **urhsrelo_h, 
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d, 
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h,
   double **solrehi_d, double **solrelo_d,
   double **solimhi_d, double **solimlo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-20;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   errsum = cmplx2_error2sum(dim,dim,
               Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
               Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = cmplx2_error2sum(dim,dim,
               Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
               Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = cmplx2_error2sum(degp1,dim,
               urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
               urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,"urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = cmplx2_error2sum(degp1,dim,
               solrehi_h,solrelo_h,solimhi_h,solimlo_h,
               solrehi_d,solrelo_d,solimhi_d,solimlo_d,"sol",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = cmplx2_error2sum(dim,degp1,
               inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
               inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
               "input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int cmplx2_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
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
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
   bool* zeroQ_h, bool *noqr_h, bool* zeroQ_d, bool *noqr_d,
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
            urhsrehi_h[i][j] = rhsrehi_h[i][j];
            urhsimhi_h[i][j] = rhsimhi_h[i][j];
            urhsrelo_h[i][j] = rhsrelo_h[i][j];
            urhsimlo_h[i][j] = rhsimlo_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solrehi_h[i][j] = 0.0; solimhi_h[i][j] = 0.0;
            solrelo_h[i][j] = 0.0; solimlo_h[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling CPU_cmplx2_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_cmplx2_qrbs_solve
         (dim,degp1,oldtail,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
          urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,
          zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_cmplx2_linear_residue ..." << endl;

         CPU_cmplx2_linear_residue
            (dim,degp1,*tailidx_h-1,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,solrehi_h,
             solrelo_h,solimhi_h,solimlo_h,
             resvecrehi,resvecrelo,resvecimhi,resvecimlo,
             resmaxhi,resmaxlo,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi << endl;
      }
      cmplx2_update_series
         (dim,degp1,*tailidx_h-1,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhsrehi_d[i][j] = rhsrehi_d[i][j];
            urhsimhi_d[i][j] = rhsimhi_d[i][j];
            urhsrelo_d[i][j] = rhsrelo_d[i][j];
            urhsimlo_d[i][j] = rhsimlo_d[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solrehi_d[i][j] = 0.0; solimhi_d[i][j] = 0.0;
            solrelo_d[i][j] = 0.0; solimlo_d[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling GPU_cmplx2_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_cmplx2_bals_solve
         (dim,degp1,szt,nbt,oldtail,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
          Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,zeroQ_d,noqr_d,
          upidx_d,bsidx_d,&newtail,totqrlapsedms,totqtblapsedms,
          totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_cmplx2_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_cmplx2_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
             solrehi_d,solrelo_d,solimhi_d,solimlo_d,
             resvecrehi,resvecrelo,resvecimhi,resvecimlo,
             resmaxhi,resmaxlo,&elapsedms,&addcnt,&mulcnt,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;

         *totreslapsedms += elapsedms;
      }
      cmplx2_update_series
         (dim,degp1,*tailidx_d-1,
          inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,vrblvl);
   }
   return 0;
}

int cmplx2_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehi, double **mbrelo, double **mbimhi, double **mbimlo,
   double dpr,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
   double **accrehi, double **accrelo, double **accimhi, double **accimlo,
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
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
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
         cout << "calling CPU_cmplx2_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx2_unit_series_vector
            (dim,deg,cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_cmplx2_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0],
             accrehi[0],accrelo[0],accimhi[0],accimlo[0],
             inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
             outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,vrblvl);
      }
      else
         CPU_cmplx2_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,
             cffrehi,cffrelo,cffimhi,cffimlo,
             accrehi,accrelo,accimhi,accimlo,
             inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
             funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_cmplx2_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx2_unit_series_vector
            (dim,deg,cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0]);
         // reset the coefficients
/*
         GPU_cmplx2_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0],
             accrehi[0],accrelo[0],accimhi[0],accimlo[0],
             inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             totcnvlapsedms,vrblvl);
 */
         GPU_cmplx2vectorized_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0],
             accrehi[0],accrelo[0],accimhi[0],accimlo[0],
             inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             totcnvlapsedms,vrblvl);
      }
      else
         GPU_cmplx2_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cffrehi,cffrelo,cffimhi,cffimlo,
             inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
             totcnvlapsedms,vrblvl);

   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;

      double errsum = cmplx2_error3sum
                         (dim,dim+1,degp1,
                          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
                          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
                          "output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum << endl;
   }
   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         cmplx2_define_rhs
            (dim,degp1,mbrehi,mbrelo,mbimhi,mbimlo,
             funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalrehi_h[i][j][k] = 0.0; jacvalimhi_h[i][j][k] = 0.0;
                  jacvalrelo_h[i][j][k] = 0.0; jacvalimlo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         cmplx2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbrehi,mbrelo,mbimhi,mbimlo,dpr,
             outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
             funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         cmplx2_define_rhs
            (dim,degp1,mbrehi,mbrelo,mbimhi,mbimlo,
             funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalrehi_d[i][j][k] = 0.0; jacvalimhi_d[i][j][k] = 0.0;
                  jacvalrelo_d[i][j][k] = 0.0; jacvalimlo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         cmplx2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbrehi,mbrelo,mbimhi,mbimlo,dpr,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = cmplx2_errors_funjacrhs(dim,deg,
                    funvalrehi_h,funvalrelo_h, funvalimhi_h,funvalimlo_h,
                    funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
                    jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
                    jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
                    rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
                    rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,vrblvl);
   }
   cmplx2_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
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
       workvecrehi,workvecrelo,workvecimhi,workvecimlo,
       resvecrehi,resvecrelo,resvecimhi,resvecimlo,resmaxhi,resmaxlo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return cmplx2_errors_inurhsQRsol(dim,deg,
                inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
                inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
                Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
                Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
                Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
                Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
                urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
                urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
                solrehi_h,solrelo_h,solimhi_h,solimlo_h,
                solrehi_d,solrelo_d,solimhi_d,solimlo_d,vrblvl);
   }
   else
      return 0;
}

int cmplx2_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cstrehi, double **cstrelo, double **cstimhi, double **cstimlo,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
   double dpr,
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
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
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
         cout << "calling CPU_cmplx2_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_cmplx2_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],
            cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
            cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
            inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
            outputrehi_h[i],outputrelo_h[i],
            outputimhi_h[i],outputimlo_h[i],&lapsed,0); // no output

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
         cout << "calling GPU_cmplx2_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ComplexConvolutionJobs cnvjobs(dim);
         ComplexIncrementJobs incjobs(cnvjobs,false);
         ComplexAdditionJobs addjobs(dim,nbr[i]);

         make_all_complex_jobs
            (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,false);

         GPU_cmplx2vectorized_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
             cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
             inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
             outputrehi_d[i],outputrelo_d[i],outputimhi_d[i],outputimlo_d[i],
             cnvjobs,incjobs,addjobs,
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

      errsum = cmplx2_error3sum
                  (dim,dim+1,degp1,
                   outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
                   outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
                   "output",vrblvl);
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
      cmplx2_map_evaldiff_output(dim,deg,
          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
          funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsrehi_h[i][j] = -funvalrehi_h[j][i];
            rhsrelo_h[i][j] = -funvalrelo_h[j][i];
            rhsimhi_h[i][j] = -funvalimhi_h[j][i];
            rhsimlo_h[i][j] = -funvalimlo_h[j][i];
         }
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx2_map_evaldiff_output(dim,deg,
          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
          funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsrehi_d[i][j] = -funvalrehi_d[j][i];
            rhsrelo_d[i][j] = -funvalrelo_d[j][i];
            rhsimhi_d[i][j] = -funvalimhi_d[j][i];
            rhsimlo_d[i][j] = -funvalimlo_d[j][i];
         }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = cmplx2_errors_funjacrhs(dim,deg,
                    funvalrehi_h,funvalrelo_h, funvalimhi_h,funvalimlo_h,
                    funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
                    jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
                    jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
                    rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
                    rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,vrblvl);
   }
   cmplx2_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
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
       workvecrehi,workvecrelo,workvecimhi,workvecimlo,
       resvecrehi,resvecrelo,resvecimhi,resvecimlo,resmaxhi,resmaxlo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return cmplx2_errors_inurhsQRsol(dim,deg,
                inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
                inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
                Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
                Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
                Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
                Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
                urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
                urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
                solrehi_h,solrelo_h,solimhi_h,solimlo_h,
                solrehi_d,solrelo_d,solimhi_d,solimlo_d,vrblvl);
   }
   else
      return 0;
}

int cmplx2_allocate_inoutfunjac
 ( int dim, int deg, int mode,
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
   double ***jacvalimhi_d, double ***jacvalimlo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputrehi_h[i] = new double[degp1];
         inputrelo_h[i] = new double[degp1];
         inputimhi_h[i] = new double[degp1];
         inputimlo_h[i] = new double[degp1];

         outputrehi_h[i] = new double*[dim+1];
         outputrelo_h[i] = new double*[dim+1];
         outputimhi_h[i] = new double*[dim+1];
         outputimlo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputrehi_h[i][j] = new double[degp1];
            outputrelo_h[i][j] = new double[degp1];
            outputimhi_h[i][j] = new double[degp1];
            outputimlo_h[i][j] = new double[degp1];
         }
         funvalrehi_h[i] = new double[degp1];
         funvalrelo_h[i] = new double[degp1];
         funvalimhi_h[i] = new double[degp1];
         funvalimlo_h[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalrehi_h[i] = new double*[dim];
         jacvalrelo_h[i] = new double*[dim];
         jacvalimhi_h[i] = new double*[dim];
         jacvalimlo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalrehi_h[i][j] = new double[dim];
            jacvalrelo_h[i][j] = new double[dim];
            jacvalimhi_h[i][j] = new double[dim];
            jacvalimlo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputrehi_d[i] = new double[degp1];
         inputrelo_d[i] = new double[degp1];
         inputimhi_d[i] = new double[degp1];
         inputimlo_d[i] = new double[degp1];

         outputrehi_d[i] = new double*[dim+1];
         outputrelo_d[i] = new double*[dim+1];
         outputimhi_d[i] = new double*[dim+1];
         outputimlo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputrehi_d[i][j] = new double[degp1];
            outputrelo_d[i][j] = new double[degp1];
            outputimhi_d[i][j] = new double[degp1];
            outputimlo_d[i][j] = new double[degp1];
         }
         funvalrehi_d[i] = new double[degp1];
         funvalrelo_d[i] = new double[degp1];
         funvalimhi_d[i] = new double[degp1];
         funvalimlo_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalrehi_d[i] = new double*[dim];
         jacvalrelo_d[i] = new double*[dim];
         jacvalimhi_d[i] = new double*[dim];
         jacvalimlo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalrehi_d[i][j] = new double[dim];
            jacvalrelo_d[i][j] = new double[dim];
            jacvalimhi_d[i][j] = new double[dim];
            jacvalimlo_d[i][j] = new double[dim];
         }
      }
   }
   return 0;
}

int cmplx2_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d,
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h, 
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d, 
   double **urhsimhi_d, double **urhsimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h, 
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d,
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h,
   double **solrehi_d, double **solrelo_d,
   double **solimhi_d, double **solimlo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) 
      {
         rhsrehi_h[i] = new double[dim];
         rhsrelo_h[i] = new double[dim];
         rhsimhi_h[i] = new double[dim];
         rhsimlo_h[i] = new double[dim];
         urhsrehi_h[i] = new double[dim];
         urhsrelo_h[i] = new double[dim];
         urhsimhi_h[i] = new double[dim];
         urhsimlo_h[i] = new double[dim];
         solrehi_h[i] = new double[dim];
         solrelo_h[i] = new double[dim];
         solimhi_h[i] = new double[dim];
         solimlo_h[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qrehi_h[i] = new double[dim];
         Qrelo_h[i] = new double[dim];
         Qimhi_h[i] = new double[dim];
         Qimlo_h[i] = new double[dim];
         Rrehi_h[i] = new double[dim];
         Rrelo_h[i] = new double[dim];
         Rimhi_h[i] = new double[dim];
         Rimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) 
      {
         rhsrehi_d[i] = new double[dim];
         rhsrelo_d[i] = new double[dim];
         rhsimhi_d[i] = new double[dim];
         rhsimlo_d[i] = new double[dim];
         urhsrehi_d[i] = new double[dim];
         urhsrelo_d[i] = new double[dim];
         urhsimhi_d[i] = new double[dim];
         urhsimlo_d[i] = new double[dim];
         solrehi_d[i] = new double[dim];
         solrelo_d[i] = new double[dim];
         solimhi_d[i] = new double[dim];
         solimlo_d[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qrehi_d[i] = new double[dim];
         Qrelo_d[i] = new double[dim];
         Qimhi_d[i] = new double[dim];
         Qimlo_d[i] = new double[dim];
         Rrehi_d[i] = new double[dim];
         Rrelo_d[i] = new double[dim];
         Rimhi_d[i] = new double[dim];
         Rimlo_d[i] = new double[dim];
      }
   }
   return 0;
}

void cmplx2_start_setup
 ( int dim, int deg,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d, int mode, int vrblvl )
{
   double *start0rehi = new double[dim];
   double *start0relo = new double[dim];
   double *start0imhi = new double[dim];
   double *start0imlo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0rehi[i] = testsolrehi[i][0];
      start0relo[i] = testsolrelo[i][0];
      start0imhi[i] = testsolimhi[i][0]; 
      start0imlo[i] = testsolimlo[i][0]; 
   }
   if((mode == 1) || (mode == 2))
      cmplx2_start_series_vector(dim,deg,
          start0rehi,start0relo,start0imhi,start0imlo,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h);
   else
      cmplx2_start_series_vector(dim,deg,
          start0rehi,start0relo,start0imhi,start0imlo,
          inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
         {
            inputrehi_d[i][j] = inputrehi_h[i][j];
            inputrelo_d[i][j] = inputrelo_h[i][j];
            inputimhi_d[i][j] = inputimhi_h[i][j];
            inputimlo_d[i][j] = inputimlo_h[i][j];
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
                 << inputrehi_h[i][0] << "  "
                 << inputrelo_h[i][0] << endl
                 << inputimhi_h[i][0] << "  "
                 << inputimlo_h[i][0] << endl;
      }
      else
      {
         for(int i=0; i<dim; i++)
            cout << i << " : "
                 << inputrehi_d[i][0] << "  "
                 << inputrelo_d[i][0] << endl
                 << inputimhi_d[i][0] << "  "
                 << inputimlo_d[i][0] << endl;
      }
   }
}

void cmplx2_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **mbrhsrehi, double **mbrhsrelo,
   double **mbrhsimhi, double **mbrhsimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolrehi[i] = new double[degp1];
      testsolrelo[i] = new double[degp1];
      testsolimhi[i] = new double[degp1];
      testsolimlo[i] = new double[degp1];
   }
   make_complex2_exponentials
      (dim,deg,testsolrehi,testsolrelo,testsolimhi,testsolimlo);

   // compute the right hand sides via evaluation

   for(int i=0; i<dim; i++)
   {
      mbrhsrehi[i] = new double[degp1];
      mbrhsrelo[i] = new double[degp1];
      mbrhsimhi[i] = new double[degp1];
      mbrhsimlo[i] = new double[degp1];

      mbrhsrehi[i][0] = 1.0;     // initialize product to one
      mbrhsrelo[i][0] = 0.0; 
      mbrhsimhi[i][0] = 0.0;
      mbrhsimlo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhsrehi[i][k] = 0.0; mbrhsrelo[i][k] = 0.0;
         mbrhsimhi[i][k] = 0.0; mbrhsimlo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_complex2_monomials
         (dim,deg,rowsA,testsolrehi,testsolrelo,testsolimhi,testsolimlo,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo);
   else
      evaluate_complex2_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffrehi,cffrelo,cffimhi,cffimlo,
          testsolrehi,testsolrelo,testsolimhi,testsolimlo,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,vrblvl);
 
   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsrehi[i][j] << "  "
                 << mbrhsrelo[i][j] << endl
                 << mbrhsimhi[i][j] << "  "
                 << mbrhsimlo[i][j] << endl;
   }
   cmplx2_start_setup(dim,deg,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,mode,vrblvl);
}

void cmplx2_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehi, double **cstrelo, double **cstimhi, double **cstimlo,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolrehi[i] = new double[degp1];
      testsolrelo[i] = new double[degp1];
      testsolimhi[i] = new double[degp1];
      testsolimlo[i] = new double[degp1];
   }
   make_complex2_exponentials(dim,deg,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_cmplx2_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],
            cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
            cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
            testsolrehi,testsolrelo,testsolimhi,testsolimlo,
            outputrehi_h[i],outputrelo_h[i],
            outputimhi_h[i],outputimlo_h[i],&lapsed,0); // no output

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
            // cstre[i][j] = cstre[i][j] - outputre_h[i][dim][j];
            ddf_dec(&cstrehi[i][j],&cstrelo[i][j],
                    outputrehi_h[i][dim][j],outputrelo_h[i][dim][j]);
            // cstim[i][j] = cstim[i][j] - outputim_h[i][dim][j];
            ddf_dec(&cstimhi[i][j],&cstimlo[i][j],
                    outputimhi_h[i][dim][j],outputimlo_h[i][dim][j]);
         }
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_cmplx2_poly_evaldiff
              (dim,nbr[i],deg,nvr[i],idx[i],
               cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
               cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
               testsolrehi,testsolrelo,testsolimhi,testsolimlo,
               outputrehi_h[i],outputrelo_h[i],
               outputimhi_h[i],outputimlo_h[i],&lapsed,0); // no output

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
               errsum = errsum + outputrehi_h[i][dim][j]
                               + outputrelo_h[i][dim][j]
                               + outputimhi_h[i][dim][j]
                               + outputimlo_h[i][dim][j]; 

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

         GPU_cmplx2vectorized_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
             cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
             testsolrehi,testsolrelo,testsolimhi,testsolimlo,
             outputrehi_d[i],outputrelo_d[i],outputimhi_d[i],outputimlo_d[i],
             cnvjobs,incjobs,addjobs,
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
            // cstre[i][j] = cstre[i][j] - outputre_d[i][dim][j];
            ddf_dec(&cstrehi[i][j],&cstrelo[i][j],
                    outputrehi_d[i][dim][j],outputrelo_d[i][dim][j]);
            // cstim[i][j] = cstim[i][j] - outputim_d[i][dim][j];
            ddf_dec(&cstimhi[i][j],&cstimlo[i][j],
                    outputimhi_d[i][dim][j],outputimlo_d[i][dim][j]);
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

            GPU_cmplx2vectorized_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],
                cstrehi[i],cstrelo[i],cstimhi[i],cstimlo[i],
                cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
                testsolrehi,testsolrelo,testsolimhi,testsolimlo,
                outputrehi_d[i],outputrelo_d[i],
                outputimhi_d[i],outputimlo_d[i],cnvjobs,incjobs,addjobs,
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
               errsum = errsum + outputrehi_d[i][dim][j]
                               + outputrelo_d[i][dim][j]
                               + outputimhi_d[i][dim][j]
                               + outputimlo_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
      }
   }
   cmplx2_start_setup(dim,deg,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,mode,vrblvl);
}

int cmplx2_error_testsol
 ( int dim, int deg, int mode,
   double **testsolrehi, double **testsolrelo,
   double **testsolimhi, double **testsolimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d )
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
                        << testsolrehi[i][j] << "  "
                        << testsolrelo[i][j] << endl << "  "
                        << testsolimhi[i][j] << "  "
                        << testsolimlo[i][j] << endl;
         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehi_d[i][j] << "  "
                           << inputrelo_d[i][j] << endl << "  "
                           << inputimhi_d[i][j] << "  "
                           << inputimlo_d[i][j] << endl;
            errsum += abs(testsolrehi[i][j] - inputrehi_d[i][j])
                    + abs(testsolrelo[i][j] - inputrelo_d[i][j])
                    + abs(testsolimhi[i][j] - inputimhi_d[i][j])
                    + abs(testsolimlo[i][j] - inputimlo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehi_h[i][j] << "  "
                           << inputrelo_h[i][j] << endl << "  "
                           << inputimhi_h[i][j] << "  "
                           << inputimlo_h[i][j] << endl;
            errsum += abs(testsolrehi[i][j] - inputrehi_h[i][j])
                    + abs(testsolrelo[i][j] - inputrelo_h[i][j])
                    + abs(testsolimhi[i][j] - inputimhi_h[i][j])
                    + abs(testsolimlo[i][j] - inputimlo_h[i][j]);
         }
      }
   }
   cout << scientific << setprecision(2)
        << "error : " << errsum << endl;

   return (errsum > 1.0e-20);
}

int test_cmplx2_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **accrehi = new double*[dim+1]; // accumulated power series
   double **accrelo = new double*[dim+1]; // in one column
   double **accimhi = new double*[dim+1];
   double **accimlo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehi[i] = new double[degp1];
      accrelo[i] = new double[degp1];
      accimhi[i] = new double[degp1];
      accimlo[i] = new double[degp1];
   }
   double ***cffrehi = new double**[nbrcol]; // coefficients of monomials
   double ***cffrelo = new double**[nbrcol];
   double ***cffimhi = new double**[nbrcol]; 
   double ***cffimlo = new double**[nbrcol]; 

   for(int i=0; i<nbrcol; i++)
   {
      cffrehi[i] = new double*[dim];
      cffrelo[i] = new double*[dim];
      cffimhi[i] = new double*[dim];
      cffimlo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffrehi[i][j] = new double[degp1];
         cffrelo[i][j] = new double[degp1];
         cffimhi[i][j] = new double[degp1];
         cffimlo[i][j] = new double[degp1];
      }
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_complex2_coefficients(nbrcol,dim,cffrehi,cffrelo,cffimhi,cffimlo);

   double **inputrehi_h;
   double **inputrelo_h;
   double **inputimhi_h;
   double **inputimlo_h;
   double **inputrehi_d;
   double **inputrelo_d;
   double **inputimhi_d;
   double **inputimlo_d;
   double ***outputrehi_h;
   double ***outputrelo_h;
   double ***outputimhi_h;
   double ***outputimlo_h;
   double ***outputrehi_d;
   double ***outputrelo_d;
   double ***outputimhi_d;
   double ***outputimlo_d;
   double **funvalrehi_h;
   double **funvalrelo_h;
   double **funvalimhi_h;
   double **funvalimlo_h;
   double **funvalrehi_d;
   double **funvalrelo_d;
   double **funvalimhi_d;
   double **funvalimlo_d;
   double ***jacvalrehi_h;
   double ***jacvalrelo_h;
   double ***jacvalimhi_h;
   double ***jacvalimlo_h;
   double ***jacvalrehi_d;
   double ***jacvalrelo_d;
   double ***jacvalimhi_d;
   double ***jacvalimlo_d;

   if((mode == 1) || (mode == 2))
   {
      inputrehi_h = new double*[dim];
      inputrelo_h = new double*[dim];
      inputimhi_h = new double*[dim];
      inputimlo_h = new double*[dim];
      outputrehi_h = new double**[dim];
      outputrelo_h = new double**[dim];
      outputimhi_h = new double**[dim];
      outputimlo_h = new double**[dim];
      funvalrehi_h = new double*[dim];
      funvalrelo_h = new double*[dim];
      funvalimhi_h = new double*[dim];
      funvalimlo_h = new double*[dim];
      jacvalrehi_h = new double**[degp1];
      jacvalrelo_h = new double**[degp1];
      jacvalimhi_h = new double**[degp1];
      jacvalimlo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputrehi_d = new double*[dim];
      inputrelo_d = new double*[dim];
      inputimhi_d = new double*[dim];
      inputimlo_d = new double*[dim];
      outputrehi_d = new double**[dim];
      outputrelo_d = new double**[dim];
      outputimhi_d = new double**[dim];
      outputimlo_d = new double**[dim];
      funvalrehi_d = new double*[dim];
      funvalrelo_d = new double*[dim];
      funvalimhi_d = new double*[dim];
      funvalimlo_d = new double*[dim];
      jacvalrehi_d = new double**[degp1];
      jacvalrelo_d = new double**[degp1];
      jacvalimhi_d = new double**[degp1];
      jacvalimlo_d = new double**[degp1];
   }
   cmplx2_allocate_inoutfunjac(dim,deg,mode,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
       outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
       outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
       funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
       funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
       jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
       jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhsrehi_h;
   double **rhsrelo_h;
   double **rhsimhi_h;
   double **rhsimlo_h;
   double **rhsrehi_d;
   double **rhsrelo_d;
   double **rhsimhi_d;
   double **rhsimlo_d;
   double **urhsrehi_h;
   double **urhsrelo_h;
   double **urhsimhi_h;
   double **urhsimlo_h;
   double **urhsrehi_d;
   double **urhsrelo_d;
   double **urhsimhi_d;
   double **urhsimlo_d;
   double **solrehi_h;
   double **solrelo_h;
   double **solimhi_h;
   double **solimlo_h;
   double **solrehi_d;
   double **solrelo_d;
   double **solimhi_d;
   double **solimlo_d;
   double **Qrehi_h;
   double **Qrelo_h;
   double **Qimhi_h;
   double **Qimlo_h;
   double **Qrehi_d;
   double **Qrelo_d;
   double **Qimhi_d;
   double **Qimlo_d;
   double **Rrehi_h;
   double **Rrelo_h;
   double **Rimhi_h;
   double **Rimlo_h;
   double **Rrehi_d;
   double **Rrelo_d;
   double **Rimhi_d;
   double **Rimlo_d;

   if((mode == 1) || (mode == 2))
   {
      urhsrehi_h = new double*[degp1];
      urhsrelo_h = new double*[degp1];
      urhsimhi_h = new double*[degp1];
      urhsimlo_h = new double*[degp1];
      rhsrehi_h = new double*[degp1];
      rhsrelo_h = new double*[degp1];
      rhsimhi_h = new double*[degp1];
      rhsimlo_h = new double*[degp1];
      solrehi_h = new double*[degp1];
      solrelo_h = new double*[degp1];
      solimhi_h = new double*[degp1];
      solimlo_h = new double*[degp1];
      Qrehi_h = new double*[dim];
      Qrelo_h = new double*[dim];
      Qimhi_h = new double*[dim];
      Qimlo_h = new double*[dim];
      Rrehi_h = new double*[dim];
      Rrelo_h = new double*[dim];
      Rimhi_h = new double*[dim];
      Rimlo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      urhsrehi_d = new double*[degp1];
      urhsrelo_d = new double*[degp1];
      urhsimhi_d = new double*[degp1];
      urhsimlo_d = new double*[degp1];
      rhsrehi_d = new double*[degp1];
      rhsrelo_d = new double*[degp1];
      rhsimhi_d = new double*[degp1];
      rhsimlo_d = new double*[degp1];
      solrehi_d = new double*[degp1];
      solrelo_d = new double*[degp1];
      solimhi_d = new double*[degp1];
      solimlo_d = new double*[degp1];
      Qrehi_d = new double*[dim];
      Qrelo_d = new double*[dim];
      Qimhi_d = new double*[dim];
      Qimlo_d = new double*[dim];
      Rrehi_d = new double*[dim];
      Rrelo_d = new double*[dim];
      Rimhi_d = new double*[dim];
      Rimlo_d = new double*[dim];
   }
   cmplx2_allocate_rhsqrsol(dim,deg,mode,
       rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
       rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
       urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
       urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
       Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
       Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
       Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
       Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
       solrehi_h,solrelo_h,solimhi_h,solimlo_h,
       solrehi_d,solrelo_d,solimhi_d,solimlo_d);
   
   double *workvecrehi = new double[dim];
   double *workvecrelo = new double[dim];
   double *workvecimhi = new double[dim];
   double *workvecimlo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

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
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the test solution and the start series.

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolrehi = new double*[dim];
   double **testsolrelo = new double*[dim];
   double **testsolimhi = new double*[dim];
   double **testsolimlo = new double*[dim];
   double **mbrhsrehi = new double*[dim];
   double **mbrhsrelo = new double*[dim];
   double **mbrhsimhi = new double*[dim];
   double **mbrhsimlo = new double*[dim];

   cmplx2_column_setup
      (dim,deg,nbrcol,nvr,idx,rowsA,cffrehi,cffrelo,cffimhi,cffimlo,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo,
       mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,mode,vrblvl);

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

      cmplx2_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,dpr,
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
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,
          resvecrehi,resvecrelo,resvecimhi,resvecimlo,&resmaxhi,&resmaxlo,
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

   cmplx2_error_testsol
      (dim,deg,mode,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_cmplx2_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   double **cstrehi = new double*[dim];
   double **cstrelo = new double*[dim];
   double **cstimhi = new double*[dim]; 
   double **cstimlo = new double*[dim]; 
   double ***cffrehi = new double**[dim];
   double ***cffrelo = new double**[dim];
   double ***cffimhi = new double**[dim]; 
   double ***cffimlo = new double**[dim]; 

   cmplx2_make_coefficients(dim,deg,nbr,nvr,idx,
      cstrehi,cstrelo,cstimhi,cstimlo,
      cffrehi,cffrelo,cffimhi,cffimlo,vrblvl);

   double **inputrehi_h;
   double **inputrelo_h;
   double **inputimhi_h;
   double **inputimlo_h;
   double **inputrehi_d;
   double **inputrelo_d;
   double **inputimhi_d;
   double **inputimlo_d;
   double ***outputrehi_h;
   double ***outputrelo_h;
   double ***outputimhi_h;
   double ***outputimlo_h;
   double ***outputrehi_d;
   double ***outputrelo_d;
   double ***outputimhi_d;
   double ***outputimlo_d;
   double **funvalrehi_h;
   double **funvalrelo_h;
   double **funvalimhi_h;
   double **funvalimlo_h;
   double **funvalrehi_d;
   double **funvalrelo_d;
   double **funvalimhi_d;
   double **funvalimlo_d;
   double ***jacvalrehi_h;
   double ***jacvalrelo_h;
   double ***jacvalimhi_h;
   double ***jacvalimlo_h;
   double ***jacvalrehi_d;
   double ***jacvalrelo_d;
   double ***jacvalimhi_d;
   double ***jacvalimlo_d;

   if((mode == 1) || (mode == 2))
   {
      inputrehi_h = new double*[dim];
      inputrelo_h = new double*[dim];
      inputimhi_h = new double*[dim];
      inputimlo_h = new double*[dim];
      outputrehi_h = new double**[dim];
      outputrelo_h = new double**[dim];
      outputimhi_h = new double**[dim];
      outputimlo_h = new double**[dim];
      funvalrehi_h = new double*[dim];
      funvalrelo_h = new double*[dim];
      funvalimhi_h = new double*[dim];
      funvalimlo_h = new double*[dim];
      jacvalrehi_h = new double**[degp1];
      jacvalrelo_h = new double**[degp1];
      jacvalimhi_h = new double**[degp1];
      jacvalimlo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputrehi_d = new double*[dim];
      inputrelo_d = new double*[dim];
      inputimhi_d = new double*[dim];
      inputimlo_d = new double*[dim];
      outputrehi_d = new double**[dim];
      outputrelo_d = new double**[dim];
      outputimhi_d = new double**[dim];
      outputimlo_d = new double**[dim];
      funvalrehi_d = new double*[dim];
      funvalrelo_d = new double*[dim];
      funvalimhi_d = new double*[dim];
      funvalimlo_d = new double*[dim];
      jacvalrehi_d = new double**[degp1];
      jacvalrelo_d = new double**[degp1];
      jacvalimhi_d = new double**[degp1];
      jacvalimlo_d = new double**[degp1];
   }
   cmplx2_allocate_inoutfunjac(dim,deg,mode,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
       outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
       outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
       funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
       funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
       jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
       jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d);

   double **rhsrehi_h;
   double **rhsrelo_h;
   double **rhsimhi_h;
   double **rhsimlo_h;
   double **rhsrehi_d;
   double **rhsrelo_d;
   double **rhsimhi_d;
   double **rhsimlo_d;
   double **urhsrehi_h;
   double **urhsrelo_h;
   double **urhsimhi_h;
   double **urhsimlo_h;
   double **urhsrehi_d;
   double **urhsrelo_d;
   double **urhsimhi_d;
   double **urhsimlo_d;
   double **solrehi_h;
   double **solrelo_h;
   double **solimhi_h;
   double **solimlo_h;
   double **solrehi_d;
   double **solrelo_d;
   double **solimhi_d;
   double **solimlo_d;
   double **Qrehi_h;
   double **Qrelo_h;
   double **Qimhi_h;
   double **Qimlo_h;
   double **Qrehi_d;
   double **Qrelo_d;
   double **Qimhi_d;
   double **Qimlo_d;
   double **Rrehi_h;
   double **Rrelo_h;
   double **Rimhi_h;
   double **Rimlo_h;
   double **Rrehi_d;
   double **Rrelo_d;
   double **Rimhi_d;
   double **Rimlo_d;

   if((mode == 1) || (mode == 2))
   {
      urhsrehi_h = new double*[degp1];
      urhsrelo_h = new double*[degp1];
      urhsimhi_h = new double*[degp1];
      urhsimlo_h = new double*[degp1];
      rhsrehi_h = new double*[degp1];
      rhsrelo_h = new double*[degp1];
      rhsimhi_h = new double*[degp1];
      rhsimlo_h = new double*[degp1];
      solrehi_h = new double*[degp1];
      solrelo_h = new double*[degp1];
      solimhi_h = new double*[degp1];
      solimlo_h = new double*[degp1];
      Qrehi_h = new double*[dim];
      Qrelo_h = new double*[dim];
      Qimhi_h = new double*[dim];
      Qimlo_h = new double*[dim];
      Rrehi_h = new double*[dim];
      Rrelo_h = new double*[dim];
      Rimhi_h = new double*[dim];
      Rimlo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      urhsrehi_d = new double*[degp1];
      urhsrelo_d = new double*[degp1];
      urhsimhi_d = new double*[degp1];
      urhsimlo_d = new double*[degp1];
      rhsrehi_d = new double*[degp1];
      rhsrelo_d = new double*[degp1];
      rhsimhi_d = new double*[degp1];
      rhsimlo_d = new double*[degp1];
      solrehi_d = new double*[degp1];
      solrelo_d = new double*[degp1];
      solimhi_d = new double*[degp1];
      solimlo_d = new double*[degp1];
      Qrehi_d = new double*[dim];
      Qrelo_d = new double*[dim];
      Qimhi_d = new double*[dim];
      Qimlo_d = new double*[dim];
      Rrehi_d = new double*[dim];
      Rrelo_d = new double*[dim];
      Rimhi_d = new double*[dim];
      Rimlo_d = new double*[dim];
   }
   cmplx2_allocate_rhsqrsol(dim,deg,mode,
       rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
       rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
       urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
       urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
       Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
       Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
       Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
       Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
       solrehi_h,solrelo_h,solimhi_h,solimlo_h,
       solrehi_d,solrelo_d,solimhi_d,solimlo_d);
   
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolrehi = new double*[dim];
   double **testsolrelo = new double*[dim];
   double **testsolimhi = new double*[dim];
   double **testsolimlo = new double*[dim];

   cmplx2_row_setup(dim,deg,nbr,nvr,idx,
       cstrehi,cstrelo,cstimhi,cstimlo,
       cffrehi,cffrelo,cffimhi,cffimlo,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
       outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
       outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,mode,vrblvl);

   double *workvecrehi = new double[dim];
   double *workvecrelo = new double[dim];
   double *workvecimhi = new double[dim];
   double *workvecimlo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

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

      cmplx2_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,
          nbr,nvr,idx,cstrehi,cstrelo,cstimhi,cstimlo,
          cffrehi,cffrelo,cffimhi,cffimlo,dpr,
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
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,
          resvecrehi,resvecrelo,resvecimhi,resvecimlo,&resmaxhi,&resmaxlo,
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

   cmplx2_error_testsol
      (dim,deg,mode,
       testsolrehi,testsolrelo,testsolimhi,testsolimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
