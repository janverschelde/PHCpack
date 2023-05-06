// The file cmplx2_path_tracker.cpp defines the functions with prototypes in
// the file cmplx2_path_tracker.h.

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
#include "double_double_functions.h"
#include "dbl2_factorizations.h"
#include "dbl2_monomial_systems.h"
#include "unimodular_matrices.h"
#include "dbl2_systems_host.h"
#include "dbl2_systems_kernels.h"
#include "dbl2_bals_host.h"
#include "dbl2_tail_kernels.h"
#include "dbl2_bals_kernels.h"
#include "dbl2_newton_testers.h"
#include "cmplx2_newton_method.h"
#include "cmplx2_path_tracker.h"

int cmplx2_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrhsrehi, double **mbrhsrelo,
   double **mbrhsimhi, double **mbrhsimlo,
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
   int vrblvl, int mode )
{
   const int degp1 = deg+1;
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

   double dpr = 1.0; // dummy value

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      cmplx2_newton_qrstep
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
          resvecrehi,resvecrelo,resvecimhi,resvecimlo,resmaxhi,resmaxlo,
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

   double errsum = 0.0;

   cout << scientific << setprecision(16); // just in case vrblvl == 0
   cout << "The solution series : " << endl;
   for(int j=0; j<degp1; j++)
   {
      cout << "coefficient of degree " << j << " :" << endl;
      for(int i=0; i<dim; i++)
      {
         if(mode == 0)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehi_d[i][j] << "  "
                           << inputrelo_d[i][j] << endl << "  "
                           << inputimhi_d[i][j] << "  "
                           << inputimlo_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehi_h[i][j] << "  "
                           << inputrelo_h[i][j] << endl << "  "
                           << inputimhi_h[i][j] << "  "
                           << inputimlo_h[i][j] << endl;
         }
         if(mode == 2)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehi_d[i][j] << "  "
                           << inputrelo_d[i][j] << endl << "  "
                           << inputimhi_d[i][j] << "  "
                           << inputimlo_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehi_h[i][j] << "  "
                           << inputrelo_h[i][j] << endl << "  "
                           << inputimhi_h[i][j] << "  "
                           << inputimlo_h[i][j] << endl;
            errsum += abs(inputrehi_h[i][j] - inputrehi_d[i][j])
                    + abs(inputrelo_h[i][j] - inputrelo_d[i][j])
                    + abs(inputimhi_h[i][j] - inputimhi_d[i][j])
                    + abs(inputimlo_h[i][j] - inputimlo_d[i][j]);
         }
      }
   }
   if(mode == 2) cout << "error : " << errsum << endl;

   cout << "Wall clock time on all " << stepcnt << " Newton steps : ";
   cout << fixed << setprecision(3) 
        << walltimesec << " seconds." << endl;
   cout << "     Time spent by all convolution kernels : "
        << totcnvlapsedms << " milliseconds." << endl;
   cout << "  Time spent by all Householder QR kernels : "
        << totqrlapsedms << " milliseconds." << endl;
   cout << "     Time spent by all Q times rhs kernels : "
        << totqtblapsedms << " milliseconds." << endl;
   cout << "Time spent by all backsubstitution kernels : "
        << totbslapsedms << " milliseconds." << endl;
   cout << "          Time spent by all update kernels : "
        << totupdlapsedms << " milliseconds." << endl;
   cout << "        Time spent by all residual kernels : "
        << totreslapsedms << " milliseconds." << endl;

   double totkerneltime = totcnvlapsedms + totqrlapsedms + totqtblapsedms
                        + totbslapsedms + totupdlapsedms + totreslapsedms;

   cout << "           Total time spent by all kernels : "
        << totkerneltime << " milliseconds." << endl;

   return 0;
}

int test_dbl2_cmplx_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
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
   double **accrehi = new double*[dim+1]; // accumulate series in one column
   double **accrelo = new double*[dim+1];
   double **accimhi = new double*[dim+1];
   double **accimlo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehi[i] = new double[degp1];
      accrelo[i] = new double[degp1];
      accimhi[i] = new double[degp1];
      accimlo[i] = new double[degp1];
   }
   double ***cffrehi = new double**[nbrcol];
   double ***cffrelo = new double**[nbrcol];
   double ***cffimhi = new double**[nbrcol];
   double ***cffimlo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffrehi[i] = new double*[dim]; // the coefficients of monomials
      cffrelo[i] = new double*[dim];
      cffimhi[i] = new double*[dim];
      cffimlo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffrehi[i][j] = new double[degp1];
         cffrelo[i][j] = new double[degp1];
         cffimhi[i][j] = new double[degp1];
         cffimlo[i][j] = new double[degp1];

         for(int k=0; k<deg; k++)
         {
            cffrehi[i][j][k] = 0.0;
            cffrelo[i][j][k] = 0.0;
            cffimhi[i][j][k] = 0.0;
            cffimlo[i][j][k] = 0.0;
         }
      }
   }
   double ***outputrehi_h;
   double ***outputrelo_h;
   double ***outputimhi_h;
   double ***outputimlo_h;
   double ***outputrehi_d;
   double ***outputrelo_d;
   double ***outputimhi_d;
   double ***outputimlo_d;

   if((mode == 1) || (mode == 2))
   {
      outputrehi_h = new double**[dim];
      outputrelo_h = new double**[dim];
      outputimhi_h = new double**[dim];
      outputimlo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
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
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputrehi_d = new double**[dim];
      outputrelo_d = new double**[dim];
      outputimhi_d = new double**[dim];
      outputimlo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
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
      }
   }
   double **funvalrehi_h;  // function values on host
   double **funvalrelo_h;
   double **funvalimhi_h;
   double **funvalimlo_h;
   double **funvalrehi_d;  // function values on device
   double **funvalrelo_d;
   double **funvalimhi_d;
   double **funvalimlo_d;

   if((mode == 1) || (mode == 2))
   {
      funvalrehi_h = new double*[dim];
      funvalrelo_h = new double*[dim];
      funvalimhi_h = new double*[dim];
      funvalimlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalrehi_h[i] = new double[degp1];
         funvalrelo_h[i] = new double[degp1];
         funvalimhi_h[i] = new double[degp1];
         funvalimlo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalrehi_d = new double*[dim];
      funvalrelo_d = new double*[dim];
      funvalimhi_d = new double*[dim];
      funvalimlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalrehi_d[i] = new double[degp1];
         funvalrelo_d[i] = new double[degp1];
         funvalimhi_d[i] = new double[degp1];
         funvalimlo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
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
      jacvalrehi_h = new double**[degp1];
      jacvalrelo_h = new double**[degp1];
      jacvalimhi_h = new double**[degp1];
      jacvalimlo_h = new double**[degp1];

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
      jacvalrehi_d = new double**[degp1];
      jacvalrelo_d = new double**[degp1];
      jacvalimhi_d = new double**[degp1];
      jacvalimlo_d = new double**[degp1];

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
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solrehi_h;
   double **solrelo_h;
   double **solimhi_h;
   double **solimlo_h;
   double **solrehi_d;
   double **solrelo_d;
   double **solimhi_d;
   double **solimlo_d;

   if((mode == 1) || (mode == 2))
   {
      solrehi_h = new double*[degp1];
      solrelo_h = new double*[degp1];
      solimhi_h = new double*[degp1];
      solimlo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solrehi_h[i] = new double[dim];
         solrelo_h[i] = new double[dim];
         solimhi_h[i] = new double[dim];
         solimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solrehi_d = new double*[degp1];
      solrelo_d = new double*[degp1];
      solimhi_d = new double*[degp1];
      solimlo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solrehi_d[i] = new double[dim];
         solrelo_d[i] = new double[dim];
         solimhi_d[i] = new double[dim];
         solimlo_d[i] = new double[dim];
      }
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

   if((mode == 1) || (mode == 2))
   {
      rhsrehi_h = new double*[degp1];
      rhsrelo_h = new double*[degp1];
      rhsimhi_h = new double*[degp1];
      rhsimlo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsrehi_h[i] = new double[dim];
         rhsrelo_h[i] = new double[dim];
         rhsimhi_h[i] = new double[dim];
         rhsimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhsrehi_d = new double*[degp1];
      rhsrelo_d = new double*[degp1];
      rhsimhi_d = new double*[degp1];
      rhsimlo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsrehi_d[i] = new double[dim];
         rhsrelo_d[i] = new double[dim];
         rhsimhi_d[i] = new double[dim];
         rhsimlo_d[i] = new double[dim];
      }
   }
   double *workvecrehi = new double[dim];
   double *workvecrelo = new double[dim];
   double *workvecimhi = new double[dim];
   double *workvecimlo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **urhsrehi_h;
   double **urhsrelo_h;
   double **urhsimhi_h;
   double **urhsimlo_h;
   double **urhsrehi_d;
   double **urhsrelo_d;
   double **urhsimhi_d;
   double **urhsimlo_d;

   if((mode == 1) || (mode == 2))
   {
      urhsrehi_h = new double*[degp1];
      urhsrelo_h = new double*[degp1];
      urhsimhi_h = new double*[degp1];
      urhsimlo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsrehi_h[i] = new double[dim];
         urhsrelo_h[i] = new double[dim];
         urhsimhi_h[i] = new double[dim];
         urhsimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhsrehi_d = new double*[degp1];
      urhsrelo_d = new double*[degp1];
      urhsimhi_d = new double*[degp1];
      urhsimlo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsrehi_d[i] = new double[dim];
         urhsrelo_d[i] = new double[dim];
         urhsimhi_d[i] = new double[dim];
         urhsimlo_d[i] = new double[dim];
      }
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
   double resmaxhi;
   double resmaxlo;

   double **Qrehi_h;
   double **Qrelo_h;
   double **Qimhi_h;
   double **Qimlo_h;
   double **Qrehi_d;
   double **Qrelo_d;
   double **Qimhi_d;
   double **Qimlo_d;

   if((mode == 1) || (mode == 2))
   {
      Qrehi_h = new double*[dim];
      Qrelo_h = new double*[dim];
      Qimhi_h = new double*[dim];
      Qimlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qrehi_h[i] = new double[dim];
         Qrelo_h[i] = new double[dim];
         Qimhi_h[i] = new double[dim];
         Qimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qrehi_d = new double*[dim];
      Qrelo_d = new double*[dim];
      Qimhi_d = new double*[dim];
      Qimlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qrehi_d[i] = new double[dim];
         Qrelo_d[i] = new double[dim];
         Qimhi_d[i] = new double[dim];
         Qimlo_d[i] = new double[dim];
      }
   }
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
      Rrehi_h = new double*[dim];
      Rrelo_h = new double*[dim];
      Rimhi_h = new double*[dim];
      Rimlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rrehi_h[i] = new double[dim];
         Rrelo_h[i] = new double[dim];
         Rimhi_h[i] = new double[dim];
         Rimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rrehi_d = new double*[dim];
      Rrelo_d = new double*[dim];
      Rimhi_d = new double*[dim];
      Rimlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rrehi_d[i] = new double[dim];
         Rrelo_d[i] = new double[dim];
         Rimhi_d[i] = new double[dim];
         Rimlo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **startsolrehi = new double*[dim];
   double **startsolrelo = new double*[dim];
   double **startsolimhi = new double*[dim];
   double **startsolimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      startsolrehi[i] = new double[degp1];
      startsolrelo[i] = new double[degp1];
      startsolimhi[i] = new double[degp1];
      startsolimlo[i] = new double[degp1];
   }
   make_complex2_exponentials
      (dim,deg,startsolrehi,startsolrelo,startsolimhi,startsolimlo);

   if(nbrcol != 1) // generate coefficients for the columns
   {
      // sets only the leading coefficient to a random double ...
      make_complex2_coefficients(nbrcol,dim,cffrehi,cffrelo,cffimhi,cffimlo);
   }
   // compute the right hand sides via evaluation
   double **mbrhsrehi = new double*[dim];
   double **mbrhsrelo = new double*[dim];
   double **mbrhsimhi = new double*[dim];
   double **mbrhsimlo = new double*[dim];

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
         mbrhsrehi[i][k] = 0.0;
         mbrhsrelo[i][k] = 0.0;
         mbrhsimhi[i][k] = 0.0;
         mbrhsimlo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_complex2_monomials
         (dim,deg,rowsA,startsolrehi,startsolrelo,startsolimhi,startsolimlo,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo);
   else
      evaluate_complex2_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffrehi,cffrelo,cffimhi,cffimlo,
          startsolrehi,startsolrelo,startsolimhi,startsolimlo,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++) // mbrhs[i][j] = mbrhs[i][j] - mbrhs[i][j-1];
      {
         double acchi,acclo;

         ddf_sub(mbrhsrehi[i][j],  mbrhsrelo[i][j],
                 mbrhsrehi[i][j-1],mbrhsrelo[i][j-1],&acchi,&acclo);
         mbrhsrehi[i][j] = acchi;
         mbrhsrelo[i][j] = acchi;
         ddf_sub(mbrhsimhi[i][j],  mbrhsimlo[i][j],
                 mbrhsimhi[i][j-1],mbrhsimlo[i][j-1],&acchi,&acclo);
         mbrhsimhi[i][j] = acchi;
         mbrhsimlo[i][j] = acchi;
      }

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsrehi[i][j] << "  " << mbrhsrelo[i][j] << endl
                 << mbrhsimhi[i][j] << "  " << mbrhsimlo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehi_h[i][j] = startsolrehi[i][j];
         inputrelo_h[i][j] = startsolrelo[i][j];
         inputimhi_h[i][j] = startsolimhi[i][j];
         inputimlo_h[i][j] = startsolimlo[i][j];
         inputrehi_d[i][j] = startsolrehi[i][j];
         inputrelo_d[i][j] = startsolrelo[i][j];
         inputimhi_d[i][j] = startsolimhi[i][j];
         inputimlo_d[i][j] = startsolimlo[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputrehi_h[i][0] << "  "
                            << inputrelo_h[i][0] << endl
                            << inputimhi_h[i][0] << "  "
                            << inputimlo_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   cmplx2_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,
       nvr,idx,exp,nbrfac,expfac,
       mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,
       cffrehi,cffrelo,cffimhi,cffimlo,
       accrehi,accrelo,accimhi,accimlo,
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
       vrblvl,mode);

   return 0;
}
