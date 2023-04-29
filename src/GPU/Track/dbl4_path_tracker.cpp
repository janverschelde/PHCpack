// The file dbl4_path_tracker.cpp defines the functions with prototypes in
// the file dbl4_path_tracker.h.

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
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"
#include "dbl4_monomial_systems.h"
#include "unimodular_matrices.h"
#include "dbl4_systems_host.h"
#include "dbl4_systems_kernels.h"
#include "dbl4_bals_host.h"
#include "dbl4_tail_kernels.h"
#include "dbl4_bals_kernels.h"
#include "dbl4_newton_testers.h"
#include "dbl4_newton_method.h"
#include "dbl_fabry_host.h"
#include "dbl4_path_tracker.h"

using namespace std;

int dbl4_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double **acchihi, double **acclohi, double **acchilo, double **acclolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h,
   double ***outputhilo_h, double ***outputlolo_h,
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double **workmathihi, double **workmatlohi,
   double **workmathilo, double **workmatlolo,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
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

      dbl4_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,
          &tailidx_h,&tailidx_d,nvr,idx,exp,nbrfac,expfac,
          rhshihi,rhslohi,rhshilo,rhslolo,dpr,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
          inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
          outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
          funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
          funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
          jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
          jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
          rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
          rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
          urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
          urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
          solhihi_h,sollohi_h,solhilo_h,sollolo_h,
          solhihi_d,sollohi_d,solhilo_d,sollolo_d,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
          Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          workvechihi,workveclohi,workvechilo,workveclolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
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
                           << inputhihi_d[i][j] << "  "
                           << inputlohi_d[i][j] << endl
                           << inputhilo_d[i][j] << "  "
                           << inputlolo_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhihi_h[i][j] << "  "
                           << inputlohi_h[i][j] << endl
                           << inputhilo_h[i][j] << "  "
                           << inputlolo_h[i][j] << endl;
         }
         if(mode == 2)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhihi_d[i][j] << "  "
                           << inputlohi_d[i][j] << endl
                           << inputhilo_d[i][j] << "  "
                           << inputlolo_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhihi_h[i][j] << "  "
                           << inputlohi_h[i][j] << endl
                           << inputhilo_h[i][j] << "  "
                           << inputlolo_h[i][j] << endl;
            errsum += abs(inputhihi_h[i][j] - inputhihi_d[i][j])
                    + abs(inputlohi_h[i][j] - inputlohi_d[i][j])
                    + abs(inputhilo_h[i][j] - inputhilo_d[i][j])
                    + abs(inputlolo_h[i][j] - inputlolo_d[i][j]);
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

int test_dbl4_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhihi_h = new double*[dim];
   double **inputlohi_h = new double*[dim];
   double **inputhilo_h = new double*[dim];
   double **inputlolo_h = new double*[dim];
   double **inputhihi_d = new double*[dim];
   double **inputlohi_d = new double*[dim];
   double **inputhilo_d = new double*[dim];
   double **inputlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihi_h[i] = new double[degp1];
      inputlohi_h[i] = new double[degp1];
      inputhilo_h[i] = new double[degp1];
      inputlolo_h[i] = new double[degp1];
      inputhihi_d[i] = new double[degp1];
      inputlohi_d[i] = new double[degp1];
      inputhilo_d[i] = new double[degp1];
      inputlolo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **acchihi = new double*[dim+1]; // accumulated power series
   double **acclohi = new double*[dim+1]; // in one column
   double **acchilo = new double*[dim+1];
   double **acclolo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchihi[i] = new double[degp1];
      acclohi[i] = new double[degp1];
      acchilo[i] = new double[degp1];
      acclolo[i] = new double[degp1];
   }

   double ***cffhihi = new double**[nbrcol]; // coefficients of monomials
   double ***cfflohi = new double**[nbrcol];
   double ***cffhilo = new double**[nbrcol];
   double ***cfflolo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffhihi[i] = new double*[dim];
      cfflohi[i] = new double*[dim];
      cffhilo[i] = new double*[dim];
      cfflolo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffhihi[i][j] = new double[degp1];
         cfflohi[i][j] = new double[degp1];
         cffhilo[i][j] = new double[degp1];
         cfflolo[i][j] = new double[degp1];

         for(int k=0; k<deg; k++)
         {
            cffhihi[i][j][k] = 0.0;
            cfflohi[i][j][k] = 0.0;
            cffhilo[i][j][k] = 0.0;
            cfflolo[i][j][k] = 0.0;
         }
      }
   }
   double ***outputhihi_h;
   double ***outputlohi_h;
   double ***outputhilo_h;
   double ***outputlolo_h;
   double ***outputhihi_d;
   double ***outputlohi_d;
   double ***outputhilo_d;
   double ***outputlolo_d;

   if((mode == 1) || (mode == 2))
   {
      outputhihi_h = new double**[dim];
      outputlohi_h = new double**[dim];
      outputhilo_h = new double**[dim];
      outputlolo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhihi_h[i] = new double*[dim+1];
         outputlohi_h[i] = new double*[dim+1];
         outputhilo_h[i] = new double*[dim+1];
         outputlolo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhihi_h[i][j] = new double[degp1];
            outputlohi_h[i][j] = new double[degp1];
            outputhilo_h[i][j] = new double[degp1];
            outputlolo_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputhihi_d = new double**[dim];
      outputlohi_d = new double**[dim];
      outputhilo_d = new double**[dim];
      outputlolo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhihi_d[i] = new double*[dim+1];
         outputlohi_d[i] = new double*[dim+1];
         outputhilo_d[i] = new double*[dim+1];
         outputlolo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhihi_d[i][j] = new double[degp1];
            outputlohi_d[i][j] = new double[degp1];
            outputhilo_d[i][j] = new double[degp1];
            outputlolo_d[i][j] = new double[degp1];
         }
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihi_h;
   double **funvallohi_h;
   double **funvalhilo_h;
   double **funvallolo_h;
   double **funvalhihi_d;
   double **funvallohi_d;
   double **funvalhilo_d;
   double **funvallolo_d;

   if((mode == 1) || (mode == 2))
   {
      funvalhihi_h = new double*[dim];
      funvallohi_h = new double*[dim];
      funvalhilo_h = new double*[dim];
      funvallolo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhihi_h[i] = new double[degp1];
         funvallohi_h[i] = new double[degp1];
         funvalhilo_h[i] = new double[degp1];
         funvallolo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalhihi_d = new double*[dim];
      funvallohi_d = new double*[dim];
      funvalhilo_d = new double*[dim];
      funvallolo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhihi_d[i] = new double[degp1];
         funvallohi_d[i] = new double[degp1];
         funvalhilo_d[i] = new double[degp1];
         funvallolo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihi_h;
   double ***jacvallohi_h;
   double ***jacvalhilo_h;
   double ***jacvallolo_h;
   double ***jacvalhihi_d;
   double ***jacvallohi_d;
   double ***jacvalhilo_d;
   double ***jacvallolo_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalhihi_h = new double**[degp1];
      jacvallohi_h = new double**[degp1];
      jacvalhilo_h = new double**[degp1];
      jacvallolo_h = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhihi_h[i] = new double*[dim];
         jacvallohi_h[i] = new double*[dim];
         jacvalhilo_h[i] = new double*[dim];
         jacvallolo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhihi_h[i][j] = new double[dim];
            jacvallohi_h[i][j] = new double[dim];
            jacvalhilo_h[i][j] = new double[dim];
            jacvallolo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      jacvalhihi_d = new double**[degp1];
      jacvallohi_d = new double**[degp1];
      jacvalhilo_d = new double**[degp1];
      jacvallolo_d = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhihi_d[i] = new double*[dim];
         jacvallohi_d[i] = new double*[dim];
         jacvalhilo_d[i] = new double*[dim];
         jacvallolo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhihi_d[i][j] = new double[dim];
            jacvallohi_d[i][j] = new double[dim];
            jacvalhilo_d[i][j] = new double[dim];
            jacvallolo_d[i][j] = new double[dim];
         }
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihi_h;
   double **sollohi_h;
   double **solhilo_h;
   double **sollolo_h;
   double **solhihi_d;
   double **sollohi_d;
   double **solhilo_d;
   double **sollolo_d;

   if((mode == 1) || (mode == 2))
   {
      solhihi_h = new double*[degp1];
      sollohi_h = new double*[degp1];
      solhilo_h = new double*[degp1];
      sollolo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhihi_h[i] = new double[dim];
         sollohi_h[i] = new double[dim];
         solhilo_h[i] = new double[dim];
         sollolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solhihi_d = new double*[degp1];
      sollohi_d = new double*[degp1];
      solhilo_d = new double*[degp1];
      sollolo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhihi_d[i] = new double[dim];
         sollohi_d[i] = new double[dim];
         solhilo_d[i] = new double[dim];
         sollolo_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihi_h;
   double **rhslohi_h;
   double **rhshilo_h;
   double **rhslolo_h;
   double **rhshihi_d;
   double **rhslohi_d;
   double **rhshilo_d;
   double **rhslolo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshihi_h = new double*[degp1];
      rhslohi_h = new double*[degp1];
      rhshilo_h = new double*[degp1];
      rhslolo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshihi_h[i] = new double[dim];
         rhslohi_h[i] = new double[dim];
         rhshilo_h[i] = new double[dim];
         rhslolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhshihi_d = new double*[degp1];
      rhslohi_d = new double*[degp1];
      rhshilo_d = new double*[degp1];
      rhslolo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshihi_d[i] = new double[dim];
         rhslohi_d[i] = new double[dim];
         rhshilo_d[i] = new double[dim];
         rhslolo_d[i] = new double[dim];
      }
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihi_h;
   double **urhslohi_h;
   double **urhshilo_h;
   double **urhslolo_h;
   double **urhshihi_d;
   double **urhslohi_d;
   double **urhshilo_d;
   double **urhslolo_d;

   if((mode == 1) || (mode == 2))
   {
      urhshihi_h = new double*[degp1];
      urhslohi_h = new double*[degp1];
      urhshilo_h = new double*[degp1];
      urhslolo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshihi_h[i] = new double[dim];
         urhslohi_h[i] = new double[dim];
         urhshilo_h[i] = new double[dim];
         urhslolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhshihi_d = new double*[degp1];
      urhslohi_d = new double*[degp1];
      urhshilo_d = new double*[degp1];
      urhslolo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshihi_d[i] = new double[dim];
         urhslohi_d[i] = new double[dim];
         urhshilo_d[i] = new double[dim];
         urhslolo_d[i] = new double[dim];
      }
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihi = new double*[dim];
   double **workmatlohi = new double*[dim];
   double **workmathilo = new double*[dim];
   double **workmatlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihi[i] = new double[dim];
      workmatlohi[i] = new double[dim];
      workmathilo[i] = new double[dim];
      workmatlolo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechihi = new double[dim];
   double *workveclohi = new double[dim];
   double *workvechilo = new double[dim];
   double *workveclolo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihi = new double*[degp1];
   double **workrhslohi = new double*[degp1];
   double **workrhshilo = new double*[degp1];
   double **workrhslolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihi[i] = new double[dim];
      workrhslohi[i] = new double[dim];
      workrhshilo[i] = new double[dim];
      workrhslolo[i] = new double[dim];
   }
   double **resvechihi = new double*[degp1];
   double **resveclohi = new double*[degp1];
   double **resvechilo = new double*[degp1];
   double **resveclolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihi[i] = new double[dim];
      resveclohi[i] = new double[dim];
      resvechilo[i] = new double[dim];
      resveclolo[i] = new double[dim];
   }
   double resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo;

   double **Qhihi_h;
   double **Qlohi_h;
   double **Qhilo_h;
   double **Qlolo_h;
   double **Qhihi_d;
   double **Qlohi_d;
   double **Qhilo_d;
   double **Qlolo_d;

   if((mode == 1) || (mode == 2))
   {
      Qhihi_h = new double*[dim];
      Qlohi_h = new double*[dim];
      Qhilo_h = new double*[dim];
      Qlolo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhihi_h[i] = new double[dim];
         Qlohi_h[i] = new double[dim];
         Qhilo_h[i] = new double[dim];
         Qlolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qhihi_d = new double*[dim];
      Qlohi_d = new double*[dim];
      Qhilo_d = new double*[dim];
      Qlolo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhihi_d[i] = new double[dim];
         Qlohi_d[i] = new double[dim];
         Qhilo_d[i] = new double[dim];
         Qlolo_d[i] = new double[dim];
      }
   }
   double **Rhihi_h;
   double **Rlohi_h;
   double **Rhilo_h;
   double **Rlolo_h;
   double **Rhihi_d;
   double **Rlohi_d;
   double **Rhilo_d;
   double **Rlolo_d;

   if((mode == 1) || (mode == 2))
   {
      Rhihi_h = new double*[dim];
      Rlohi_h = new double*[dim];
      Rhilo_h = new double*[dim];
      Rlolo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhihi_h[i] = new double[dim];
         Rlohi_h[i] = new double[dim];
         Rhilo_h[i] = new double[dim];
         Rlolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rhihi_d = new double*[dim];
      Rlohi_d = new double*[dim];
      Rhilo_d = new double*[dim];
      Rlolo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhihi_d[i] = new double[dim];
         Rlohi_d[i] = new double[dim];
         Rhilo_d[i] = new double[dim];
         Rlolo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **startsolhihi = new double*[dim];
   double **startsollohi = new double*[dim];
   double **startsolhilo = new double*[dim];
   double **startsollolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      startsolhihi[i] = new double[degp1];
      startsollohi[i] = new double[degp1];
      startsolhilo[i] = new double[degp1];
      startsollolo[i] = new double[degp1];
   }
   make_real4_exponentials
      (dim,deg,startsolhihi,startsollohi,startsolhilo,startsollolo);

   if(nbrcol != 1) // generate coefficients for the columns
      make_real4_coefficients(nbrcol,dim,cffhihi,cfflohi,cffhilo,cfflolo);

   // compute the right hand sides via evaluation

   double **mbrhshihi = new double*[dim];
   double **mbrhslohi = new double*[dim];
   double **mbrhshilo = new double*[dim];
   double **mbrhslolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhshihi[i] = new double[degp1];
      mbrhslohi[i] = new double[degp1];
      mbrhshilo[i] = new double[degp1];
      mbrhslolo[i] = new double[degp1];

      mbrhshihi[i][0] = 1.0;     // initialize product to one
      mbrhslohi[i][0] = 0.0;
      mbrhshilo[i][0] = 0.0;
      mbrhslolo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhshihi[i][k] = 0.0;
         mbrhslohi[i][k] = 0.0;
         mbrhshilo[i][k] = 0.0;
         mbrhslolo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_real4_monomials
         (dim,deg,rowsA,startsolhihi,startsollohi,startsolhilo,startsollolo,
          mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo);
   else
      evaluate_real4_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffhihi,cfflohi,cffhilo,cfflolo,
          startsolhihi,startsollohi,startsolhilo,startsollolo,
          mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++) // mbrhs[i][j] = mbrhs[i][j] - mbrhs[i][j-1];
      {
         double acchihi,acclohi,acchilo,acclolo;

         qdf_sub(mbrhshihi[i][j],  mbrhslohi[i][j],
                 mbrhshilo[i][j],  mbrhslolo[i][j],
                 mbrhshihi[i][j-1],mbrhslohi[i][j-1],
                 mbrhshilo[i][j-1],mbrhslolo[i][j-1],
                 &acchihi,&acclohi,&acchilo,&acclolo);
         mbrhshihi[i][j] = acchihi;
         mbrhslohi[i][j] = acclohi;
         mbrhshilo[i][j] = acchilo;
         mbrhslolo[i][j] = acclolo;
      }

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshihi[i][j] << "  " << mbrhslohi[i][j] << endl
                 << "  "
                 << mbrhshilo[i][j] << "  " << mbrhslolo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihi_h[i][j] = startsolhihi[i][j];
         inputlohi_h[i][j] = startsollohi[i][j];
         inputhilo_h[i][j] = startsolhilo[i][j];
         inputlolo_h[i][j] = startsollolo[i][j];
         inputhihi_d[i][j] = startsolhihi[i][j];
         inputlohi_d[i][j] = startsollohi[i][j];
         inputhilo_d[i][j] = startsolhilo[i][j];
         inputlolo_d[i][j] = startsollolo[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihi_h[i][0] << "  "
                            << inputlohi_h[i][0] << endl;
         cout << "     " << inputhilo_h[i][0] << "  "
                         << inputlolo_h[i][0] << endl;
      }
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   double *ratios_d = new double[dim];
   double *ratios_h = new double[dim];
   double step_d,step_h;

   dbl4_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,nvr,idx,exp,nbrfac,expfac,
       mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo,
       cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
       outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
       outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
       funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
       funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
       jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
       jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
       rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
       rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
       urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
       urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
       solhihi_h,sollohi_h,solhilo_h,sollolo_h,
       solhihi_d,sollohi_d,solhilo_d,sollolo_d,
       Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
       Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
       workmathihi,workmatlohi,workmathilo,workmatlolo,
       workvechihi,workveclohi,workvechilo,workveclolo,
       resvechihi,resveclohi,resvechilo,resveclolo,
       &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,vrblvl,mode);

   if((mode == 0) || (mode == 2))
      dbl_fabry_step(dim,deg,inputhihi_d,ratios_d,&step_d,1); // vrblvl);

   if((mode == 1) || (mode == 2))
      dbl_fabry_step(dim,deg,inputhihi_h,ratios_h,&step_h,1); // vrblvl);

   return 0;
}
