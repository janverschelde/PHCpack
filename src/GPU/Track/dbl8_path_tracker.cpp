// The file dbl8_path_tracker.cpp defines the functions with prototypes in
// the file dbl8_path_tracker.h.

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
#include "dbl8_factorizations.h"
#include "dbl8_monomial_systems.h"
#include "octo_double_functions.h"
#include "unimodular_matrices.h"
#include "dbl8_systems_host.h"
#include "dbl8_systems_kernels.h"
#include "dbl8_bals_host.h"
#include "dbl8_tail_kernels.h"
#include "dbl8_bals_kernels.h"
#include "dbl8_newton_testers.h"
#include "dbl8_newton_method.h"
#include "dbl_fabry_host.h"
#include "dbl8_path_tracker.h"

using namespace std;

int dbl8_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **acchihihi, double **acclohihi,
   double **acchilohi, double **acclolohi,
   double **acchihilo, double **acclohilo,
   double **acchilolo, double **acclololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
   double ***outputhihihi_h, double ***outputlohihi_h,
   double ***outputhilohi_h, double ***outputlolohi_h,
   double ***outputhihilo_h, double ***outputlohilo_h,
   double ***outputhilolo_h, double ***outputlololo_h,
   double ***outputhihihi_d, double ***outputlohihi_d,
   double ***outputhilohi_d, double ***outputlolohi_d,
   double ***outputhihilo_d, double ***outputlohilo_d,
   double ***outputhilolo_d, double ***outputlololo_d,
   double **funvalhihihi_h, double **funvallohihi_h,
   double **funvalhilohi_h, double **funvallolohi_h,
   double **funvalhihilo_h, double **funvallohilo_h,
   double **funvalhilolo_h, double **funvallololo_h,
   double **funvalhihihi_d, double **funvallohihi_d,
   double **funvalhilohi_d, double **funvallolohi_d,
   double **funvalhihilo_d, double **funvallohilo_d,
   double **funvalhilolo_d, double **funvallololo_d,
   double ***jacvalhihihi_h, double ***jacvallohihi_h,
   double ***jacvalhilohi_h, double ***jacvallolohi_h,
   double ***jacvalhihilo_h, double ***jacvallohilo_h,
   double ***jacvalhilolo_h, double ***jacvallololo_h,
   double ***jacvalhihihi_d, double ***jacvallohihi_d,
   double ***jacvalhilohi_d, double ***jacvallolohi_d,
   double ***jacvalhihilo_d, double ***jacvallohilo_d,
   double ***jacvalhilolo_d, double ***jacvallololo_d,
   double **rhshihihi_h, double **rhslohihi_h,
   double **rhshilohi_h, double **rhslolohi_h,
   double **rhshihilo_h, double **rhslohilo_h,
   double **rhshilolo_h, double **rhslololo_h,
   double **rhshihihi_d, double **rhslohihi_d,
   double **rhshilohi_d, double **rhslolohi_d,
   double **rhshihilo_d, double **rhslohilo_d,
   double **rhshilolo_d, double **rhslololo_d,
   double **urhshihihi_h, double **urhslohihi_h,
   double **urhshilohi_h, double **urhslolohi_h,
   double **urhshihilo_h, double **urhslohilo_h,
   double **urhshilolo_h, double **urhslololo_h,
   double **urhshihihi_d, double **urhslohihi_d,
   double **urhshilohi_d, double **urhslolohi_d,
   double **urhshihilo_d, double **urhslohilo_d,
   double **urhshilolo_d, double **urhslololo_d,
   double **solhihihi_h, double **sollohihi_h,
   double **solhilohi_h, double **sollolohi_h,
   double **solhihilo_h, double **sollohilo_h,
   double **solhilolo_h, double **sollololo_h,
   double **solhihihi_d, double **sollohihi_d,
   double **solhilohi_d, double **sollolohi_d,
   double **solhihilo_d, double **sollohilo_d,
   double **solhilolo_d, double **sollololo_d,
   double **Qhihihi_h, double **Qlohihi_h,
   double **Qhilohi_h, double **Qlolohi_h,
   double **Qhihilo_h, double **Qlohilo_h,
   double **Qhilolo_h, double **Qlololo_h,
   double **Qhihihi_d, double **Qlohihi_d,
   double **Qhilohi_d, double **Qlolohi_d,
   double **Qhihilo_d, double **Qlohilo_d,
   double **Qhilolo_d, double **Qlololo_d,
   double **Rhihihi_h, double **Rlohihi_h,
   double **Rhilohi_h, double **Rlolohi_h,
   double **Rhihilo_h, double **Rlohilo_h,
   double **Rhilolo_h, double **Rlololo_h,
   double **Rhihihi_d, double **Rlohihi_d,
   double **Rhilohi_d, double **Rlolohi_d,
   double **Rhihilo_d, double **Rlohilo_d,
   double **Rhilolo_d, double **Rlololo_d,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **resvechihihi, double **resveclohihi, 
   double **resvechilohi, double **resveclolohi, 
   double **resvechihilo, double **resveclohilo, 
   double **resvechilolo, double **resveclololo, 
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo,
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

      dbl8_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,
          &tailidx_h,&tailidx_d,nvr,idx,exp,nbrfac,expfac,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,dpr,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
          funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
          funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
          funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
          funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
          rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
          rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
          rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
          rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
          urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
          urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
          urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
          urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
          resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,
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
                           << inputhihihi_d[i][j] << "  "
                           << inputlohihi_d[i][j] << endl
                           << inputhilohi_d[i][j] << "  "
                           << inputlolohi_d[i][j] << endl
                           << inputhihilo_d[i][j] << "  "
                           << inputlohilo_d[i][j] << endl
                           << inputhilolo_d[i][j] << "  "
                           << inputlololo_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhihihi_h[i][j] << "  "
                           << inputlohihi_h[i][j] << endl
                           << inputhilohi_h[i][j] << "  "
                           << inputlolohi_h[i][j] << endl
                           << inputhihilo_h[i][j] << "  "
                           << inputlohilo_h[i][j] << endl
                           << inputhilolo_h[i][j] << "  "
                           << inputlololo_h[i][j] << endl;
         }
         if(mode == 2)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhihihi_d[i][j] << "  "
                           << inputlohihi_d[i][j] << endl
                           << inputhilohi_d[i][j] << "  "
                           << inputlolohi_d[i][j] << endl
                           << inputhihilo_d[i][j] << "  "
                           << inputlohilo_d[i][j] << endl
                           << inputhilolo_d[i][j] << "  "
                           << inputlololo_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhihihi_h[i][j] << "  "
                           << inputlohihi_h[i][j] << endl
                           << inputhilohi_h[i][j] << "  "
                           << inputlolohi_h[i][j] << endl
                           << inputhihilo_h[i][j] << "  "
                           << inputlohilo_h[i][j] << endl
                           << inputhilolo_h[i][j] << "  "
                           << inputlololo_h[i][j] << endl;
            errsum += abs(inputhihihi_h[i][j] - inputhihihi_d[i][j])
                    + abs(inputlohihi_h[i][j] - inputlohihi_d[i][j])
                    + abs(inputhilohi_h[i][j] - inputhilohi_d[i][j])
                    + abs(inputlolohi_h[i][j] - inputlolohi_d[i][j])
                    + abs(inputhihilo_h[i][j] - inputhihilo_d[i][j])
                    + abs(inputlohilo_h[i][j] - inputlohilo_d[i][j])
                    + abs(inputhilolo_h[i][j] - inputhilolo_d[i][j])
                    + abs(inputlololo_h[i][j] - inputlololo_d[i][j]);
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

int test_dbl8_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhihihi_h = new double*[dim];
   double **inputlohihi_h = new double*[dim];
   double **inputhilohi_h = new double*[dim];
   double **inputlolohi_h = new double*[dim];
   double **inputhihilo_h = new double*[dim];
   double **inputlohilo_h = new double*[dim];
   double **inputhilolo_h = new double*[dim];
   double **inputlololo_h = new double*[dim];
   double **inputhihihi_d = new double*[dim];
   double **inputlohihi_d = new double*[dim];
   double **inputhilohi_d = new double*[dim];
   double **inputlolohi_d = new double*[dim];
   double **inputhihilo_d = new double*[dim];
   double **inputlohilo_d = new double*[dim];
   double **inputhilolo_d = new double*[dim];
   double **inputlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihihi_h[i] = new double[degp1];
      inputlohihi_h[i] = new double[degp1];
      inputhilohi_h[i] = new double[degp1];
      inputlolohi_h[i] = new double[degp1];
      inputhihilo_h[i] = new double[degp1];
      inputlohilo_h[i] = new double[degp1];
      inputhilolo_h[i] = new double[degp1];
      inputlololo_h[i] = new double[degp1];
      inputhihihi_d[i] = new double[degp1];
      inputlohihi_d[i] = new double[degp1];
      inputhilohi_d[i] = new double[degp1];
      inputlolohi_d[i] = new double[degp1];
      inputhihilo_d[i] = new double[degp1];
      inputlohilo_d[i] = new double[degp1];
      inputhilolo_d[i] = new double[degp1];
      inputlololo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **acchihihi = new double*[dim+1]; // accumulated power series
   double **acclohihi = new double*[dim+1];
   double **acchilohi = new double*[dim+1];
   double **acclolohi = new double*[dim+1];
   double **acchihilo = new double*[dim+1];
   double **acclohilo = new double*[dim+1];
   double **acchilolo = new double*[dim+1];
   double **acclololo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchihihi[i] = new double[degp1];
      acclohihi[i] = new double[degp1];
      acchilohi[i] = new double[degp1];
      acclolohi[i] = new double[degp1];
      acchihilo[i] = new double[degp1];
      acclohilo[i] = new double[degp1];
      acchilolo[i] = new double[degp1];
      acclololo[i] = new double[degp1];
   }

   double ***cffhihihi = new double**[nbrcol]; // coefficients of monomials
   double ***cfflohihi = new double**[nbrcol];
   double ***cffhilohi = new double**[nbrcol];
   double ***cfflolohi = new double**[nbrcol];
   double ***cffhihilo = new double**[nbrcol];
   double ***cfflohilo = new double**[nbrcol];
   double ***cffhilolo = new double**[nbrcol];
   double ***cfflololo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffhihihi[i] = new double*[dim];
      cfflohihi[i] = new double*[dim];
      cffhilohi[i] = new double*[dim];
      cfflolohi[i] = new double*[dim];
      cffhihilo[i] = new double*[dim];
      cfflohilo[i] = new double*[dim];
      cffhilolo[i] = new double*[dim];
      cfflololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffhihihi[i][j] = new double[degp1];
         cfflohihi[i][j] = new double[degp1];
         cffhilohi[i][j] = new double[degp1];
         cfflolohi[i][j] = new double[degp1];
         cffhihilo[i][j] = new double[degp1];
         cfflohilo[i][j] = new double[degp1];
         cffhilolo[i][j] = new double[degp1];
         cfflololo[i][j] = new double[degp1];

         for(int k=0; k<deg; k++)
         {
            cffhihihi[i][j][k] = 0.0;
            cfflohihi[i][j][k] = 0.0;
            cffhilohi[i][j][k] = 0.0;
            cfflolohi[i][j][k] = 0.0;
            cffhihilo[i][j][k] = 0.0;
            cfflohilo[i][j][k] = 0.0;
            cffhilolo[i][j][k] = 0.0;
            cfflololo[i][j][k] = 0.0;
         }
      }
   }
   double ***outputhihihi_h;
   double ***outputlohihi_h;
   double ***outputhilohi_h;
   double ***outputlolohi_h;
   double ***outputhihilo_h;
   double ***outputlohilo_h;
   double ***outputhilolo_h;
   double ***outputlololo_h;
   double ***outputhihihi_d;
   double ***outputlohihi_d;
   double ***outputhilohi_d;
   double ***outputlolohi_d;
   double ***outputhihilo_d;
   double ***outputlohilo_d;
   double ***outputhilolo_d;
   double ***outputlololo_d;

   if((mode == 1) || (mode == 2))
   {
      outputhihihi_h = new double**[dim];
      outputlohihi_h = new double**[dim];
      outputhilohi_h = new double**[dim];
      outputlolohi_h = new double**[dim];
      outputhihilo_h = new double**[dim];
      outputlohilo_h = new double**[dim];
      outputhilolo_h = new double**[dim];
      outputlololo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhihihi_h[i] = new double*[dim+1];
         outputlohihi_h[i] = new double*[dim+1];
         outputhilohi_h[i] = new double*[dim+1];
         outputlolohi_h[i] = new double*[dim+1];
         outputhihilo_h[i] = new double*[dim+1];
         outputlohilo_h[i] = new double*[dim+1];
         outputhilolo_h[i] = new double*[dim+1];
         outputlololo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhihihi_h[i][j] = new double[degp1];
            outputlohihi_h[i][j] = new double[degp1];
            outputhilohi_h[i][j] = new double[degp1];
            outputlolohi_h[i][j] = new double[degp1];
            outputhihilo_h[i][j] = new double[degp1];
            outputlohilo_h[i][j] = new double[degp1];
            outputhilolo_h[i][j] = new double[degp1];
            outputlololo_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputhihihi_d = new double**[dim];
      outputlohihi_d = new double**[dim];
      outputhilohi_d = new double**[dim];
      outputlolohi_d = new double**[dim];
      outputhihilo_d = new double**[dim];
      outputlohilo_d = new double**[dim];
      outputhilolo_d = new double**[dim];
      outputlololo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhihihi_d[i] = new double*[dim+1];
         outputlohihi_d[i] = new double*[dim+1];
         outputhilohi_d[i] = new double*[dim+1];
         outputlolohi_d[i] = new double*[dim+1];
         outputhihilo_d[i] = new double*[dim+1];
         outputlohilo_d[i] = new double*[dim+1];
         outputhilolo_d[i] = new double*[dim+1];
         outputlololo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhihihi_d[i][j] = new double[degp1];
            outputlohihi_d[i][j] = new double[degp1];
            outputhilohi_d[i][j] = new double[degp1];
            outputlolohi_d[i][j] = new double[degp1];
            outputhihilo_d[i][j] = new double[degp1];
            outputlohilo_d[i][j] = new double[degp1];
            outputhilolo_d[i][j] = new double[degp1];
            outputlololo_d[i][j] = new double[degp1];
         }
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihihi_h;
   double **funvallohihi_h;
   double **funvalhilohi_h;
   double **funvallolohi_h;
   double **funvalhihilo_h;
   double **funvallohilo_h;
   double **funvalhilolo_h;
   double **funvallololo_h;
   double **funvalhihihi_d;
   double **funvallohihi_d;
   double **funvalhilohi_d;
   double **funvallolohi_d;
   double **funvalhihilo_d;
   double **funvallohilo_d;
   double **funvalhilolo_d;
   double **funvallololo_d;

   if((mode == 1) || (mode == 2))
   {
      funvalhihihi_h = new double*[dim];
      funvallohihi_h = new double*[dim];
      funvalhilohi_h = new double*[dim];
      funvallolohi_h = new double*[dim];
      funvalhihilo_h = new double*[dim];
      funvallohilo_h = new double*[dim];
      funvalhilolo_h = new double*[dim];
      funvallololo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhihihi_h[i] = new double[degp1];
         funvallohihi_h[i] = new double[degp1];
         funvalhilohi_h[i] = new double[degp1];
         funvallolohi_h[i] = new double[degp1];
         funvalhihilo_h[i] = new double[degp1];
         funvallohilo_h[i] = new double[degp1];
         funvalhilolo_h[i] = new double[degp1];
         funvallololo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalhihihi_d = new double*[dim];
      funvallohihi_d = new double*[dim];
      funvalhilohi_d = new double*[dim];
      funvallolohi_d = new double*[dim];
      funvalhihilo_d = new double*[dim];
      funvallohilo_d = new double*[dim];
      funvalhilolo_d = new double*[dim];
      funvallololo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhihihi_d[i] = new double[degp1];
         funvallohihi_d[i] = new double[degp1];
         funvalhilohi_d[i] = new double[degp1];
         funvallolohi_d[i] = new double[degp1];
         funvalhihilo_d[i] = new double[degp1];
         funvallohilo_d[i] = new double[degp1];
         funvalhilolo_d[i] = new double[degp1];
         funvallololo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihihi_h;
   double ***jacvallohihi_h;
   double ***jacvalhilohi_h;
   double ***jacvallolohi_h;
   double ***jacvalhihilo_h;
   double ***jacvallohilo_h;
   double ***jacvalhilolo_h;
   double ***jacvallololo_h;
   double ***jacvalhihihi_d;
   double ***jacvallohihi_d;
   double ***jacvalhilohi_d;
   double ***jacvallolohi_d;
   double ***jacvalhihilo_d;
   double ***jacvallohilo_d;
   double ***jacvalhilolo_d;
   double ***jacvallololo_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalhihihi_h = new double**[degp1];
      jacvallohihi_h = new double**[degp1];
      jacvalhilohi_h = new double**[degp1];
      jacvallolohi_h = new double**[degp1];
      jacvalhihilo_h = new double**[degp1];
      jacvallohilo_h = new double**[degp1];
      jacvalhilolo_h = new double**[degp1];
      jacvallololo_h = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhihihi_h[i] = new double*[dim];
         jacvallohihi_h[i] = new double*[dim];
         jacvalhilohi_h[i] = new double*[dim];
         jacvallolohi_h[i] = new double*[dim];
         jacvalhihilo_h[i] = new double*[dim];
         jacvallohilo_h[i] = new double*[dim];
         jacvalhilolo_h[i] = new double*[dim];
         jacvallololo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhihihi_h[i][j] = new double[dim];
            jacvallohihi_h[i][j] = new double[dim];
            jacvalhilohi_h[i][j] = new double[dim];
            jacvallolohi_h[i][j] = new double[dim];
            jacvalhihilo_h[i][j] = new double[dim];
            jacvallohilo_h[i][j] = new double[dim];
            jacvalhilolo_h[i][j] = new double[dim];
            jacvallololo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      jacvalhihihi_d = new double**[degp1];
      jacvallohihi_d = new double**[degp1];
      jacvalhilohi_d = new double**[degp1];
      jacvallolohi_d = new double**[degp1];
      jacvalhihilo_d = new double**[degp1];
      jacvallohilo_d = new double**[degp1];
      jacvalhilolo_d = new double**[degp1];
      jacvallololo_d = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhihihi_d[i] = new double*[dim];
         jacvallohihi_d[i] = new double*[dim];
         jacvalhilohi_d[i] = new double*[dim];
         jacvallolohi_d[i] = new double*[dim];
         jacvalhihilo_d[i] = new double*[dim];
         jacvallohilo_d[i] = new double*[dim];
         jacvalhilolo_d[i] = new double*[dim];
         jacvallololo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhihihi_d[i][j] = new double[dim];
            jacvallohihi_d[i][j] = new double[dim];
            jacvalhilohi_d[i][j] = new double[dim];
            jacvallolohi_d[i][j] = new double[dim];
            jacvalhihilo_d[i][j] = new double[dim];
            jacvallohilo_d[i][j] = new double[dim];
            jacvalhilolo_d[i][j] = new double[dim];
            jacvallololo_d[i][j] = new double[dim];
         }
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihihi_h;
   double **sollohihi_h;
   double **solhilohi_h;
   double **sollolohi_h;
   double **solhihilo_h;
   double **sollohilo_h;
   double **solhilolo_h;
   double **sollololo_h;
   double **solhihihi_d;
   double **sollohihi_d;
   double **solhilohi_d;
   double **sollolohi_d;
   double **solhihilo_d;
   double **sollohilo_d;
   double **solhilolo_d;
   double **sollololo_d;

   if((mode == 1) || (mode == 2))
   {
      solhihihi_h = new double*[degp1];
      sollohihi_h = new double*[degp1];
      solhilohi_h = new double*[degp1];
      sollolohi_h = new double*[degp1];
      solhihilo_h = new double*[degp1];
      sollohilo_h = new double*[degp1];
      solhilolo_h = new double*[degp1];
      sollololo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhihihi_h[i] = new double[dim];
         sollohihi_h[i] = new double[dim];
         solhilohi_h[i] = new double[dim];
         sollolohi_h[i] = new double[dim];
         solhihilo_h[i] = new double[dim];
         sollohilo_h[i] = new double[dim];
         solhilolo_h[i] = new double[dim];
         sollololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solhihihi_d = new double*[degp1];
      sollohihi_d = new double*[degp1];
      solhilohi_d = new double*[degp1];
      sollolohi_d = new double*[degp1];
      solhihilo_d = new double*[degp1];
      sollohilo_d = new double*[degp1];
      solhilolo_d = new double*[degp1];
      sollololo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhihihi_d[i] = new double[dim];
         sollohihi_d[i] = new double[dim];
         solhilohi_d[i] = new double[dim];
         sollolohi_d[i] = new double[dim];
         solhihilo_d[i] = new double[dim];
         sollohilo_d[i] = new double[dim];
         solhilolo_d[i] = new double[dim];
         sollololo_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihihi_h;
   double **rhslohihi_h;
   double **rhshilohi_h;
   double **rhslolohi_h;
   double **rhshihilo_h;
   double **rhslohilo_h;
   double **rhshilolo_h;
   double **rhslololo_h;
   double **rhshihihi_d;
   double **rhslohihi_d;
   double **rhshilohi_d;
   double **rhslolohi_d;
   double **rhshihilo_d;
   double **rhslohilo_d;
   double **rhshilolo_d;
   double **rhslololo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshihihi_h = new double*[degp1];
      rhslohihi_h = new double*[degp1];
      rhshilohi_h = new double*[degp1];
      rhslolohi_h = new double*[degp1];
      rhshihilo_h = new double*[degp1];
      rhslohilo_h = new double*[degp1];
      rhshilolo_h = new double*[degp1];
      rhslololo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshihihi_h[i] = new double[dim];
         rhslohihi_h[i] = new double[dim];
         rhshilohi_h[i] = new double[dim];
         rhslolohi_h[i] = new double[dim];
         rhshihilo_h[i] = new double[dim];
         rhslohilo_h[i] = new double[dim];
         rhshilolo_h[i] = new double[dim];
         rhslololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhshihihi_d = new double*[degp1];
      rhslohihi_d = new double*[degp1];
      rhshilohi_d = new double*[degp1];
      rhslolohi_d = new double*[degp1];
      rhshihilo_d = new double*[degp1];
      rhslohilo_d = new double*[degp1];
      rhshilolo_d = new double*[degp1];
      rhslololo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshihihi_d[i] = new double[dim];
         rhslohihi_d[i] = new double[dim];
         rhshilohi_d[i] = new double[dim];
         rhslolohi_d[i] = new double[dim];
         rhshihilo_d[i] = new double[dim];
         rhslohilo_d[i] = new double[dim];
         rhshilolo_d[i] = new double[dim];
         rhslololo_d[i] = new double[dim];
      }
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihihi_h;
   double **urhslohihi_h;
   double **urhshilohi_h;
   double **urhslolohi_h;
   double **urhshihilo_h;
   double **urhslohilo_h;
   double **urhshilolo_h;
   double **urhslololo_h;
   double **urhshihihi_d;
   double **urhslohihi_d;
   double **urhshilohi_d;
   double **urhslolohi_d;
   double **urhshihilo_d;
   double **urhslohilo_d;
   double **urhshilolo_d;
   double **urhslololo_d;

   if((mode == 1) || (mode == 2))
   {
      urhshihihi_h = new double*[degp1];
      urhslohihi_h = new double*[degp1];
      urhshilohi_h = new double*[degp1];
      urhslolohi_h = new double*[degp1];
      urhshihilo_h = new double*[degp1];
      urhslohilo_h = new double*[degp1];
      urhshilolo_h = new double*[degp1];
      urhslololo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshihihi_h[i] = new double[dim];
         urhslohihi_h[i] = new double[dim];
         urhshilohi_h[i] = new double[dim];
         urhslolohi_h[i] = new double[dim];
         urhshihilo_h[i] = new double[dim];
         urhslohilo_h[i] = new double[dim];
         urhshilolo_h[i] = new double[dim];
         urhslololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhshihihi_d = new double*[degp1];
      urhslohihi_d = new double*[degp1];
      urhshilohi_d = new double*[degp1];
      urhslolohi_d = new double*[degp1];
      urhshihilo_d = new double*[degp1];
      urhslohilo_d = new double*[degp1];
      urhshilolo_d = new double*[degp1];
      urhslololo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshihihi_d[i] = new double[dim];
         urhslohihi_d[i] = new double[dim];
         urhshilohi_d[i] = new double[dim];
         urhslolohi_d[i] = new double[dim];
         urhshihilo_d[i] = new double[dim];
         urhslohilo_d[i] = new double[dim];
         urhshilolo_d[i] = new double[dim];
         urhslololo_d[i] = new double[dim];
      }
   }
   double *workvechihihi = new double[dim];
   double *workveclohihi = new double[dim];
   double *workvechilohi = new double[dim];
   double *workveclolohi = new double[dim];
   double *workvechihilo = new double[dim];
   double *workveclohilo = new double[dim];
   double *workvechilolo = new double[dim];
   double *workveclololo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihihi = new double*[degp1];
   double **workrhslohihi = new double*[degp1];
   double **workrhshilohi = new double*[degp1];
   double **workrhslolohi = new double*[degp1];
   double **workrhshihilo = new double*[degp1];
   double **workrhslohilo = new double*[degp1];
   double **workrhshilolo = new double*[degp1];
   double **workrhslololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihihi[i] = new double[dim];
      workrhslohihi[i] = new double[dim];
      workrhshilohi[i] = new double[dim];
      workrhslolohi[i] = new double[dim];
      workrhshihilo[i] = new double[dim];
      workrhslohilo[i] = new double[dim];
      workrhshilolo[i] = new double[dim];
      workrhslololo[i] = new double[dim];
   }
   double **resvechihihi = new double*[degp1];
   double **resveclohihi = new double*[degp1];
   double **resvechilohi = new double*[degp1];
   double **resveclolohi = new double*[degp1];
   double **resvechihilo = new double*[degp1];
   double **resveclohilo = new double*[degp1];
   double **resvechilolo = new double*[degp1];
   double **resveclololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihihi[i] = new double[dim];
      resveclohihi[i] = new double[dim];
      resvechilohi[i] = new double[dim];
      resveclolohi[i] = new double[dim];
      resvechihilo[i] = new double[dim];
      resveclohilo[i] = new double[dim];
      resvechilolo[i] = new double[dim];
      resveclololo[i] = new double[dim];
   }
   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

   double **Qhihihi_h;
   double **Qlohihi_h;
   double **Qhilohi_h;
   double **Qlolohi_h;
   double **Qhihilo_h;
   double **Qlohilo_h;
   double **Qhilolo_h;
   double **Qlololo_h;
   double **Qhihihi_d;
   double **Qlohihi_d;
   double **Qhilohi_d;
   double **Qlolohi_d;
   double **Qhihilo_d;
   double **Qlohilo_d;
   double **Qhilolo_d;
   double **Qlololo_d;

   if((mode == 1) || (mode == 2))
   {
      Qhihihi_h = new double*[dim];
      Qlohihi_h = new double*[dim];
      Qhilohi_h = new double*[dim];
      Qlolohi_h = new double*[dim];
      Qhihilo_h = new double*[dim];
      Qlohilo_h = new double*[dim];
      Qhilolo_h = new double*[dim];
      Qlololo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhihihi_h[i] = new double[dim];
         Qlohihi_h[i] = new double[dim];
         Qhilohi_h[i] = new double[dim];
         Qlolohi_h[i] = new double[dim];
         Qhihilo_h[i] = new double[dim];
         Qlohilo_h[i] = new double[dim];
         Qhilolo_h[i] = new double[dim];
         Qlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qhihihi_d = new double*[dim];
      Qlohihi_d = new double*[dim];
      Qhilohi_d = new double*[dim];
      Qlolohi_d = new double*[dim];
      Qhihilo_d = new double*[dim];
      Qlohilo_d = new double*[dim];
      Qhilolo_d = new double*[dim];
      Qlololo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhihihi_d[i] = new double[dim];
         Qlohihi_d[i] = new double[dim];
         Qhilohi_d[i] = new double[dim];
         Qlolohi_d[i] = new double[dim];
         Qhihilo_d[i] = new double[dim];
         Qlohilo_d[i] = new double[dim];
         Qhilolo_d[i] = new double[dim];
         Qlololo_d[i] = new double[dim];
      }
   }
   double **Rhihihi_h;
   double **Rlohihi_h;
   double **Rhilohi_h;
   double **Rlolohi_h;
   double **Rhihilo_h;
   double **Rlohilo_h;
   double **Rhilolo_h;
   double **Rlololo_h;
   double **Rhihihi_d;
   double **Rlohihi_d;
   double **Rhilohi_d;
   double **Rlolohi_d;
   double **Rhihilo_d;
   double **Rlohilo_d;
   double **Rhilolo_d;
   double **Rlololo_d;

   if((mode == 1) || (mode == 2))
   {
      Rhihihi_h = new double*[dim];
      Rlohihi_h = new double*[dim];
      Rhilohi_h = new double*[dim];
      Rlolohi_h = new double*[dim];
      Rhihilo_h = new double*[dim];
      Rlohilo_h = new double*[dim];
      Rhilolo_h = new double*[dim];
      Rlololo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhihihi_h[i] = new double[dim];
         Rlohihi_h[i] = new double[dim];
         Rhilohi_h[i] = new double[dim];
         Rlolohi_h[i] = new double[dim];
         Rhihilo_h[i] = new double[dim];
         Rlohilo_h[i] = new double[dim];
         Rhilolo_h[i] = new double[dim];
         Rlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rhihihi_d = new double*[dim];
      Rlohihi_d = new double*[dim];
      Rhilohi_d = new double*[dim];
      Rlolohi_d = new double*[dim];
      Rhihilo_d = new double*[dim];
      Rlohilo_d = new double*[dim];
      Rhilolo_d = new double*[dim];
      Rlololo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhihihi_d[i] = new double[dim];
         Rlohihi_d[i] = new double[dim];
         Rhilohi_d[i] = new double[dim];
         Rlolohi_d[i] = new double[dim];
         Rhihilo_d[i] = new double[dim];
         Rlohilo_d[i] = new double[dim];
         Rhilolo_d[i] = new double[dim];
         Rlololo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **startsolhihihi = new double*[dim];
   double **startsollohihi = new double*[dim];
   double **startsolhilohi = new double*[dim];
   double **startsollolohi = new double*[dim];
   double **startsolhihilo = new double*[dim];
   double **startsollohilo = new double*[dim];
   double **startsolhilolo = new double*[dim];
   double **startsollololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      startsolhihihi[i] = new double[degp1];
      startsollohihi[i] = new double[degp1];
      startsolhilohi[i] = new double[degp1];
      startsollolohi[i] = new double[degp1];
      startsolhihilo[i] = new double[degp1];
      startsollohilo[i] = new double[degp1];
      startsolhilolo[i] = new double[degp1];
      startsollololo[i] = new double[degp1];
   }
   make_real8_exponentials
      (dim,deg,startsolhihihi,startsollohihi,startsolhilohi,startsollolohi,
               startsolhihilo,startsollohilo,startsolhilolo,startsollololo);

   if(nbrcol != 1) // generate coefficients for the columns
      make_real8_coefficients
         (nbrcol,dim,cffhihihi,cfflohihi,cffhilohi,cfflolohi,
                     cffhihilo,cfflohilo,cffhilolo,cfflololo);

   // compute the right hand sides via evaluation

   double **mbrhshihihi = new double*[dim];
   double **mbrhslohihi = new double*[dim];
   double **mbrhshilohi = new double*[dim];
   double **mbrhslolohi = new double*[dim];
   double **mbrhshihilo = new double*[dim];
   double **mbrhslohilo = new double*[dim];
   double **mbrhshilolo = new double*[dim];
   double **mbrhslololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhshihihi[i] = new double[degp1];
      mbrhslohihi[i] = new double[degp1];
      mbrhshilohi[i] = new double[degp1];
      mbrhslolohi[i] = new double[degp1];
      mbrhshihilo[i] = new double[degp1];
      mbrhslohilo[i] = new double[degp1];
      mbrhshilolo[i] = new double[degp1];
      mbrhslololo[i] = new double[degp1];

      mbrhshihihi[i][0] = 1.0;     // initialize product to one
      mbrhslohihi[i][0] = 0.0;
      mbrhshilohi[i][0] = 0.0;
      mbrhslolohi[i][0] = 0.0;
      mbrhshihilo[i][0] = 0.0; 
      mbrhslohilo[i][0] = 0.0;
      mbrhshilolo[i][0] = 0.0;
      mbrhslololo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhshihihi[i][k] = 0.0; mbrhslohihi[i][k] = 0.0;
         mbrhshilohi[i][k] = 0.0; mbrhslolohi[i][k] = 0.0;
         mbrhshihilo[i][k] = 0.0; mbrhslohilo[i][k] = 0.0;
         mbrhshilolo[i][k] = 0.0; mbrhslololo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_real8_monomials
         (dim,deg,rowsA,
          startsolhihihi,startsollohihi,startsolhilohi,startsollolohi,
          startsolhihilo,startsollohilo,startsolhilolo,startsollololo,
          mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
          mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo);
   else
      evaluate_real8_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          startsolhihihi,startsollohihi,startsolhilohi,startsollolohi,
          startsolhihilo,startsollohilo,startsolhilolo,startsollololo,
          mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
          mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++) // mbrhs[i][j] = mbrhs[i][j] - mbrhs[i][j-1];
      {
         double acchihihi,acclohihi,acchilohi,acclolohi;
         double acchihilo,acclohilo,acchilolo,acclololo;

         odf_sub(mbrhshihihi[i][j],  mbrhslohihi[i][j],
                 mbrhshilohi[i][j],  mbrhslolohi[i][j],
                 mbrhshihilo[i][j],  mbrhslohilo[i][j],
                 mbrhshilolo[i][j],  mbrhslololo[i][j],
                 mbrhshihihi[i][j-1],mbrhslohihi[i][j-1],
                 mbrhshilohi[i][j-1],mbrhslolohi[i][j-1],
                 mbrhshihilo[i][j-1],mbrhslohilo[i][j-1],
                 mbrhshilolo[i][j-1],mbrhslololo[i][j-1],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         mbrhshihihi[i][j] = acchihihi;
         mbrhslohihi[i][j] = acclohihi;
         mbrhshilohi[i][j] = acchilohi;
         mbrhslolohi[i][j] = acclolohi;
         mbrhshihilo[i][j] = acchihilo;
         mbrhslohilo[i][j] = acclohilo;
         mbrhshilolo[i][j] = acchilolo;
         mbrhslololo[i][j] = acclololo;
      }

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshihihi[i][j] << "  " << mbrhslohihi[i][j] << endl
                 << "  "
                 << mbrhshilohi[i][j] << "  " << mbrhslolohi[i][j] << endl
                 << "  "
                 << mbrhshihilo[i][j] << "  " << mbrhslohilo[i][j] << endl
                 << "  "
                 << mbrhshilolo[i][j] << "  " << mbrhslololo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihihi_h[i][j] = startsolhihihi[i][j];
         inputlohihi_h[i][j] = startsollohihi[i][j];
         inputhilohi_h[i][j] = startsolhilohi[i][j];
         inputlolohi_h[i][j] = startsollolohi[i][j];
         inputhihilo_h[i][j] = startsolhihilo[i][j];
         inputlohilo_h[i][j] = startsollohilo[i][j];
         inputhilolo_h[i][j] = startsolhilolo[i][j];
         inputlololo_h[i][j] = startsollololo[i][j];
         inputhihihi_d[i][j] = startsolhihihi[i][j];
         inputlohihi_d[i][j] = startsollohihi[i][j];
         inputhilohi_d[i][j] = startsolhilohi[i][j];
         inputlolohi_d[i][j] = startsollolohi[i][j];
         inputhihilo_d[i][j] = startsolhihilo[i][j];
         inputlohilo_d[i][j] = startsollohilo[i][j];
         inputhilolo_d[i][j] = startsolhilolo[i][j];
         inputlololo_d[i][j] = startsollololo[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihihi_h[i][0] << "  "
                            << inputlohihi_h[i][0] << endl;
         cout << "     " << inputhilohi_h[i][0] << "  "
                         << inputlolohi_h[i][0] << endl;
         cout << "     " << inputhihilo_h[i][0] << "  "
                         << inputlohilo_h[i][0] << endl;
         cout << "     " << inputhilolo_h[i][0] << "  "
                         << inputlololo_h[i][0] << endl;
      }
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   double *ratios_d = new double[dim];
   double *ratios_h = new double[dim];
   double step_d,step_h;

   dbl8_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,nvr,idx,exp,nbrfac,expfac,
       mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
       mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       acchihihi,acclohihi,acchilohi,acclolohi,
       acchihilo,acclohilo,acchilolo,acclololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
       outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
       outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
       outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
       outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
       funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
       funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
       funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
       funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
       jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
       jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
       jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
       jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
       rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
       rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
       rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
       rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
       urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
       urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
       urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
       urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
       solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
       solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
       solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
       solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
       Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
       Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
       Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
       Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
       Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
       Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
       Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
       Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
       workvechihihi,workveclohihi,workvechilohi,workveclolohi,
       workvechihilo,workveclohilo,workvechilolo,workveclololo,
       resvechihihi,resveclohihi,resvechilohi,resveclolohi,
       resvechihilo,resveclohilo,resvechilolo,resveclololo,
       &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
       &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
       vrblvl,mode);

   if((mode == 0) || (mode == 2))
      dbl_fabry_step(dim,deg,inputhihihi_d,ratios_d,&step_d,1); // vrblvl);

   if((mode == 1) || (mode == 2))
      dbl_fabry_step(dim,deg,inputhihihi_h,ratios_h,&step_h,1); // vrblvl);

   return 0;
}
