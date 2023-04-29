// The file dbl2_path_tracker.cpp defines the functions with prototypes in
// the file dbl2_path_tracker.h.

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
#include "dbl2_newton_method.h"
#include "dbl_fabry_host.h"
#include "dbl2_path_tracker.h"

using namespace std;

int dbl2_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **rhshi, double **rhslo,
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

      dbl2_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,rhshi,rhslo,dpr,cffhi,cfflo,acchi,acclo,
          inputhi_h,inputlo_h,inputhi_d,inputlo_d,
          outputhi_h,outputlo_h,outputhi_d,outputlo_d,
          funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,urhshi_h,urhslo_h,urhshi_d,urhslo_d,
          solhi_h,sollo_h,solhi_d,sollo_d,
          Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
          workvechi,workveclo,resvechi,resveclo,resmaxhi,resmaxlo,
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
                           << inputhi_d[i][j] << "  "
                           << inputlo_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhi_h[i][j] << "  "
                           << inputlo_h[i][j] << endl;
         }
         if(mode == 2)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhi_d[i][j] << "  "
                           << inputlo_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhi_h[i][j] << "  "
                           << inputlo_h[i][j] << endl;
            errsum += abs(inputhi_h[i][j] - inputhi_d[i][j])
                    + abs(inputlo_h[i][j] - inputlo_d[i][j]);
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

int test_dbl2_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhi_h = new double*[dim];
   double **inputlo_h = new double*[dim];
   double **inputhi_d = new double*[dim];
   double **inputlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputhi_h[i] = new double[degp1];
       inputlo_h[i] = new double[degp1];
       inputhi_d[i] = new double[degp1];
       inputlo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **acchi = new double*[dim+1]; // accumulate series in one column
   double **acclo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchi[i] = new double[degp1];
      acclo[i] = new double[degp1];
   }
   double ***cffhi = new double**[nbrcol];
   double ***cfflo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffhi[i] = new double*[dim]; // the coefficients of monomials
      cfflo[i] = new double*[dim];
      for(int j=0; j<dim; j++)
      {
         cffhi[i][j] = new double[degp1];
         cfflo[i][j] = new double[degp1];
         for(int k=0; k<deg; k++)
         {
            cffhi[i][j][k] = 0.0;
            cfflo[i][j][k] = 0.0;
         }
      }
   }
   double ***outputhi_h;
   double ***outputlo_h;
   double ***outputhi_d;
   double ***outputlo_d;

   if((mode == 1) || (mode == 2))
   {
      outputhi_h = new double**[dim];
      outputlo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhi_h[i] = new double*[dim+1];
         outputlo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhi_h[i][j] = new double[degp1];
            outputlo_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputhi_d = new double**[dim];
      outputlo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhi_d[i] = new double*[dim+1];
         outputlo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhi_d[i][j] = new double[degp1];
            outputlo_d[i][j] = new double[degp1];
         }
      }
   }
   double **funvalhi_h;  // high doubles of function values on host
   double **funvallo_h;  // low doubles of function values on host
   double **funvalhi_d;  // high doubles of function values on device
   double **funvallo_d;  // low doubles of function values on device

   if((mode == 1) || (mode == 2))
   {
      funvalhi_h = new double*[dim];
      funvallo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhi_h[i] = new double[degp1];
         funvallo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalhi_d = new double*[dim];
      funvallo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhi_d[i] = new double[degp1];
         funvallo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h;
   double ***jacvallo_h;
   double ***jacvalhi_d;
   double ***jacvallo_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalhi_h = new double**[degp1];
      jacvallo_h = new double**[degp1];

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
      jacvalhi_d = new double**[degp1];
      jacvallo_d = new double**[degp1];

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
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhi_h;
   double **sollo_h;
   double **solhi_d;
   double **sollo_d;

   if((mode == 1) || (mode == 2))
   {
      solhi_h = new double*[degp1];
      sollo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhi_h[i] = new double[dim];
         sollo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solhi_d = new double*[degp1];
      sollo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhi_d[i] = new double[dim];
         sollo_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshi_h = new double*[degp1];
   double **rhslo_h = new double*[degp1];
   double **rhshi_d = new double*[degp1];
   double **rhslo_d = new double*[degp1];

   if((mode == 1) || (mode == 2))
   {
      rhshi_h = new double*[degp1];
      rhslo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshi_h[i] = new double[dim];
         rhslo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhshi_d = new double*[degp1];
      rhslo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshi_d[i] = new double[dim];
         rhslo_d[i] = new double[dim];
      }
   }
   double *workvechi = new double[dim];
   double *workveclo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **urhshi_h;
   double **urhslo_h;
   double **urhshi_d;
   double **urhslo_d;

   if((mode == 1) || (mode == 2))
   {
      urhshi_h = new double*[degp1];
      urhslo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshi_h[i] = new double[dim];
         urhslo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhshi_d = new double*[degp1];
      urhslo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshi_d[i] = new double[dim];
         urhslo_d[i] = new double[dim];
      }
   }
   double **resvechi = new double*[degp1];
   double **resveclo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechi[i] = new double[dim];
      resveclo[i] = new double[dim];
   }
   double resmaxhi;
   double resmaxlo;

   double **Qhi_h;
   double **Qlo_h;
   double **Qhi_d;
   double **Qlo_d;

   if((mode == 1) || (mode == 2))
   {
      Qhi_h = new double*[dim];
      Qlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhi_h[i] = new double[dim];
         Qlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qhi_d = new double*[dim];
      Qlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhi_d[i] = new double[dim];
         Qlo_d[i] = new double[dim];
      }
   }
   double **Rhi_h;
   double **Rlo_h;
   double **Rhi_d;
   double **Rlo_d;

   if((mode == 1) || (mode == 2))
   {
      Rhi_h = new double*[dim];
      Rlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhi_h[i] = new double[dim];
         Rlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rhi_d = new double*[dim];
      Rlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhi_d[i] = new double[dim];
         Rlo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **startsolhi = new double*[dim];
   double **startsollo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      startsolhi[i] = new double[degp1];
      startsollo[i] = new double[degp1];
   }
   make_real2_exponentials(dim,deg,startsolhi,startsollo);

   if(nbrcol != 1) // generate coefficients for the columns
   {
      // sets only the leading coefficient to a random double ...
      make_real2_coefficients(nbrcol,dim,cffhi,cfflo);
   }
   // compute the right hand sides via evaluation

   double **mbrhshi = new double*[dim];
   double **mbrhslo = new double*[dim];

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
         (dim,deg,rowsA,startsolhi,startsollo,mbrhshi,mbrhslo);
   else
      evaluate_real2_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffhi,cfflo,startsolhi,startsollo,
          mbrhshi,mbrhslo,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++) // mbrhs[i][j] = mbrhs[i][j] - mbrhs[i][j-1];
      {
         double acchi,acclo;

         ddf_sub(mbrhshi[i][j],  mbrhslo[i][j],
                 mbrhshi[i][j-1],mbrhslo[i][j-1],&acchi,&acclo);
         mbrhshi[i][j] = acchi;
         mbrhslo[i][j] = acchi;
      }

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshi[i][j] << "  " << mbrhslo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhi_h[i][j] = startsolhi[i][j];
         inputlo_h[i][j] = startsollo[i][j];
         inputhi_d[i][j] = startsolhi[i][j];
         inputlo_d[i][j] = startsollo[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputhi_h[i][0] << "  "
                            << inputlo_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   double *ratios_d = new double[dim];
   double *ratios_h = new double[dim];
   double step_d,step_h;

   dbl2_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,
       nvr,idx,exp,nbrfac,expfac,mbrhshi,mbrhslo,cffhi,cfflo,acchi,acclo,
       inputhi_h,inputlo_h,inputhi_d,inputlo_d,
       outputhi_h,outputlo_h,outputhi_d,outputlo_d,
       funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
       jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
       rhshi_h,rhslo_h,rhshi_d,rhslo_d,urhshi_h,urhslo_h,urhshi_d,urhslo_d,
       solhi_h,sollo_h,solhi_d,sollo_d,
       Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
       workvechi,workveclo,resvechi,resveclo,&resmaxhi,&resmaxlo,
       vrblvl,mode);

   if((mode == 0) || (mode == 2))
      dbl_fabry_step(dim,deg,inputhi_d,ratios_d,&step_d,1); // vrblvl);

   if((mode == 1) || (mode == 2))
      dbl_fabry_step(dim,deg,inputhi_h,ratios_h,&step_h,1); // vrblvl);

   return 0;
}
