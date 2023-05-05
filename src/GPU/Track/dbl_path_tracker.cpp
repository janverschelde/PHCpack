// The file dbl_path_tracker.cpp defines the functions with prototypes in
// the file dbl_path_tracker.h.

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
#include "dbl_factorizations.h"
#include "dbl_monomial_systems.h"
#include "unimodular_matrices.h"
#include "dbl_systems_host.h"
#include "dbl_systems_kernels.h"
#include "dbl_bals_host.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_kernels.h"
#include "dbl_newton_testers.h"
#include "dbl_newton_method.h"
#include "dbl_fabry_host.h"
#include "dbl_path_tracker.h"

using namespace std;

int dbl_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrhs, double ***cff, double **acc,
   double **input_h, double **input_d, double ***output_h, double ***output_d,
   double **funval_h, double **funval_d,
   double ***jacval_h, double ***jacval_d, double **rhs_h, double **rhs_d,
   double **urhs_h, double **urhs_d, double **sol_h, double **sol_d,
   double **Q_h, double **Q_d, double **R_h, double **R_d,
   double *workvec, double **resvec, double *resmax,
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

      dbl_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,mbrhs,dpr,cff,acc,
          input_h,input_d,output_h,output_d,funval_h,funval_d,
          jacval_h,jacval_d,rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
          Q_h,Q_d,R_h,R_d,workvec,resvec,resmax,
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
                           << input_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << input_h[i][j] << endl;
         }
         if(mode == 2)
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << input_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << input_h[i][j] << endl;
            errsum += abs(input_h[i][j] - input_d[i][j]);
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

int test_dbl_real_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl )
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
      for(int j=0; j<dim; j++)
      {
         cff[i][j] = new double[degp1];
         for(int k=0; k<deg; k++) cff[i][j][k] = 0.0;
      }
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

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **startsol = new double*[dim];

   for(int i=0; i<dim; i++) startsol[i] = new double[degp1];
   make_real_exponentials(dim,deg,startsol);

   if(nbrcol != 1) // generate coefficients for the columns
   {
      // sets only the leading coefficient to a random double ...
      make_real_coefficients(nbrcol,dim,cff);
   }
   // compute the right hand sides via evaluation

   double **mbrhs = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhs[i] = new double[degp1];

      mbrhs[i][0] = 1.0;     // initialize product to one

      for(int k=1; k<degp1; k++) mbrhs[i][k] = 0.0;
   }
   if(nbrcol == 1)
      evaluate_real_monomials(dim,deg,rowsA,startsol,mbrhs);
   else
      evaluate_real_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cff,startsol,mbrhs,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++) mbrhs[i][j] = mbrhs[i][j] - mbrhs[i][j-1];

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhs[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         input_h[i][j] = startsol[i][j];
         input_d[i][j] = startsol[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << input_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   double *ratios_d = new double[dim];
   double *ratios_h = new double[dim];
   double step_d,step_h;
   double *predval_d = new double[dim];
   double *predval_h = new double[dim];
   int fail;

   dbl_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,
       nvr,idx,exp,nbrfac,expfac,mbrhs,cff,acc,
       input_h,input_d,output_h,output_d,funval_h,funval_d,
       jacval_h,jacval_d,rhs_h,rhs_d,urhs_h,urhs_d,sol_h,sol_d,
       Q_h,Q_d,R_h,R_d,workvec,resvec,&resmax,vrblvl,mode);

   if((mode == 0) || (mode == 2))
   {
      fail = dbl_fabry_step(dim,deg,input_d,ratios_d,&step_d,1); // vrblvl);
      if(step_d > 0.2) step_d = 0.2;
      fail = dbl_fabry_predictor(dim,deg,input_d,0.5*step_d,predval_d,2);
   }
   if((mode == 1) || (mode == 2))
   {
      fail = dbl_fabry_step(dim,deg,input_h,ratios_h,&step_h,1); // vrblvl);
      if(step_h > 0.2) step_h = 0.2;
      fail = dbl_fabry_predictor(dim,deg,input_h,0.5*step_h,predval_h,2);
   }
   return 0;
}
