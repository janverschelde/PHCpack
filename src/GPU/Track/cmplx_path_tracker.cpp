// The file cmplx_path_tracker.cpp defines the functions with prototypes in
// the file cmplx_path_tracker.h.

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
#include "cmplx_newton_method.h"
#include "cmplx_path_tracker.h"

using namespace std;

int cmplx_run_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol, int nbsteps,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrhsre, double **mbrhsim,
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

      cmplx_newton_qrstep
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
          workvecre,workvecim,resvecre,resvecim,resmax,
          &zeroQ_h,&noqr_h,&zeroQ_d,&noqr_d,
          &upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
          &totcnvlapsedms,&totqrlapsedms,&totqtblapsedms,
          &totbslapsedms,&totupdlapsedms,&totreslapsedms,vrblvl,mode);

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
                           << inputre_d[i][j] << "  "
                           << inputim_d[i][j] << endl;
         }
         if(mode == 1)
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputre_h[i][j] << "  "
                           << inputim_h[i][j] << endl;
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputre_d[i][j] << "  "
                           << inputim_d[i][j] << endl;
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputre_h[i][j] << "  "
                           << inputim_h[i][j] << endl;
            errsum += abs(inputre_h[i][j] - inputre_d[i][j])
                    + abs(inputim_h[i][j] - inputim_d[i][j]);
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

int test_dbl_cmplx_track
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   int nbsteps, int mode, int vrblvl )
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
   double **accre = new double*[dim+1]; // accumulate series in one column
   double **accim = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accre[i] = new double[degp1];
      accim[i] = new double[degp1];
   }
   double ***cffre = new double**[nbrcol];
   double ***cffim = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffre[i] = new double*[dim]; // the coefficients of monomials
      cffim[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffre[i][j] = new double[degp1];
         cffim[i][j] = new double[degp1];

         for(int k=0; k<deg; k++)
         {
            cffre[i][j][k] = 0.0;
            cffim[i][j][k] = 0.0;
         }
      }
   }
   double ***outputre_h;
   double ***outputim_h;
   double ***outputre_d;
   double ***outputim_d;

   if((mode == 1) || (mode == 2))
   {
      outputre_h = new double**[dim];
      outputim_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputre_h[i] = new double*[dim+1];
         outputim_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputre_h[i][j] = new double[degp1];
            outputim_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputre_d = new double**[dim];
      outputim_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputre_d[i] = new double*[dim+1];
         outputim_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputre_d[i][j] = new double[degp1];
            outputim_d[i][j] = new double[degp1];
         }
      }
   }
   double **funvalre_h;  // function values on host
   double **funvalim_h;
   double **funvalre_d;  // function values on device
   double **funvalim_d;

   if((mode == 1) || (mode == 2))
   {
      funvalre_h = new double*[dim];
      funvalim_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalre_h[i] = new double[degp1];
         funvalim_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalre_d = new double*[dim];
      funvalim_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalre_d[i] = new double[degp1];
         funvalim_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalre_h;
   double ***jacvalim_h;
   double ***jacvalre_d;
   double ***jacvalim_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalre_h = new double**[degp1];
      jacvalim_h = new double**[degp1];

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
      jacvalre_d = new double**[degp1];
      jacvalim_d = new double**[degp1];

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
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solre_h;
   double **solim_h;
   double **solre_d;
   double **solim_d;

   if((mode == 1) || (mode == 2))
   {
      solre_h = new double*[degp1];
      solim_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solre_h[i] = new double[dim];
         solim_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solre_d = new double*[degp1];
      solim_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solre_d[i] = new double[dim];
         solim_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsre_h = new double*[degp1];
   double **rhsim_h = new double*[degp1];
   double **rhsre_d = new double*[degp1];
   double **rhsim_d = new double*[degp1];

   if((mode == 1) || (mode == 2))
   {
      rhsre_h = new double*[degp1];
      rhsim_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsre_h[i] = new double[dim];
         rhsim_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhsre_d = new double*[degp1];
      rhsim_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsre_d[i] = new double[dim];
         rhsim_d[i] = new double[dim];
      }
   }
   double *workvecre = new double[dim];
   double *workvecim = new double[dim];
   // Copy the rhs vector into work space for inplace solver.

   double **urhsre_h;
   double **urhsim_h;
   double **urhsre_d;
   double **urhsim_d;

   if((mode == 1) || (mode == 2))
   {
      urhsre_h = new double*[degp1];
      urhsim_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsre_h[i] = new double[dim];
         urhsim_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhsre_d = new double*[degp1];
      urhsim_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsre_d[i] = new double[dim];
         urhsim_d[i] = new double[dim];
      }
   }
   double **resvecre = new double*[degp1];
   double **resvecim = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecre[i] = new double[dim];
      resvecim[i] = new double[dim];
   }
   double resmax;

   double **Qre_h;
   double **Qim_h;
   double **Qre_d;
   double **Qim_d;

   if((mode == 1) || (mode == 2))
   {
      Qre_h = new double*[dim];
      Qim_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qre_h[i] = new double[dim];
         Qim_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qre_d = new double*[dim];
      Qim_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qre_d[i] = new double[dim];
         Qim_d[i] = new double[dim];
      }
   }
   double **Rre_h;
   double **Rim_h;
   double **Rre_d;
   double **Rim_d;

   if((mode == 1) || (mode == 2))
   {
      Rre_h = new double*[dim];
      Rim_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rre_h[i] = new double[dim];
         Rim_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rre_d = new double*[dim];
      Rim_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rre_d[i] = new double[dim];
         Rim_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double *angles = new double[dim];
   double **startsolre = new double*[dim];
   double **startsolim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      startsolre[i] = new double[degp1];
      startsolim[i] = new double[degp1];
   }
   make_complex_exponentials(dim,deg,angles,startsolre,startsolim);

   if(nbrcol != 1) // generate coefficients for the columns
   {
      // sets only the leading coefficient to a random double ...
      make_complex_coefficients(nbrcol,dim,cffre,cffim);
   }
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
   if(nbrcol == 1)
      evaluate_complex_monomials
         (dim,deg,rowsA,startsolre,startsolim,mbrhsre,mbrhsim);
   else
      evaluate_complex_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffre,cffim,
          startsolre,startsolim,mbrhsre,mbrhsim,vrblvl);

   // rhs coefficients are c(t) = (1-t)*c(t) = c(t) - t*c(t)
   for(int i=0; i<dim; i++)
      for(int j=1; j<degp1; j++)
      {
         mbrhsre[i][j] = mbrhsre[i][j] - mbrhsre[i][j-1];
         mbrhsim[i][j] = mbrhsim[i][j] - mbrhsim[i][j-1];
      }

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsre[i][j] << "  " << mbrhsim[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputre_h[i][j] = startsolre[i][j];
         inputim_h[i][j] = startsolim[i][j];
         inputre_d[i][j] = startsolre[i][j];
         inputim_d[i][j] = startsolim[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputre_h[i][0] << "  "
                            << inputim_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   cmplx_run_newton
      (szt,nbt,dim,deg,nbrcol,nbsteps,
       nvr,idx,exp,nbrfac,expfac,mbrhsre,mbrhsim,cffre,cffim,accre,accim,
       inputre_h,inputim_h,inputre_d,inputim_d,
       outputre_h,outputim_h,outputre_d,outputim_d,
       funvalre_h,funvalim_h,funvalre_d,funvalim_d,
       jacvalre_h,jacvalim_h,jacvalre_d,jacvalim_d,
       rhsre_h,rhsim_h,rhsre_d,rhsim_d,
       urhsre_h,urhsim_h,urhsre_d,urhsim_d,
       solre_h,solim_h,solre_d,solim_d,
       Qre_h,Qim_h,Qre_d,Qim_d,Rre_h,Rim_h,Rre_d,Rim_d,
       workvecre,workvecim,resvecre,resvecim,&resmax,vrblvl,mode);

   return 0;
}
