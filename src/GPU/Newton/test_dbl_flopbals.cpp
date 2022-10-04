/* Tests flops on the blocked accelerated linear series
 * in double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <vector_types.h>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "random_numbers.h"
#include "dbl_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

using namespace std;

int test_dbl_flopbals ( int dim, int deg, int szt, int nbt );
/*
 * DESCRIPTION :
 *   Tests the flops counts on real random data.
 *
 * ON ENTRY :
 *   dim      dimension of the linear system;
 *   deg      degree of truncation;
 *   szt      number of threads in a block;
 *   nbt      number of blocks, szt*nbt must equal dim */

int test_cmplx_flopbals ( int dim, int deg, int szt, int nbt );
/*
 * DESCRIPTION :
 *   Tests the flops counts on complex random data.
 *
 * ON ENTRY :
 *   dim      dimension of the linear system;
 *   deg      degree of truncation;
 *   szt      number of threads in a block;
 *   nbt      number of blocks, szt*nbt must equal dim */

int main ( void )
{
   cout << "Testing flops of blocked accelerated linear series"
        << " in double precision ..." << endl;

   int seed,dim,deg,szt,nbt,cdata;

   prompt_flopbals_setup(&seed,&dim,&deg,&szt,&nbt,&cdata);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   cout << "Launching kernels ... " << endl << endl;

   if(cdata == 0)
      test_dbl_flopbals(dim,deg,szt,nbt);
   else
      test_cmplx_flopbals(dim,deg,szt,nbt);

   return 0;
}

int test_dbl_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmax;

   double **rhs = new double*[degp1];
   double **sol = new double*[degp1];
   double **res = new double*[degp1];

   double ***mat = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhs[i] = new double[dim];
      sol[i] = new double[dim];
      res[i] = new double[dim];

      mat[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         rhs[i][j] = random_double();
         sol[i][j] = random_double();
         mat[i][j] = new double[dim];
         for(int k=0; k<dim; k++)
            mat[i][j][k] = random_double();
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_dbl_linear_residue
      (dim,degp1,szt,nbt,mat,rhs,sol,res,&resmax,
       &elapsedms,&addcnt,&mulcnt,0);
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   
   double walltimesec = seconds + microseconds*1.0e-6;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << walltimesec << " seconds." << endl;

   cout << "                        Time spent by kernels : "
        << elapsedms << " milliseconds." << endl;

   cout << "        Number of additions/substractions : " << addcnt << endl;
   cout << "                Number of multiplications : " << mulcnt << endl;
   long long int flopcnt = addcnt + mulcnt;
   cout << "Total number of floating-point operations : " << flopcnt << endl;

   double kernflops = 1000.0*((double) flopcnt)/elapsedms;
   double wallflops = ((double) flopcnt)/walltimesec;
   const int gigacnt = pow(2.0,30);
   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;

   return 0;
}

int test_cmplx_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmax;

   double **rhsre = new double*[degp1];
   double **rhsim = new double*[degp1];
   double **solre = new double*[degp1];
   double **solim = new double*[degp1];
   double **resre = new double*[degp1];
   double **resim = new double*[degp1];

   double ***matre = new double**[degp1];
   double ***matim = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsre[i] = new double[dim];
      rhsim[i] = new double[dim];
      solre[i] = new double[dim];
      solim[i] = new double[dim];
      resre[i] = new double[dim];
      resim[i] = new double[dim];

      matre[i] = new double*[dim];
      matim[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         rhsre[i][j] = random_double();
         rhsim[i][j] = random_double();
         solre[i][j] = random_double();
         solim[i][j] = random_double();
         matre[i][j] = new double[dim];
         matim[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            matre[i][j][k] = random_double();
            matim[i][j][k] = random_double();
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_cmplx_linear_residue
      (dim,degp1,szt,nbt,matre,matim,rhsre,rhsim,solre,solim,
       resre,resim,&resmax,&elapsedms,&addcnt,&mulcnt,0);
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;

   double walltimesec = seconds + microseconds*1.0e-6;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << walltimesec << " seconds." << endl;

   cout << "                        Time spent by kernels : "
        << elapsedms << " milliseconds." << endl;

   cout << "        Number of additions/substractions : " << addcnt << endl;
   cout << "                Number of multiplications : " << mulcnt << endl;
   long long int flopcnt = addcnt + mulcnt;
   cout << "Total number of floating-point operations : " << flopcnt << endl;

   double kernflops = 1000.0*((double) flopcnt)/elapsedms;
   double wallflops = ((double) flopcnt)/walltimesec;
   const int gigacnt = pow(2.0,30);
   cout << "Kernel Time Flops : "
        << scientific << setprecision(3) << kernflops;
   cout << fixed << setprecision(3)
        << " = " << kernflops/gigacnt << " Gigaflops" << endl;
   cout << " Wall Clock Flops : "
        << scientific << setprecision(3) << wallflops;
   cout << fixed << setprecision(3)
        << " = " << wallflops/gigacnt << " Gigaflops" << endl;
 
   return 0;
}
