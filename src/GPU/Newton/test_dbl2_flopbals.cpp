/* Tests flops on the blocked accelerated linear series
 * in double double precision. */

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
#include "random2_vectors.h"
#include "dbl2_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

using namespace std;

int test_dbl2_flopbals ( int dim, int deg, int szt, int nbt );
/*
 * DESCRIPTION :
 *   Tests the flops counts on real random data.
 *
 * ON ENTRY :
 *   dim      dimension of the linear system;
 *   deg      degree of truncation;
 *   szt      number of threads in a block;
 *   nbt      number of blocks, szt*nbt must equal dim */

int test_cmplx2_flopbals ( int dim, int deg, int szt, int nbt );
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
        << " with double doubles ..." << endl;

   int seed,dim,deg,szt,nbt,cdata;

   prompt_flopbals_setup(&seed,&dim,&deg,&szt,&nbt,&cdata);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   cout << "Launching kernels ... " << endl << endl;

   if(cdata == 0)
      test_dbl2_flopbals(dim,deg,szt,nbt);
   else
      test_cmplx2_flopbals(dim,deg,szt,nbt);

   return 0;
}

int test_dbl2_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmaxhi,resmaxlo;

   double **rhshi = new double*[degp1];
   double **rhslo = new double*[degp1];
   double **solhi = new double*[degp1];
   double **sollo = new double*[degp1];
   double **reshi = new double*[degp1];
   double **reslo = new double*[degp1];

   double ***mathi = new double**[degp1];
   double ***matlo = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshi[i] = new double[dim];
      rhslo[i] = new double[dim];
      solhi[i] = new double[dim];
      sollo[i] = new double[dim];
      reshi[i] = new double[dim];
      reslo[i] = new double[dim];

      mathi[i] = new double*[dim];
      matlo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         double rndhi,rndlo;

         random_double_double(&rndhi,&rndlo);
         rhshi[i][j] = rndhi;
         rhslo[i][j] = rndlo;

         random_double_double(&rndhi,&rndlo);
         solhi[i][j] = rndhi;
         sollo[i][j] = rndlo;

         mathi[i][j] = new double[dim];
         matlo[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            random_double_double(&rndhi,&rndlo);
            mathi[i][j][k] = rndhi;
            matlo[i][j][k] = rndlo;
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_dbl2_linear_residue
      (dim,degp1,szt,nbt,mathi,matlo,rhshi,rhslo,solhi,sollo,
       reshi,reslo,&resmaxhi,&resmaxlo,&elapsedms,&addcnt,&mulcnt,0);
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   
   double walltimesec = seconds + microseconds*1.0e-6;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << walltimesec << " seconds." << endl;

   cout << "                        Time spent by kernels : "
        << elapsedms << " milliseconds." << endl;

   cout << "        Number of additions/substractions : "
        << addcnt << " x 20 " << endl;
   cout << "                Number of multiplications : "
        << mulcnt << " x 23" << endl;
   long long int flopcnt = 20*addcnt + 23*mulcnt;
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

int test_cmplx2_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmaxhi,resmaxlo;

   double **rhsrehi = new double*[degp1];
   double **rhsrelo = new double*[degp1];
   double **rhsimhi = new double*[degp1];
   double **rhsimlo = new double*[degp1];
   double **solrehi = new double*[degp1];
   double **solrelo = new double*[degp1];
   double **solimhi = new double*[degp1];
   double **solimlo = new double*[degp1];
   double **resrehi = new double*[degp1];
   double **resrelo = new double*[degp1];
   double **resimhi = new double*[degp1];
   double **resimlo = new double*[degp1];

   double ***matrehi = new double**[degp1];
   double ***matrelo = new double**[degp1];
   double ***matimhi = new double**[degp1];
   double ***matimlo = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsrehi[i] = new double[dim];
      rhsrelo[i] = new double[dim];
      rhsimhi[i] = new double[dim];
      rhsimlo[i] = new double[dim];
      solrehi[i] = new double[dim];
      solrelo[i] = new double[dim];
      solimhi[i] = new double[dim];
      solimlo[i] = new double[dim];
      resrehi[i] = new double[dim];
      resrelo[i] = new double[dim];
      resimhi[i] = new double[dim];
      resimlo[i] = new double[dim];

      matrehi[i] = new double*[dim];
      matrelo[i] = new double*[dim];
      matimhi[i] = new double*[dim];
      matimlo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         double rndhi,rndlo;

         random_double_double(&rndhi,&rndlo);
         rhsrehi[i][j] = rndhi;
         rhsrelo[i][j] = rndlo;
         random_double_double(&rndhi,&rndlo);
         rhsimhi[i][j] = rndhi;
         rhsimlo[i][j] = rndlo;
         random_double_double(&rndhi,&rndlo);
         solrehi[i][j] = rndhi;
         solrelo[i][j] = rndlo;
         random_double_double(&rndhi,&rndlo);
         solimhi[i][j] = rndhi;
         solimlo[i][j] = rndlo;

         matrehi[i][j] = new double[dim];
         matrelo[i][j] = new double[dim];
         matimhi[i][j] = new double[dim];
         matimlo[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            random_double_double(&rndhi,&rndlo);
            matrehi[i][j][k] = rndhi;
            matrelo[i][j][k] = rndlo;
            random_double_double(&rndhi,&rndlo);
            matimhi[i][j][k] = rndhi;
            matimlo[i][j][k] = rndlo;
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_cmplx2_linear_residue
      (dim,degp1,szt,nbt,matrehi,matrelo,matimhi,matimlo,
                         rhsrehi,rhsrelo,rhsimhi,rhsimlo,
                         solrehi,solrelo,solimhi,solimlo,
                         resrehi,resrelo,resimhi,resimlo,
       &resmaxhi,&resmaxlo,&elapsedms,&addcnt,&mulcnt,0);
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;

   double walltimesec = seconds + microseconds*1.0e-6;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << walltimesec << " seconds." << endl;

   cout << "                        Time spent by kernels : "
        << elapsedms << " milliseconds." << endl;

   cout << "        Number of additions/substractions : "
        << addcnt << " x 20 " << endl;
   cout << "                Number of multiplications : "
        << mulcnt << " x 23" << endl;
   long long int flopcnt = 20*addcnt + 23*mulcnt;
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
