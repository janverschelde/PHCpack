/* Tests flops on the blocked accelerated linear series
 * in quad double precision. */

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
#include "random4_vectors.h"
#include "dbl4_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

using namespace std;

int test_dbl4_flopbals ( int dim, int deg, int szt, int nbt );
/*
 * DESCRIPTION :
 *   Tests the flops counts on real random data.
 *
 * ON ENTRY :
 *   dim      dimension of the linear system;
 *   deg      degree of truncation;
 *   szt      number of threads in a block;
 *   nbt      number of blocks, szt*nbt must equal dim */

int test_cmplx4_flopbals ( int dim, int deg, int szt, int nbt );
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
        << " with quad doubles ..." << endl;

   int seed,dim,deg,szt,nbt,cdata;

   prompt_flopbals_setup(&seed,&dim,&deg,&szt,&nbt,&cdata);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   cout << "Launching kernels ... " << endl << endl;

   if(cdata == 0)
      test_dbl4_flopbals(dim,deg,szt,nbt);
   else
      test_cmplx4_flopbals(dim,deg,szt,nbt);

   return 0;
}

int test_dbl4_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo;

   double **rhshihi = new double*[degp1];
   double **rhslohi = new double*[degp1];
   double **rhshilo = new double*[degp1];
   double **rhslolo = new double*[degp1];
   double **solhihi = new double*[degp1];
   double **sollohi = new double*[degp1];
   double **solhilo = new double*[degp1];
   double **sollolo = new double*[degp1];
   double **reshihi = new double*[degp1];
   double **reslohi = new double*[degp1];
   double **reshilo = new double*[degp1];
   double **reslolo = new double*[degp1];

   double ***mathihi = new double**[degp1];
   double ***matlohi = new double**[degp1];
   double ***mathilo = new double**[degp1];
   double ***matlolo = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihi[i] = new double[dim];
      rhslohi[i] = new double[dim];
      rhshilo[i] = new double[dim];
      rhslolo[i] = new double[dim];
      solhihi[i] = new double[dim];
      sollohi[i] = new double[dim];
      solhilo[i] = new double[dim];
      sollolo[i] = new double[dim];
      reshihi[i] = new double[dim];
      reslohi[i] = new double[dim];
      reshilo[i] = new double[dim];
      reslolo[i] = new double[dim];

      mathihi[i] = new double*[dim];
      matlohi[i] = new double*[dim];
      mathilo[i] = new double*[dim];
      matlolo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         double rndhihi,rndlohi,rndhilo,rndlolo;

         random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
         rhshihi[i][j] = rndhihi;
         rhslohi[i][j] = rndlohi;
         rhshilo[i][j] = rndhilo;
         rhslolo[i][j] = rndlolo;

         random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
         solhihi[i][j] = rndhihi;
         sollohi[i][j] = rndlohi;
         solhilo[i][j] = rndhilo;
         sollolo[i][j] = rndlolo;

         mathihi[i][j] = new double[dim];
         matlohi[i][j] = new double[dim];
         mathilo[i][j] = new double[dim];
         matlolo[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
            mathihi[i][j][k] = rndhihi;
            matlohi[i][j][k] = rndlohi;
            mathilo[i][j][k] = rndhilo;
            matlolo[i][j][k] = rndlolo;
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_dbl4_linear_residue(dim,degp1,szt,nbt,
      mathihi,matlohi,mathilo,matlolo,rhshihi,rhslohi,rhshilo,rhslolo,
      solhihi,sollohi,solhilo,sollolo,reshihi,reslohi,reshilo,reslolo,
      &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,
      &elapsedms,&addcnt,&mulcnt,0);
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   
   double walltimesec = seconds + microseconds*1.0e-6;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << walltimesec << " seconds." << endl;

   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : "
        << walltimesec << " seconds." << endl;

   cout << "                        Time spent by kernels : "
        << elapsedms << " milliseconds." << endl;

   cout << "        Number of additions/substractions : "
        << addcnt << " x 89 " << endl;
   cout << "                Number of multiplications : "
        << mulcnt << " x 336" << endl;
   long long int flopcnt = 89*addcnt + 336*mulcnt;
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

int test_cmplx4_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo;

   double **rhsrehihi = new double*[degp1];
   double **rhsrelohi = new double*[degp1];
   double **rhsrehilo = new double*[degp1];
   double **rhsrelolo = new double*[degp1];
   double **rhsimhihi = new double*[degp1];
   double **rhsimlohi = new double*[degp1];
   double **rhsimhilo = new double*[degp1];
   double **rhsimlolo = new double*[degp1];
   double **solrehihi = new double*[degp1];
   double **solrelohi = new double*[degp1];
   double **solrehilo = new double*[degp1];
   double **solrelolo = new double*[degp1];
   double **solimhihi = new double*[degp1];
   double **solimlohi = new double*[degp1];
   double **solimhilo = new double*[degp1];
   double **solimlolo = new double*[degp1];
   double **resrehihi = new double*[degp1];
   double **resrelohi = new double*[degp1];
   double **resrehilo = new double*[degp1];
   double **resrelolo = new double*[degp1];
   double **resimhihi = new double*[degp1];
   double **resimlohi = new double*[degp1];
   double **resimhilo = new double*[degp1];
   double **resimlolo = new double*[degp1];

   double ***matrehihi = new double**[degp1];
   double ***matrelohi = new double**[degp1];
   double ***matrehilo = new double**[degp1];
   double ***matrelolo = new double**[degp1];
   double ***matimhihi = new double**[degp1];
   double ***matimlohi = new double**[degp1];
   double ***matimhilo = new double**[degp1];
   double ***matimlolo = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsrehihi[i] = new double[dim];
      rhsrelohi[i] = new double[dim];
      rhsrehilo[i] = new double[dim];
      rhsrelolo[i] = new double[dim];
      rhsimhihi[i] = new double[dim];
      rhsimlohi[i] = new double[dim];
      rhsimhilo[i] = new double[dim];
      rhsimlolo[i] = new double[dim];
      solrehihi[i] = new double[dim];
      solrelohi[i] = new double[dim];
      solrehilo[i] = new double[dim];
      solrelolo[i] = new double[dim];
      solimhihi[i] = new double[dim];
      solimlohi[i] = new double[dim];
      solimhilo[i] = new double[dim];
      solimlolo[i] = new double[dim];
      resrehihi[i] = new double[dim];
      resrelohi[i] = new double[dim];
      resrehilo[i] = new double[dim];
      resrelolo[i] = new double[dim];
      resimhihi[i] = new double[dim];
      resimlohi[i] = new double[dim];
      resimhilo[i] = new double[dim];
      resimlolo[i] = new double[dim];

      matrehihi[i] = new double*[dim];
      matrelohi[i] = new double*[dim];
      matrehilo[i] = new double*[dim];
      matrelolo[i] = new double*[dim];
      matimhihi[i] = new double*[dim];
      matimlohi[i] = new double*[dim];
      matimhilo[i] = new double*[dim];
      matimlolo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         double rndhihi,rndlohi,rndhilo,rndlolo;

         random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
         rhsrehihi[i][j] = rndhihi;
         rhsrelohi[i][j] = rndlohi;
         rhsrehilo[i][j] = rndhilo;
         rhsrelolo[i][j] = rndlolo;
         random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
         rhsimhihi[i][j] = rndhihi;
         rhsimlohi[i][j] = rndlohi;
         rhsimhilo[i][j] = rndhilo;
         rhsimlolo[i][j] = rndlolo;
         random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
         solrehihi[i][j] = rndhihi;
         solrelohi[i][j] = rndlohi;
         solrehilo[i][j] = rndhilo;
         solrelolo[i][j] = rndlolo;
         random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
         solimhihi[i][j] = rndhihi;
         solimlohi[i][j] = rndlohi;
         solimhilo[i][j] = rndhilo;
         solimlolo[i][j] = rndlolo;

         matrehihi[i][j] = new double[dim];
         matrelohi[i][j] = new double[dim];
         matrehilo[i][j] = new double[dim];
         matrelolo[i][j] = new double[dim];
         matimhihi[i][j] = new double[dim];
         matimlohi[i][j] = new double[dim];
         matimhilo[i][j] = new double[dim];
         matimlolo[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
            matrehihi[i][j][k] = rndhihi;
            matrelohi[i][j][k] = rndlohi;
            matrehilo[i][j][k] = rndhilo;
            matrelolo[i][j][k] = rndlolo;
            random_quad_double(&rndhihi,&rndlohi,&rndhilo,&rndlolo);
            matimhihi[i][j][k] = rndhihi;
            matimlohi[i][j][k] = rndlohi;
            matimhilo[i][j][k] = rndhilo;
            matimlolo[i][j][k] = rndlolo;
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_cmplx4_linear_residue(dim,degp1,szt,nbt,
      matrehihi,matrelohi,matrehilo,matrelolo,
      matimhihi,matimlohi,matimhilo,matimlolo,
      rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
      rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
      solrehihi,solrelohi,solrehilo,solrelolo,
      solimhihi,solimlohi,solimhilo,solimlolo,
      resrehihi,resrelohi,resrehilo,resrelolo,
      resimhihi,resimlohi,resimhilo,resimlolo,
      &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,
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

   cout << "        Number of additions/substractions : "
        << addcnt << " x 89 " << endl;
   cout << "                Number of multiplications : "
        << mulcnt << " x 336" << endl;
   long long int flopcnt = 89*addcnt + 336*mulcnt;
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
