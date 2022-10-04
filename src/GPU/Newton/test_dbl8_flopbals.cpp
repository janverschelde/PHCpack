/* Tests flops on the blocked accelerated linear series
 * in octo double precision. */

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
#include "random8_vectors.h"
#include "dbl8_tail_kernels.h"
#include "dbl_bals_flopcounts.h"

using namespace std;

int test_dbl8_flopbals ( int dim, int deg, int szt, int nbt );
/*
 * DESCRIPTION :
 *   Tests the flops counts on real random data.
 *
 * ON ENTRY :
 *   dim      dimension of the linear system;
 *   deg      degree of truncation;
 *   szt      number of threads in a block;
 *   nbt      number of blocks, szt*nbt must equal dim */

int test_cmplx8_flopbals ( int dim, int deg, int szt, int nbt );
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
        << " with octo doubles ..." << endl;

   int seed,dim,deg,szt,nbt,cdata;

   prompt_flopbals_setup(&seed,&dim,&deg,&szt,&nbt,&cdata);

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

   cout << "Launching kernels ... " << endl << endl;

   if(cdata == 0)
      test_dbl8_flopbals(dim,deg,szt,nbt);
   else
      test_cmplx8_flopbals(dim,deg,szt,nbt);

   return 0;
}

int test_dbl8_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

   double **rhshihihi = new double*[degp1];
   double **rhslohihi = new double*[degp1];
   double **rhshilohi = new double*[degp1];
   double **rhslolohi = new double*[degp1];
   double **rhshihilo = new double*[degp1];
   double **rhslohilo = new double*[degp1];
   double **rhshilolo = new double*[degp1];
   double **rhslololo = new double*[degp1];
   double **solhihihi = new double*[degp1];
   double **sollohihi = new double*[degp1];
   double **solhilohi = new double*[degp1];
   double **sollolohi = new double*[degp1];
   double **solhihilo = new double*[degp1];
   double **sollohilo = new double*[degp1];
   double **solhilolo = new double*[degp1];
   double **sollololo = new double*[degp1];
   double **reshihihi = new double*[degp1];
   double **reslohihi = new double*[degp1];
   double **reshilohi = new double*[degp1];
   double **reslolohi = new double*[degp1];
   double **reshihilo = new double*[degp1];
   double **reslohilo = new double*[degp1];
   double **reshilolo = new double*[degp1];
   double **reslololo = new double*[degp1];

   double ***mathihihi = new double**[degp1];
   double ***matlohihi = new double**[degp1];
   double ***mathilohi = new double**[degp1];
   double ***matlolohi = new double**[degp1];
   double ***mathihilo = new double**[degp1];
   double ***matlohilo = new double**[degp1];
   double ***mathilolo = new double**[degp1];
   double ***matlololo = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihihi[i] = new double[dim];
      rhslohihi[i] = new double[dim];
      rhshilohi[i] = new double[dim];
      rhslolohi[i] = new double[dim];
      rhshihilo[i] = new double[dim];
      rhslohilo[i] = new double[dim];
      rhshilolo[i] = new double[dim];
      rhslololo[i] = new double[dim];
      solhihihi[i] = new double[dim];
      sollohihi[i] = new double[dim];
      solhilohi[i] = new double[dim];
      sollolohi[i] = new double[dim];
      solhihilo[i] = new double[dim];
      sollohilo[i] = new double[dim];
      solhilolo[i] = new double[dim];
      sollololo[i] = new double[dim];
      reshihihi[i] = new double[dim];
      reslohihi[i] = new double[dim];
      reshilohi[i] = new double[dim];
      reslolohi[i] = new double[dim];
      reshihilo[i] = new double[dim];
      reslohilo[i] = new double[dim];
      reshilolo[i] = new double[dim];
      reslololo[i] = new double[dim];

      mathihihi[i] = new double*[dim];
      matlohihi[i] = new double*[dim];
      mathilohi[i] = new double*[dim];
      matlolohi[i] = new double*[dim];
      mathihilo[i] = new double*[dim];
      matlohilo[i] = new double*[dim];
      mathilolo[i] = new double*[dim];
      matlololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
         double rndhihilo,rndlohilo,rndhilolo,rndlololo;

         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
         rhshihihi[i][j] = rndhihihi;
         rhslohihi[i][j] = rndlohihi;
         rhshilohi[i][j] = rndhilohi;
         rhslolohi[i][j] = rndlolohi;
         rhshihilo[i][j] = rndhihilo;
         rhslohilo[i][j] = rndlohilo;
         rhshilolo[i][j] = rndhilolo;
         rhslololo[i][j] = rndlololo;

         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
         solhihihi[i][j] = rndhihihi;
         sollohihi[i][j] = rndlohihi;
         solhilohi[i][j] = rndhilohi;
         sollolohi[i][j] = rndlolohi;
         solhihilo[i][j] = rndhihilo;
         sollohilo[i][j] = rndlohilo;
         solhilolo[i][j] = rndhilolo;
         sollololo[i][j] = rndlololo;

         mathihihi[i][j] = new double[dim];
         matlohihi[i][j] = new double[dim];
         mathilohi[i][j] = new double[dim];
         matlolohi[i][j] = new double[dim];
         mathihilo[i][j] = new double[dim];
         matlohilo[i][j] = new double[dim];
         mathilolo[i][j] = new double[dim];
         matlololo[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            random_octo_double
               (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
            mathihihi[i][j][k] = rndhihihi;
            matlohihi[i][j][k] = rndlohihi;
            mathilohi[i][j][k] = rndhilohi;
            matlolohi[i][j][k] = rndlolohi;
            mathihilo[i][j][k] = rndhihilo;
            matlohilo[i][j][k] = rndlohilo;
            mathilolo[i][j][k] = rndhilolo;
            matlololo[i][j][k] = rndlololo;
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_dbl8_linear_residue(dim,degp1,szt,nbt,
      mathihihi,matlohihi,mathilohi,matlolohi,
      mathihilo,matlohilo,mathilolo,matlololo,
      rhshihihi,rhslohihi,rhshilohi,rhslolohi,
      rhshihilo,rhslohilo,rhshilolo,rhslololo,
      solhihihi,sollohihi,solhilohi,sollolohi,
      solhihilo,sollohilo,solhilolo,sollololo,
      reshihihi,reslohihi,reshilohi,reslolohi,
      reshihilo,reslohilo,reshilolo,reslololo,
      &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
      &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
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
        << addcnt << " x 270 " << endl;
   cout << "                Number of multiplications : "
        << mulcnt << " x 1742" << endl;
   long long int flopcnt = 270*addcnt + 1742*mulcnt;
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

int test_cmplx8_flopbals ( int dim, int deg, int szt, int nbt )
{
   const int degp1 = deg+1;

   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

   double **rhsrehihihi = new double*[degp1];
   double **rhsrelohihi = new double*[degp1];
   double **rhsrehilohi = new double*[degp1];
   double **rhsrelolohi = new double*[degp1];
   double **rhsrehihilo = new double*[degp1];
   double **rhsrelohilo = new double*[degp1];
   double **rhsrehilolo = new double*[degp1];
   double **rhsrelololo = new double*[degp1];
   double **rhsimhihihi = new double*[degp1];
   double **rhsimlohihi = new double*[degp1];
   double **rhsimhilohi = new double*[degp1];
   double **rhsimlolohi = new double*[degp1];
   double **rhsimhihilo = new double*[degp1];
   double **rhsimlohilo = new double*[degp1];
   double **rhsimhilolo = new double*[degp1];
   double **rhsimlololo = new double*[degp1];
   double **solrehihihi = new double*[degp1];
   double **solrelohihi = new double*[degp1];
   double **solrehilohi = new double*[degp1];
   double **solrelolohi = new double*[degp1];
   double **solrehihilo = new double*[degp1];
   double **solrelohilo = new double*[degp1];
   double **solrehilolo = new double*[degp1];
   double **solrelololo = new double*[degp1];
   double **solimhihihi = new double*[degp1];
   double **solimlohihi = new double*[degp1];
   double **solimhilohi = new double*[degp1];
   double **solimlolohi = new double*[degp1];
   double **solimhihilo = new double*[degp1];
   double **solimlohilo = new double*[degp1];
   double **solimhilolo = new double*[degp1];
   double **solimlololo = new double*[degp1];
   double **resrehihihi = new double*[degp1];
   double **resrelohihi = new double*[degp1];
   double **resrehilohi = new double*[degp1];
   double **resrelolohi = new double*[degp1];
   double **resrehihilo = new double*[degp1];
   double **resrelohilo = new double*[degp1];
   double **resrehilolo = new double*[degp1];
   double **resrelololo = new double*[degp1];
   double **resimhihihi = new double*[degp1];
   double **resimlohihi = new double*[degp1];
   double **resimhilohi = new double*[degp1];
   double **resimlolohi = new double*[degp1];
   double **resimhihilo = new double*[degp1];
   double **resimlohilo = new double*[degp1];
   double **resimhilolo = new double*[degp1];
   double **resimlololo = new double*[degp1];

   double ***matrehihihi = new double**[degp1];
   double ***matrelohihi = new double**[degp1];
   double ***matrehilohi = new double**[degp1];
   double ***matrelolohi = new double**[degp1];
   double ***matrehihilo = new double**[degp1];
   double ***matrelohilo = new double**[degp1];
   double ***matrehilolo = new double**[degp1];
   double ***matrelololo = new double**[degp1];
   double ***matimhihihi = new double**[degp1];
   double ***matimlohihi = new double**[degp1];
   double ***matimhilohi = new double**[degp1];
   double ***matimlolohi = new double**[degp1];
   double ***matimhihilo = new double**[degp1];
   double ***matimlohilo = new double**[degp1];
   double ***matimhilolo = new double**[degp1];
   double ***matimlololo = new double**[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsrehihihi[i] = new double[dim];
      rhsrelohihi[i] = new double[dim];
      rhsrehilohi[i] = new double[dim];
      rhsrelolohi[i] = new double[dim];
      rhsrehihilo[i] = new double[dim];
      rhsrelohilo[i] = new double[dim];
      rhsrehilolo[i] = new double[dim];
      rhsrelololo[i] = new double[dim];
      rhsimhihihi[i] = new double[dim];
      rhsimlohihi[i] = new double[dim];
      rhsimhilohi[i] = new double[dim];
      rhsimlolohi[i] = new double[dim];
      rhsimhihilo[i] = new double[dim];
      rhsimlohilo[i] = new double[dim];
      rhsimhilolo[i] = new double[dim];
      rhsimlololo[i] = new double[dim];
      solrehihihi[i] = new double[dim];
      solrelohihi[i] = new double[dim];
      solrehilohi[i] = new double[dim];
      solrelolohi[i] = new double[dim];
      solrehihilo[i] = new double[dim];
      solrelohilo[i] = new double[dim];
      solrehilolo[i] = new double[dim];
      solrelololo[i] = new double[dim];
      solimhihihi[i] = new double[dim];
      solimlohihi[i] = new double[dim];
      solimhilohi[i] = new double[dim];
      solimlolohi[i] = new double[dim];
      solimhihilo[i] = new double[dim];
      solimlohilo[i] = new double[dim];
      solimhilolo[i] = new double[dim];
      solimlololo[i] = new double[dim];
      resrehihihi[i] = new double[dim];
      resrelohihi[i] = new double[dim];
      resrehilohi[i] = new double[dim];
      resrelolohi[i] = new double[dim];
      resrehihilo[i] = new double[dim];
      resrelohilo[i] = new double[dim];
      resrehilolo[i] = new double[dim];
      resrelololo[i] = new double[dim];
      resimhihihi[i] = new double[dim];
      resimlohihi[i] = new double[dim];
      resimhilohi[i] = new double[dim];
      resimlolohi[i] = new double[dim];
      resimhihilo[i] = new double[dim];
      resimlohilo[i] = new double[dim];
      resimhilolo[i] = new double[dim];
      resimlololo[i] = new double[dim];

      matrehihihi[i] = new double*[dim];
      matrelohihi[i] = new double*[dim];
      matrehilohi[i] = new double*[dim];
      matrelolohi[i] = new double*[dim];
      matrehihilo[i] = new double*[dim];
      matrelohilo[i] = new double*[dim];
      matrehilolo[i] = new double*[dim];
      matrelololo[i] = new double*[dim];
      matimhihihi[i] = new double*[dim];
      matimlohihi[i] = new double*[dim];
      matimhilohi[i] = new double*[dim];
      matimlolohi[i] = new double*[dim];
      matimhihilo[i] = new double*[dim];
      matimlohilo[i] = new double*[dim];
      matimhilolo[i] = new double*[dim];
      matimlololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
         double rndhihilo,rndlohilo,rndhilolo,rndlololo;

         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
         rhsrehihihi[i][j] = rndhihihi;
         rhsrelohihi[i][j] = rndlohihi;
         rhsrehilohi[i][j] = rndhilohi;
         rhsrelolohi[i][j] = rndlolohi;
         rhsrehihilo[i][j] = rndhihilo;
         rhsrelohilo[i][j] = rndlohilo;
         rhsrehilolo[i][j] = rndhilolo;
         rhsrelololo[i][j] = rndlololo;
         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
         rhsimhihihi[i][j] = rndhihihi;
         rhsimlohihi[i][j] = rndlohihi;
         rhsimhilohi[i][j] = rndhilohi;
         rhsimlolohi[i][j] = rndlolohi;
         rhsimhihilo[i][j] = rndhihilo;
         rhsimlohilo[i][j] = rndlohilo;
         rhsimhilolo[i][j] = rndhilolo;
         rhsimlololo[i][j] = rndlololo;
         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
         solrehihihi[i][j] = rndhihihi;
         solrelohihi[i][j] = rndlohihi;
         solrehilohi[i][j] = rndhilohi;
         solrelolohi[i][j] = rndlolohi;
         solrehihilo[i][j] = rndhihilo;
         solrelohilo[i][j] = rndlohilo;
         solrehilolo[i][j] = rndhilolo;
         solrelololo[i][j] = rndlololo;
         random_octo_double
            (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
             &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
         solimhihihi[i][j] = rndhihihi;
         solimlohihi[i][j] = rndlohihi;
         solimhilohi[i][j] = rndhilohi;
         solimlolohi[i][j] = rndlolohi;
         solimhihilo[i][j] = rndhihilo;
         solimlohilo[i][j] = rndlohilo;
         solimhilolo[i][j] = rndhilolo;
         solimlololo[i][j] = rndlololo;

         matrehihihi[i][j] = new double[dim];
         matrelohihi[i][j] = new double[dim];
         matrehilohi[i][j] = new double[dim];
         matrelolohi[i][j] = new double[dim];
         matrehihilo[i][j] = new double[dim];
         matrelohilo[i][j] = new double[dim];
         matrehilolo[i][j] = new double[dim];
         matrelololo[i][j] = new double[dim];
         matimhihihi[i][j] = new double[dim];
         matimlohihi[i][j] = new double[dim];
         matimhilohi[i][j] = new double[dim];
         matimlolohi[i][j] = new double[dim];
         matimhihilo[i][j] = new double[dim];
         matimlohilo[i][j] = new double[dim];
         matimhilolo[i][j] = new double[dim];
         matimlololo[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            random_octo_double
               (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
            matrehihihi[i][j][k] = rndhihihi;
            matrelohihi[i][j][k] = rndlohihi;
            matrehilohi[i][j][k] = rndhilohi;
            matrelolohi[i][j][k] = rndlolohi;
            matrehihilo[i][j][k] = rndhihilo;
            matrelohilo[i][j][k] = rndlohilo;
            matrehilolo[i][j][k] = rndhilolo;
            matrelololo[i][j][k] = rndlololo;
            random_octo_double
               (&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
            matimhihihi[i][j][k] = rndhihihi;
            matimlohihi[i][j][k] = rndlohihi;
            matimhilohi[i][j][k] = rndhilohi;
            matimlolohi[i][j][k] = rndlolohi;
            matimhihilo[i][j][k] = rndhihilo;
            matimlohilo[i][j][k] = rndlohilo;
            matimhilolo[i][j][k] = rndhilolo;
            matimlololo[i][j][k] = rndlololo;
         }
      }
   }
   struct timeval begintime,endtime; // wall clock time of computations
   double elapsedms;
   long long int addcnt = 0;
   long long int mulcnt = 0;

   gettimeofday(&begintime,0);
   GPU_cmplx8_linear_residue(dim,degp1,szt,nbt,
      matrehihihi,matrelohihi,matrehilohi,matrelolohi,
      matrehihilo,matrelohilo,matrehilolo,matrelololo,
      matimhihihi,matimlohihi,matimhilohi,matimlolohi,
      matimhihilo,matimlohilo,matimhilolo,matimlololo,
      rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
      rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
      rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
      rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
      solrehihihi,solrelohihi,solrehilohi,solrelolohi,
      solrehihilo,solrelohilo,solrehilolo,solrelololo,
      solimhihihi,solimlohihi,solimhilohi,solimlolohi,
      solimhihilo,solimlohilo,solimhilolo,solimlololo,
      resrehihihi,resrelohihi,resrehilohi,resrelolohi,
      resrehihilo,resrelohilo,resrehilolo,resrelololo,
      resimhihihi,resimlohihi,resimhilohi,resimlolohi,
      resimhihilo,resimlohilo,resimhilolo,resimlololo,
      &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
      &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
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
        << addcnt << " x 270 " << endl;
   cout << "                Number of multiplications : "
        << mulcnt << " x 1742" << endl;
   long long int flopcnt = 270*addcnt + 1742*mulcnt;
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
