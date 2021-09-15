/* The file dbl8_tabs_kernels.cu defines the functions specified in
 * the file dbl8_tabs_kernels.h. */

#include <iostream>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#include "octo_double_gpufun.cu"
#endif
#include "dbl8_tabs_kernels.h"
#include "dbl_tabs_flopcounts.h"

using namespace std;

__global__ void dbl8_small_invert_upper 
( int dim,
  double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi,
  double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo,
  double *invUhihihi, double *invUlohihi,
  double *invUhilohi, double *invUlolohi,
  double *invUhihilo, double *invUlohilo,
  double *invUhilolo, double *invUlololo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolhihihi[tabsod_shmemsize];
   __shared__ double Ucollohihi[tabsod_shmemsize];
   __shared__ double Ucolhilohi[tabsod_shmemsize];
   __shared__ double Ucollolohi[tabsod_shmemsize];
   __shared__ double Ucolhihilo[tabsod_shmemsize];
   __shared__ double Ucollohilo[tabsod_shmemsize];
   __shared__ double Ucolhilolo[tabsod_shmemsize];
   __shared__ double Ucollololo[tabsod_shmemsize];
   __shared__ double invUrowshihihi[tabsod_shmemsize];
   __shared__ double invUrowslohihi[tabsod_shmemsize];
   __shared__ double invUrowshilohi[tabsod_shmemsize];
   __shared__ double invUrowslolohi[tabsod_shmemsize];
   __shared__ double invUrowshihilo[tabsod_shmemsize];
   __shared__ double invUrowslohilo[tabsod_shmemsize];
   __shared__ double invUrowshilolo[tabsod_shmemsize];
   __shared__ double invUrowslololo[tabsod_shmemsize];

   double rhshihihi,rhslohihi,rhshilohi,rhslolohi;
   double rhshihilo,rhslohilo,rhshilolo,rhslololo;
   double xvalhihihi,xvallohihi,xvalhilohi,xvallolohi;
   double xvalhihilo,xvallohilo,xvalhilolo,xvallololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   int colidx = dim*(dim-1);          // start with the last column

   Ucolhihihi[k] = Uhihihi[colidx+k]; // load the last column
   Ucollohihi[k] = Ulohihi[colidx+k];
   Ucolhilohi[k] = Uhilohi[colidx+k];
   Ucollolohi[k] = Ulolohi[colidx+k];
   Ucolhihilo[k] = Uhihilo[colidx+k];
   Ucollohilo[k] = Ulohilo[colidx+k];
   Ucolhilolo[k] = Uhilolo[colidx+k];
   Ucollololo[k] = Ulololo[colidx+k];

   rhshihihi = ((double) int(k == dim-1)); // right hand side for each thread
   rhslohihi = 0.0;
   rhshilohi = 0.0;
   rhslolohi = 0.0;
   rhshihilo = 0.0;
   rhslohilo = 0.0;
   rhshilolo = 0.0;
   rhslololo = 0.0;
   int rowidx = (dim - 1)*dim + k;      // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k]; // last row of the inverse
   odg_div(rhshihihi,rhslohihi,rhshilohi,rhslolohi,
           rhshihilo,rhslohilo,rhshilolo,rhslololo,
           Ucolhihihi[k],Ucollohihi[k],Ucolhilohi[k],Ucollolohi[k],
           Ucolhihilo[k],Ucollohilo[k],Ucolhilolo[k],Ucollololo[k],
           &invUrowshihihi[rowidx],&invUrowslohihi[rowidx],
           &invUrowshilohi[rowidx],&invUrowslolohi[rowidx],
           &invUrowshihilo[rowidx],&invUrowslohilo[rowidx],
           &invUrowshilolo[rowidx],&invUrowslololo[rowidx]);

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshihihi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslohihi = 0.0;
      rhshilohi = 0.0;
      rhslolohi = 0.0;
      rhshihilo = 0.0;
      rhslohilo = 0.0;
      rhshilolo = 0.0;
      rhslololo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolhihihi[k] = Uhihihi[colidx+k];
         Ucollohihi[k] = Ulohihi[colidx+k];
         Ucolhilohi[k] = Uhilohi[colidx+k];
         Ucollolohi[k] = Ulolohi[colidx+k];
         Ucolhihilo[k] = Uhihilo[colidx+k];
         Ucollohilo[k] = Ulohilo[colidx+k];
         Ucolhilolo[k] = Uhilolo[colidx+k];
         Ucollololo[k] = Ulololo[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalhihihi = invUrowshihihi[rowidx];
         xvallohihi = invUrowslohihi[rowidx];
         xvalhilohi = invUrowshilohi[rowidx];
         xvallolohi = invUrowslolohi[rowidx];
         xvalhihilo = invUrowshihilo[rowidx];
         xvallohilo = invUrowslohilo[rowidx];
         xvalhilolo = invUrowshilolo[rowidx];
         xvallololo = invUrowslololo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         odg_mul(Ucolhihihi[i],Ucollohihi[i],Ucolhilohi[i],Ucollolohi[i],
                 Ucolhihilo[i],Ucollohilo[i],Ucolhilolo[i],Ucollololo[i],
                 xvalhihihi,   xvallohihi,   xvalhilohi,   xvallolohi,
                 xvalhihilo,   xvallohilo,   xvalhilolo,   xvallololo,
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         __syncthreads();
         odg_dec(&rhshihihi,&rhslohihi,&rhshilohi,&rhslolohi,
                 &rhshihilo,&rhslohilo,&rhshilolo,&rhslololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolhihihi[k] = Uhihihi[colidx+k];
      Ucolhilohi[k] = Uhilohi[colidx+k];
      Ucollohihi[k] = Ulohihi[colidx+k];
      Ucollolohi[k] = Ulolohi[colidx+k];
      Ucolhihilo[k] = Uhihilo[colidx+k];
      Ucolhilolo[k] = Uhilolo[colidx+k];
      Ucollohilo[k] = Ulohilo[colidx+k];
      Ucollololo[k] = Ulololo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      odg_div(rhshihihi,rhslohihi,rhshilohi,rhslolohi,
              rhshihilo,rhslohilo,rhshilolo,rhslololo,
              Ucolhihihi[i],Ucollohihi[i],Ucolhilohi[i],Ucollolohi[i],
              Ucolhihilo[i],Ucollohilo[i],Ucolhilolo[i],Ucollololo[i],
              &invUrowshihihi[rowidx],&invUrowslohihi[rowidx],
              &invUrowshilohi[rowidx],&invUrowslolohi[rowidx],
              &invUrowshihilo[rowidx],&invUrowslohilo[rowidx],
              &invUrowshilolo[rowidx],&invUrowslololo[rowidx]);
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUhihihi[rowidx+k] = invUrowshihihi[rowidx+k];
      invUlohihi[rowidx+k] = invUrowslohihi[rowidx+k];
      invUhilohi[rowidx+k] = invUrowshilohi[rowidx+k];
      invUlolohi[rowidx+k] = invUrowslolohi[rowidx+k];
      invUhihilo[rowidx+k] = invUrowshihilo[rowidx+k];
      invUlohilo[rowidx+k] = invUrowslohilo[rowidx+k];
      invUhilolo[rowidx+k] = invUrowshilolo[rowidx+k];
      invUlololo[rowidx+k] = invUrowslololo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void cmplx8_small_invert_upper
 ( int dim,
   double *Urehihihi, double *Urelohihi, double *Urehilohi, double *Urelolohi,
   double *Urehihilo, double *Urelohilo, double *Urehilolo, double *Urelololo,
   double *Uimhihihi, double *Uimlohihi, double *Uimhilohi, double *Uimlolohi,
   double *Uimhihilo, double *Uimlohilo, double *Uimhilolo, double *Uimlololo,
   double *invUrehihihi, double *invUrelohihi,
   double *invUrehilohi, double *invUrelolohi,
   double *invUrehihilo, double *invUrelohilo,
   double *invUrehilolo, double *invUrelololo,
   double *invUimhihihi, double *invUimlohihi,
   double *invUimhilohi, double *invUimlolohi,
   double *invUimhihilo, double *invUimlohilo,
   double *invUimhilolo, double *invUimlololo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolrehihihi[tabsod_shmemsize];
   __shared__ double Ucolrelohihi[tabsod_shmemsize];
   __shared__ double Ucolrehilohi[tabsod_shmemsize];
   __shared__ double Ucolrelolohi[tabsod_shmemsize];
   __shared__ double Ucolrehihilo[tabsod_shmemsize];
   __shared__ double Ucolrelohilo[tabsod_shmemsize];
   __shared__ double Ucolrehilolo[tabsod_shmemsize];
   __shared__ double Ucolrelololo[tabsod_shmemsize];
   __shared__ double Ucolimhihihi[tabsod_shmemsize];
   __shared__ double Ucolimlohihi[tabsod_shmemsize];
   __shared__ double Ucolimhilohi[tabsod_shmemsize];
   __shared__ double Ucolimlolohi[tabsod_shmemsize];
   __shared__ double Ucolimhihilo[tabsod_shmemsize];
   __shared__ double Ucolimlohilo[tabsod_shmemsize];
   __shared__ double Ucolimhilolo[tabsod_shmemsize];
   __shared__ double Ucolimlololo[tabsod_shmemsize];
   __shared__ double invUrowsrehihihi[tabsod_shmemsize];
   __shared__ double invUrowsrelohihi[tabsod_shmemsize];
   __shared__ double invUrowsrehilohi[tabsod_shmemsize];
   __shared__ double invUrowsrelolohi[tabsod_shmemsize];
   __shared__ double invUrowsrehihilo[tabsod_shmemsize];
   __shared__ double invUrowsrelohilo[tabsod_shmemsize];
   __shared__ double invUrowsrehilolo[tabsod_shmemsize];
   __shared__ double invUrowsrelololo[tabsod_shmemsize];
   __shared__ double invUrowsimhihihi[tabsod_shmemsize];
   __shared__ double invUrowsimlohihi[tabsod_shmemsize];
   __shared__ double invUrowsimhilohi[tabsod_shmemsize];
   __shared__ double invUrowsimlolohi[tabsod_shmemsize];
   __shared__ double invUrowsimhihilo[tabsod_shmemsize];
   __shared__ double invUrowsimlohilo[tabsod_shmemsize];
   __shared__ double invUrowsimhilolo[tabsod_shmemsize];
   __shared__ double invUrowsimlololo[tabsod_shmemsize];

   double rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi;
   double rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo;
   double rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi;
   double rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo;
   double xvalrehihihi,xvalrelohihi,xvalrehilohi,xvalrelolohi;
   double xvalrehihilo,xvalrelohilo,xvalrehilolo,xvalrelololo;
   double xvalimhihihi,xvalimlohihi,xvalimhilohi,xvalimlolohi;
   double xvalimhihilo,xvalimlohilo,xvalimhilolo,xvalimlololo;
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;
   double acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi;
   double acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo;
   double acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi;
   double acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo;
   double invrehihihi,invrelohihi,invrehilohi,invrelolohi;
   double invrehihilo,invrelohilo,invrehilolo,invrelololo;
   double invimhihihi,invimlohihi,invimhilohi,invimlolohi;
   double invimhihilo,invimlohilo,invimhilolo,invimlololo;
   double denhihihi,denlohihi,denhilohi,denlolohi;
   double denhihilo,denlohilo,denhilolo,denlololo;

   int colidx = dim*(dim-1);               // start with the last column

   Ucolrehihihi[k] = Urehihihi[colidx+k];  // load the last column
   Ucolrelohihi[k] = Urelohihi[colidx+k];
   Ucolrehilohi[k] = Urehilohi[colidx+k];
   Ucolrelolohi[k] = Urelolohi[colidx+k];
   Ucolrehihilo[k] = Urehihilo[colidx+k];
   Ucolrelohilo[k] = Urelohilo[colidx+k];
   Ucolrehilolo[k] = Urehilolo[colidx+k];
   Ucolrelololo[k] = Urelololo[colidx+k];
   Ucolimhihihi[k] = Uimhihihi[colidx+k];
   Ucolimlohihi[k] = Uimlohihi[colidx+k];
   Ucolimhilohi[k] = Uimhilohi[colidx+k];
   Ucolimlolohi[k] = Uimlolohi[colidx+k];
   Ucolimhihilo[k] = Uimhihilo[colidx+k];
   Ucolimlohilo[k] = Uimlohilo[colidx+k];
   Ucolimhilolo[k] = Uimhilolo[colidx+k];
   Ucolimlololo[k] = Uimlololo[colidx+k];
   rhsrehihihi = ((double) int(k == dim-1)); // right hand side for threads
   rhsrelohihi = 0.0;
   rhsrehilohi = 0.0;
   rhsrelolohi = 0.0;
   rhsrehihilo = 0.0;
   rhsrelohilo = 0.0;
   rhsrehilolo = 0.0;
   rhsrelololo = 0.0;
   rhsimhihihi = 0.0;
   rhsimlohihi = 0.0;
   rhsimhilohi = 0.0;
   rhsimlolohi = 0.0;
   rhsimhihilo = 0.0;
   rhsimlohilo = 0.0;
   rhsimhilolo = 0.0;
   rhsimlololo = 0.0;
   int rowidx = (dim - 1)*dim + k;       // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse
   odg_mul(Ucolrehihihi[k],Ucolrelohihi[k],Ucolrehilohi[k],Ucolrelolohi[k],
           Ucolrehihilo[k],Ucolrelohilo[k],Ucolrehilolo[k],Ucolrelololo[k],
           Ucolrehihihi[k],Ucolrelohihi[k],Ucolrehilohi[k],Ucolrelolohi[k],
           Ucolrehihilo[k],Ucolrelohilo[k],Ucolrehilolo[k],Ucolrelololo[k],
             &denhihihi,     &denlohihi,     &denhilohi,     &denlolohi,
             &denhihilo,     &denlohilo,     &denhilolo,     &denlololo);
   odg_mul(Ucolimhihihi[k],Ucolimlohihi[k],Ucolimhilohi[k],Ucolimlolohi[k],
           Ucolimhihilo[k],Ucolimlohilo[k],Ucolimhilolo[k],Ucolimlololo[k],
           Ucolimhihihi[k],Ucolimlohihi[k],Ucolimhilohi[k],Ucolimlolohi[k],
           Ucolimhihilo[k],Ucolimlohilo[k],Ucolimhilolo[k],Ucolimlololo[k],
            &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
            &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
   odg_inc(&denhihihi,&denlohihi,&denhilohi,&denlolohi,
           &denhihilo,&denlohilo,&denhilolo,&denlololo,
           acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
           acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);
   odg_div(Ucolrehihihi[k],Ucolrelohihi[k],Ucolrehilohi[k],Ucolrelolohi[k],
           Ucolrehihilo[k],Ucolrelohilo[k],Ucolrehilolo[k],Ucolrelololo[k],
              denhihihi,      denlohihi,      denhilohi,      denlolohi,
              denhihilo,      denlohilo,      denhilolo,      denlololo,
           &invrehihihi,   &invrelohihi,   &invrehilohi,   &invrelolohi,
           &invrehihilo,   &invrelohilo,   &invrehilolo,   &invrelololo);
   odg_div(Ucolimhihihi[k],Ucolimlohihi[k],Ucolimhilohi[k],Ucolimlolohi[k],
           Ucolimhihilo[k],Ucolimlohilo[k],Ucolimhilolo[k],Ucolimlololo[k],
              denhihihi,      denlohihi,      denhilohi,      denlolohi,
              denhihilo,      denlohilo,      denhilolo,      denlololo,
           &invimhihihi,   &invimlohihi,   &invimhilohi,   &invimlolohi,
           &invimhihilo,   &invimlohilo,   &invimhilolo,   &invimlololo);
   odg_minus(&invimhihihi,&invimlohihi,&invimhilohi,&invimlolohi,
             &invimhihilo,&invimlohilo,&invimhilolo,&invimlololo);
   odg_mul(rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
           rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
           invrehihihi,invrelohihi,invrehilohi,invrelolohi,
           invrehihilo,invrelohilo,invrehilolo,invrelololo,
           &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
           &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
   odg_mul(rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
           rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
           invimhihihi,invimlohihi,invimhilohi,invimlolohi,
           invimhihilo,invimlohilo,invimhilolo,invimlololo,
           &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
           &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
   odg_mul(rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
           rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
           invrehihihi,invrelohihi,invrehilohi,invrelolohi,
           invrehihilo,invrelohilo,invrehilolo,invrelololo,
           &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
           &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
   odg_mul(rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
           rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
           invimhihihi,invimlohihi,invimhilohi,invimlolohi,
           invimhihilo,invimlohilo,invimhilolo,invimlololo,
           &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
           &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);
   odg_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
           &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
            acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
            acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
   invUrowsrehihihi[rowidx] = acc1hihihi;
   invUrowsrelohihi[rowidx] = acc1lohihi;
   invUrowsrehilohi[rowidx] = acc1hilohi;
   invUrowsrelolohi[rowidx] = acc1lolohi;
   invUrowsrehihilo[rowidx] = acc1hihilo;
   invUrowsrelohilo[rowidx] = acc1lohilo;
   invUrowsrehilolo[rowidx] = acc1hilolo;
   invUrowsrelololo[rowidx] = acc1lololo;
   odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
           &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
            acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
            acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
   invUrowsimhihihi[rowidx] = acc3hihihi;
   invUrowsimlohihi[rowidx] = acc3lohihi;
   invUrowsimhilohi[rowidx] = acc3hilohi;
   invUrowsimlolohi[rowidx] = acc3lolohi;
   invUrowsimhihilo[rowidx] = acc3hihilo;
   invUrowsimlohilo[rowidx] = acc3lohilo;
   invUrowsimhilolo[rowidx] = acc3hilolo;
   invUrowsimlololo[rowidx] = acc3lololo;

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehihihi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhsrelohihi = 0.0;
      rhsrehilohi = 0.0;
      rhsrelolohi = 0.0;
      rhsrehihilo = 0.0;
      rhsrelohilo = 0.0;
      rhsrehilolo = 0.0;
      rhsrelololo = 0.0;
      rhsimhihihi = 0.0;
      rhsimlohihi = 0.0;
      rhsimhilohi = 0.0;
      rhsimlolohi = 0.0;
      rhsimhihilo = 0.0;
      rhsimlohilo = 0.0;
      rhsimhilolo = 0.0;
      rhsimlololo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolrehihihi[k] = Urehihihi[colidx+k];
         Ucolrelohihi[k] = Urelohihi[colidx+k];
         Ucolrehilohi[k] = Urehilohi[colidx+k];
         Ucolrelolohi[k] = Urelolohi[colidx+k];
         Ucolrehihilo[k] = Urehihilo[colidx+k];
         Ucolrelohilo[k] = Urelohilo[colidx+k];
         Ucolrehilolo[k] = Urehilolo[colidx+k];
         Ucolrelololo[k] = Urelololo[colidx+k];
         Ucolimhihihi[k] = Uimhihihi[colidx+k];
         Ucolimlohihi[k] = Uimlohihi[colidx+k];
         Ucolimhilohi[k] = Uimhilohi[colidx+k];
         Ucolimlolohi[k] = Uimlolohi[colidx+k];
         Ucolimhihilo[k] = Uimhihilo[colidx+k];
         Ucolimlohilo[k] = Uimlohilo[colidx+k];
         Ucolimhilolo[k] = Uimhilolo[colidx+k];
         Ucolimlololo[k] = Uimlololo[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalrehihihi = invUrowsrehihihi[rowidx];
         xvalrelohihi = invUrowsrelohihi[rowidx];
         xvalrehilohi = invUrowsrehilohi[rowidx];
         xvalrelolohi = invUrowsrelolohi[rowidx];
         xvalrehihilo = invUrowsrehihilo[rowidx];
         xvalrelohilo = invUrowsrelohilo[rowidx];
         xvalrehilolo = invUrowsrehilolo[rowidx];
         xvalrelololo = invUrowsrelololo[rowidx];
         xvalimhihihi = invUrowsimhihihi[rowidx];
         xvalimlohihi = invUrowsimlohihi[rowidx];
         xvalimhilohi = invUrowsimhilohi[rowidx];
         xvalimlolohi = invUrowsimlolohi[rowidx];
         xvalimhihilo = invUrowsimhihilo[rowidx];
         xvalimlohilo = invUrowsimlohilo[rowidx];
         xvalimhilolo = invUrowsimhilolo[rowidx];
         xvalimlololo = invUrowsimlololo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                 xvalrehihihi,   xvalrelohihi,  xvalrehilohi,  xvalrelolohi,
                 xvalrehihilo,   xvalrelohilo,  xvalrehilolo,  xvalrelololo,
                  &acc1hihihi,    &acc1lohihi,   &acc1hilohi,   &acc1lolohi,
                  &acc1hihilo,    &acc1lohilo,   &acc1hilolo,   &acc1lololo);
         odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
                 xvalimhihihi,   xvalimlohihi,  xvalimhilohi,  xvalimlolohi,
                 xvalimhihilo,   xvalimlohilo,  xvalimhilolo,  xvalimlololo,
                  &acc2hihihi,    &acc2lohihi,   &acc2hilohi,   &acc2lolohi,
                  &acc2hihilo,    &acc2lohilo,   &acc2hilolo,   &acc2lololo);
         odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
                 xvalrehihihi,   xvalrelohihi,  xvalrehilohi,  xvalrelolohi,
                 xvalrehihilo,   xvalrelohilo,  xvalrehilolo,  xvalrelololo,
                  &acc3hihihi,    &acc3lohihi,   &acc3hilohi,   &acc3lolohi,
                  &acc3hihilo,    &acc3lohilo,   &acc3hilolo,   &acc3lololo);
         odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                 xvalimhihihi,   xvalimlohihi,  xvalimhilohi,  xvalimlolohi,
                 xvalimhihilo,   xvalimlohilo,  xvalimhilolo,  xvalimlololo,
                  &acc4hihihi,    &acc4lohihi,   &acc4hilohi,   &acc4lolohi,
                  &acc4hihilo,    &acc4lohilo,   &acc4hilolo,   &acc4lololo);
         odg_dec(&rhsrehihihi,&rhsrelohihi,&rhsrehilohi,&rhsrelolohi,
                 &rhsrehihilo,&rhsrelohilo,&rhsrehilolo,&rhsrelololo,
                   acc1hihihi,  acc1lohihi,  acc1hilohi,  acc1lolohi,
                   acc1hihilo,  acc1lohilo,  acc1hilolo,  acc1lololo);
         odg_inc(&rhsrehihihi,&rhsrelohihi,&rhsrehilohi,&rhsrelolohi,
                 &rhsrehihilo,&rhsrelohilo,&rhsrehilolo,&rhsrelololo,
                   acc2hihihi,  acc2lohihi,  acc2hilohi,  acc2lolohi,
                   acc2hihilo,  acc2lohilo,  acc2hilolo,  acc2lololo);
         odg_dec(&rhsimhihihi,&rhsimlohihi,&rhsimhilohi,&rhsimlolohi,
                 &rhsimhihilo,&rhsimlohilo,&rhsimhilolo,&rhsimlololo,
                   acc3hihihi,  acc3lohihi,  acc3hilohi,  acc3lolohi,
                   acc3hihilo,  acc3lohilo,  acc3hilolo,  acc3lololo);
         odg_dec(&rhsimhihihi,&rhsimlohihi,&rhsimhilohi,&rhsimlolohi,
                 &rhsimhihilo,&rhsimlohilo,&rhsimhilolo,&rhsimlololo,
                   acc4hihihi,  acc4lohihi,  acc4hilohi,  acc4lolohi,
                   acc4hihilo,  acc4lohilo,  acc4hilolo,  acc4lololo);
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolrehihihi[k] = Urehihihi[colidx+k];
      Ucolrelohihi[k] = Urelohihi[colidx+k];
      Ucolrehilohi[k] = Urehilohi[colidx+k];
      Ucolrelolohi[k] = Urelolohi[colidx+k];
      Ucolrehihilo[k] = Urehihilo[colidx+k];
      Ucolrelohilo[k] = Urelohilo[colidx+k];
      Ucolrehilolo[k] = Urehilolo[colidx+k];
      Ucolrelololo[k] = Urelololo[colidx+k];
      Ucolimhihihi[k] = Uimhihihi[colidx+k];
      Ucolimlohihi[k] = Uimlohihi[colidx+k];
      Ucolimhihihi[k] = Uimhilohi[colidx+k];
      Ucolimlolohi[k] = Uimlolohi[colidx+k];
      Ucolimhihilo[k] = Uimhihilo[colidx+k];
      Ucolimlohilo[k] = Uimlohilo[colidx+k];
      Ucolimhihilo[k] = Uimhilolo[colidx+k];
      Ucolimlololo[k] = Uimlololo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],Ucolrehilohi[i],Ucolrelolohi[i],
              Ucolrehihilo[i],Ucolrelohilo[i],Ucolrehilolo[i],Ucolrelololo[i],
              Ucolrehihihi[i],Ucolrelohihi[i],Ucolrehilohi[i],Ucolrelolohi[i],
              Ucolrehihilo[i],Ucolrelohilo[i],Ucolrehilolo[i],Ucolrelololo[i],
                &denhihihi,     &denlohihi,     &denhilohi,     &denlolohi,
                &denhihilo,     &denlohilo,     &denhilolo,     &denlololo);
      odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],Ucolimhilohi[i],Ucolimlolohi[i],
              Ucolimhihilo[i],Ucolimlohilo[i],Ucolimhilolo[i],Ucolimlololo[i],
              Ucolimhihihi[i],Ucolimlohihi[i],Ucolimhilohi[i],Ucolimlolohi[i],
              Ucolimhihilo[i],Ucolimlohilo[i],Ucolimhilolo[i],Ucolimlololo[i],
               &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
               &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
      __syncthreads();
      odg_inc(&denhihihi,&denlohihi,&denhilohi,&denlolohi,
              &denhihilo,&denlohilo,&denhilolo,&denlololo,
              acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
              acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);
      odg_div(Ucolrehihihi[i],Ucolrelohihi[i],Ucolrehilohi[i],Ucolrelolohi[i],
              Ucolrehihilo[i],Ucolrelohilo[i],Ucolrehilolo[i],Ucolrelololo[i],
                 denhihihi,      denlohihi,      denhilohi,      denlolohi,
                 denhihilo,      denlohilo,      denhilolo,      denlololo,
              &invrehihihi,   &invrelohihi,   &invrehilohi,   &invrelolohi,
              &invrehihilo,   &invrelohilo,   &invrehilolo,   &invrelololo);
      odg_div(Ucolimhihihi[i],Ucolimlohihi[i],Ucolimhilohi[i],Ucolimlolohi[i],
              Ucolimhihilo[i],Ucolimlohilo[i],Ucolimhilolo[i],Ucolimlololo[i],
                 denhihihi,      denlohihi,      denhilohi,      denlolohi,
                 denhihilo,      denlohilo,      denhilolo,      denlololo,
              &invimhihihi,   &invimlohihi,   &invimhilohi,   &invimlolohi,
              &invimhihilo,   &invimlohilo,   &invimhilolo,   &invimlololo);
      odg_minus(&invimhihihi,&invimlohihi,&invimhilohi,&invimlolohi,
                &invimhihilo,&invimlohilo,&invimhilolo,&invimlololo);
      odg_mul(rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
              rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
              invrehihihi,invrelohihi,invrehilohi,invrelolohi,
              invrehihilo,invrelohilo,invrehilolo,invrelololo,
              &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
              &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
      odg_mul(rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
              rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
              invimhihihi,invimlohihi,invimhilohi,invimlolohi,
              invimhihilo,invimlohilo,invimhilolo,invimlololo,
              &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
              &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
      odg_mul(rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
              rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
              invrehihihi,invrelohihi,invrehilohi,invrelolohi,
              invrehihilo,invrelohilo,invrehilolo,invrelololo,
              &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
              &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
      odg_mul(rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
              rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
              invimhihihi,invimlohihi,invimhilohi,invimlolohi,
              invimhihilo,invimlohilo,invimhilolo,invimlololo,
              &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
              &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);
      __syncthreads();
      odg_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
              &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
               acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
               acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
      invUrowsrehihihi[rowidx] = acc1hihihi;
      invUrowsrelohihi[rowidx] = acc1lohihi;
      invUrowsrehilohi[rowidx] = acc1hilohi;
      invUrowsrelolohi[rowidx] = acc1lolohi;
      invUrowsrehihilo[rowidx] = acc1hihilo;
      invUrowsrelohilo[rowidx] = acc1lohilo;
      invUrowsrehilolo[rowidx] = acc1hilolo;
      invUrowsrelololo[rowidx] = acc1lololo;
      odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
              &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
               acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
               acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
      invUrowsimhihihi[rowidx] = acc3hihihi;
      invUrowsimlohihi[rowidx] = acc3lohihi;
      invUrowsimhilohi[rowidx] = acc3hilohi;
      invUrowsimlolohi[rowidx] = acc3lolohi;
      invUrowsimhihilo[rowidx] = acc3hihilo;
      invUrowsimlohilo[rowidx] = acc3lohilo;
      invUrowsimhilolo[rowidx] = acc3hilolo;
      invUrowsimlololo[rowidx] = acc3lololo;
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUrehihihi[rowidx+k] = invUrowsrehihihi[rowidx+k];
      invUrelohihi[rowidx+k] = invUrowsrelohihi[rowidx+k];
      invUrehilohi[rowidx+k] = invUrowsrehilohi[rowidx+k];
      invUrelolohi[rowidx+k] = invUrowsrelolohi[rowidx+k];
      invUrehihilo[rowidx+k] = invUrowsrehihilo[rowidx+k];
      invUrelohilo[rowidx+k] = invUrowsrelohilo[rowidx+k];
      invUrehilolo[rowidx+k] = invUrowsrehilolo[rowidx+k];
      invUrelololo[rowidx+k] = invUrowsrelololo[rowidx+k];
      invUimhihihi[rowidx+k] = invUrowsimhihihi[rowidx+k];
      invUimlohihi[rowidx+k] = invUrowsimlohihi[rowidx+k];
      invUimhilohi[rowidx+k] = invUrowsimhilohi[rowidx+k];
      invUimlolohi[rowidx+k] = invUrowsimlolohi[rowidx+k];
      invUimhihilo[rowidx+k] = invUrowsimhihilo[rowidx+k];
      invUimlohilo[rowidx+k] = invUrowsimlohilo[rowidx+k];
      invUimhilolo[rowidx+k] = invUrowsimhilolo[rowidx+k];
      invUimlololo[rowidx+k] = invUrowsimlololo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void dbl8_medium_invert_upper
 ( int dim,
   double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi,
   double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo,
   double *invUhihihi, double *invUlohihi,
   double *invUhilohi, double *invUlolohi,
   double *invUhihilo, double *invUlohilo,
   double *invUhilolo, double *invUlololo )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolhihihi[tabsod_shmemsize];    // one column of U
   __shared__ double Ucollohihi[tabsod_shmemsize];
   __shared__ double Ucolhilohi[tabsod_shmemsize];
   __shared__ double Ucollolohi[tabsod_shmemsize];
   __shared__ double Ucolhihilo[tabsod_shmemsize];
   __shared__ double Ucollohilo[tabsod_shmemsize];
   __shared__ double Ucolhilolo[tabsod_shmemsize];
   __shared__ double Ucollololo[tabsod_shmemsize];
   __shared__ double invUrowhihihi[tabsod_shmemsize]; // one row of invU
   __shared__ double invUrowlohihi[tabsod_shmemsize];
   __shared__ double invUrowhilohi[tabsod_shmemsize]; 
   __shared__ double invUrowlolohi[tabsod_shmemsize];
   __shared__ double invUrowhihilo[tabsod_shmemsize];
   __shared__ double invUrowlohilo[tabsod_shmemsize];
   __shared__ double invUrowhilolo[tabsod_shmemsize]; 
   __shared__ double invUrowlololo[tabsod_shmemsize];

   double rhshihihi,rhslohihi,rhshilohi,rhslolohi;
   double rhshihilo,rhslohilo,rhshilolo,rhslololo;
   double xvalhihihi,xvallohihi,xvalhilohi,xvallolohi;
   double xvalhihilo,xvallohilo,xvalhilolo,xvallololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolhihihi[k] = Uhihihi[colidx+k];  // load the last column
   Ucollohihi[k] = Ulohihi[colidx+k];
   Ucolhilohi[k] = Uhilohi[colidx+k];
   Ucollolohi[k] = Ulolohi[colidx+k];
   Ucolhihilo[k] = Uhihilo[colidx+k];
   Ucollohilo[k] = Ulohilo[colidx+k];
   Ucolhilolo[k] = Uhilolo[colidx+k];
   Ucollololo[k] = Ulololo[colidx+k];
   rhshihihi = ((double) int(k == dim-1)); // right hand side for each thread
   rhslohihi = 0.0;
   rhshilohi = 0.0;
   rhslolohi = 0.0;
   rhshihilo = 0.0;
   rhslohilo = 0.0;
   rhshilolo = 0.0;
   rhslololo = 0.0;
   int rowidx = (dim - 1)*dim + k;       // the row index in the inverse

   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   odg_div( rhshihihi,        rhslohihi,       rhshilohi,       rhslolohi,
            rhshihilo,        rhslohilo,       rhshilolo,       rhslololo,
           Ucolhihihi[k],    Ucollohihi[k],   Ucolhilohi[k],   Ucollolohi[k],
           Ucolhihilo[k],    Ucollohilo[k],   Ucolhilolo[k],   Ucollololo[k],
       &invUrowhihihi[k],&invUrowlohihi[k],
       &invUrowhilohi[k],&invUrowlolohi[k],
       &invUrowhihilo[k],&invUrowlohilo[k],
       &invUrowhilolo[k],&invUrowlololo[k]);
   invUhihihi[rowidx] = invUrowhihihi[k];    // store the last row into invU
   invUlohihi[rowidx] = invUrowlohihi[k]; 
   invUhilohi[rowidx] = invUrowhilohi[k];
   invUlolohi[rowidx] = invUrowlolohi[k]; 
   invUhihilo[rowidx] = invUrowhihilo[k];
   invUlohilo[rowidx] = invUrowlohilo[k]; 
   invUhilolo[rowidx] = invUrowhilolo[k];
   invUlololo[rowidx] = invUrowlololo[k]; 

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshihihi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslohihi = 0.0;
      rhshilohi = 0.0;
      rhslolohi = 0.0;
      rhshihilo = 0.0;
      rhslohilo = 0.0;
      rhshilolo = 0.0;
      rhslololo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucolhihihi[k] = Uhihihi[colidx+k];
         Ucollohihi[k] = Ulohihi[colidx+k];
         Ucolhilohi[k] = Uhilohi[colidx+k];
         Ucollolohi[k] = Ulolohi[colidx+k];
         Ucolhihilo[k] = Uhihilo[colidx+k];
         Ucollohilo[k] = Ulohilo[colidx+k];
         Ucolhilolo[k] = Uhilolo[colidx+k];
         Ucollololo[k] = Ulololo[colidx+k];

         rowidx = j*dim + k;                // need solution value
         invUrowhihihi[k] = invUhihihi[rowidx]; // load invU row into invUrow
         invUrowlohihi[k] = invUlohihi[rowidx];
         invUrowhilohi[k] = invUhilohi[rowidx]; 
         invUrowlolohi[k] = invUlolohi[rowidx];
         invUrowhihilo[k] = invUhihilo[rowidx];
         invUrowlohilo[k] = invUlohilo[rowidx];
         invUrowhilolo[k] = invUhilolo[rowidx]; 
         invUrowlololo[k] = invUlololo[rowidx];
         xvalhihihi = invUrowhihihi[k];
         xvallohihi = invUrowlohihi[k];
         xvalhilohi = invUrowhilohi[k];
         xvallolohi = invUrowlolohi[k];
         xvalhihilo = invUrowhihilo[k];
         xvallohilo = invUrowlohilo[k];
         xvalhilolo = invUrowhilolo[k];
         xvallololo = invUrowlololo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         odg_mul(Ucolhihihi[i],Ucollohihi[i],Ucolhilohi[i],Ucollolohi[i],
                 Ucolhihilo[i],Ucollohilo[i],Ucolhilolo[i],Ucollololo[i],
                 xvalhihihi,   xvallohihi,   xvalhilohi,   xvallolohi,
                 xvalhihilo,   xvallohilo,   xvalhilolo,   xvallololo,
                 &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
                 &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
         odg_dec(&rhshihihi,&rhslohihi,&rhshilohi,&rhslolohi,
                 &rhshihilo,&rhslohilo,&rhshilolo,&rhslololo,
                  acchihihi, acclohihi, acchilohi, acclolohi,
                  acchihilo, acclohilo, acchilolo, acclololo);
      }
      colidx = dim*i;                 // need column i of U
      Ucolhihihi[k] = Uhihihi[colidx+k];
      Ucollohihi[k] = Ulohihi[colidx+k];
      Ucolhilohi[k] = Uhilohi[colidx+k];
      Ucollolohi[k] = Ulolohi[colidx+k];
      Ucolhihilo[k] = Uhihilo[colidx+k];
      Ucollohilo[k] = Ulohilo[colidx+k];
      Ucolhilolo[k] = Uhilolo[colidx+k];
      Ucollololo[k] = Ulololo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      odg_div(rhshihihi,rhslohihi,rhshilohi,rhslolohi,
              rhshihilo,rhslohilo,rhshilolo,rhslololo,
              Ucolhihihi[i],Ucollohihi[i],Ucolhilohi[i],Ucollolohi[i],
              Ucolhihilo[i],Ucollohilo[i],Ucolhilolo[i],Ucollololo[i],
              &invUrowhihihi[k],&invUrowlohihi[k],
              &invUrowhilohi[k],&invUrowlolohi[k],
              &invUrowhihilo[k],&invUrowlohilo[k],
              &invUrowhilolo[k],&invUrowlololo[k]);
      invUhihihi[rowidx] = invUrowhihihi[k];
      invUlohihi[rowidx] = invUrowlohihi[k];
      invUhilohi[rowidx] = invUrowhilohi[k];
      invUlolohi[rowidx] = invUrowlolohi[k];
      invUhihilo[rowidx] = invUrowhihilo[k];
      invUlohilo[rowidx] = invUrowlohilo[k];
      invUhilolo[rowidx] = invUrowhilolo[k];
      invUlololo[rowidx] = invUrowlololo[k];
   }
}

/*

__global__ void cmplx4_medium_invert_upper
 ( int dim, 
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo )
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolrehihi[tabsqd_shmemsize];    // one column of U
   __shared__ double Ucolrelohi[tabsqd_shmemsize]; 
   __shared__ double Ucolrehilo[tabsqd_shmemsize];
   __shared__ double Ucolrelolo[tabsqd_shmemsize]; 
   __shared__ double Ucolimhihi[tabsqd_shmemsize];
   __shared__ double Ucolimlohi[tabsqd_shmemsize]; 
   __shared__ double Ucolimhilo[tabsqd_shmemsize];
   __shared__ double Ucolimlolo[tabsqd_shmemsize]; 
   __shared__ double invUrowrehihi[tabsqd_shmemsize]; // one row of invU
   __shared__ double invUrowrelohi[tabsqd_shmemsize]; 
   __shared__ double invUrowrehilo[tabsqd_shmemsize];
   __shared__ double invUrowrelolo[tabsqd_shmemsize]; 
   __shared__ double invUrowimhihi[tabsqd_shmemsize]; 
   __shared__ double invUrowimlohi[tabsqd_shmemsize]; 
   __shared__ double invUrowimhilo[tabsqd_shmemsize]; 
   __shared__ double invUrowimlolo[tabsqd_shmemsize]; 

   double rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo;
   double rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo;
   double xvalrehihi,xvalrelohi,xvalrehilo,xvalrelolo;
   double xvalimhihi,xvalimlohi,xvalimhilo,xvalimlolo;
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;
   double invrehihi,invrelohi,invrehilo,invrelolo;
   double invimhihi,invimlohi,invimhilo,invimlolo;
   double denhihi,denlohi,denhilo,denlolo;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolrehihi[k] = Urehihi[colidx+k];      // load the last column
   Ucolrelohi[k] = Urelohi[colidx+k];
   Ucolrehilo[k] = Urehilo[colidx+k]; 
   Ucolrelolo[k] = Urelolo[colidx+k];
   Ucolimhihi[k] = Uimhihi[colidx+k];
   Ucolimlohi[k] = Uimlohi[colidx+k];
   Ucolimhilo[k] = Uimhilo[colidx+k];
   Ucolimlolo[k] = Uimlolo[colidx+k];
   rhsrehihi = ((double) int(k == dim-1)); // right hand side for each thread
   rhsrelohi = 0.0;
   rhsrehilo = 0.0;
   rhsrelolo = 0.0;
   rhsimhihi = 0.0;
   rhsimlohi = 0.0;
   rhsimhilo = 0.0;
   rhsimlolo = 0.0;
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   __syncthreads();
   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   qdg_mul(Ucolrehihi[k],Ucolrelohi[k],Ucolrehilo[k],Ucolrelolo[k],
           Ucolrehihi[k],Ucolrelohi[k],Ucolrehilo[k],Ucolrelolo[k],
             &denhihi,     &denlohi,     &denhilo,     &denlolo);
   qdg_mul(Ucolimhihi[k],Ucolimlohi[k],Ucolimhilo[k],Ucolimlolo[k],
           Ucolimhihi[k],Ucolimlohi[k],Ucolimhilo[k],Ucolimlolo[k],
            &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
   qdg_inc(&denhihi,&denlohi,&denhilo,&denlolo,
           acc1hihi,acc1lohi,acc1hilo,acc1lolo);
   qdg_div(Ucolrehihi[k],Ucolrelohi[k],Ucolrehilo[k],Ucolrelolo[k],
              denhihi,      denlohi,      denhilo,      denlolo,
           &invrehihi,   &invrelohi,   &invrehilo,   &invrelolo);
   qdg_div(Ucolimhihi[k],Ucolimlohi[k],Ucolimhilo[k],Ucolimlolo[k],
              denhihi,      denlohi,      denhilo,      denlolo,
           &invimhihi,   &invimlohi,   &invimhilo,   &invimlolo);
   qdg_minus(&invimhihi,&invimlohi,
             &invimhilo,&invimlolo);
   qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
           invrehihi,invrelohi,invrehilo,invrelolo,
           &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
   qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
           invimhihi,invimlohi,invimhilo,invimlolo,
           &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
   qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
           invrehihi,invrelohi,invrehilo,invrelolo,
           &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
   qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
           invimhihi,invimlohi,invimhilo,invimlolo,
           &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
   qdg_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
            acc2hihi, acc2lohi, acc2hilo, acc2lolo);
   invUrowrehihi[k] = acc1hihi;
   invUrowrelohi[k] = acc1lohi;
   invUrowrehilo[k] = acc1hilo;
   invUrowrelolo[k] = acc1lolo;
   qdg_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
            acc4hihi, acc4lohi, acc4hilo, acc4lolo);
   invUrowimhihi[k] = acc3hihi;
   invUrowimlohi[k] = acc3lohi;
   invUrowimhilo[k] = acc3hilo;
   invUrowimlolo[k] = acc3lolo;
   invUrehihi[rowidx] = invUrowrehihi[k];  // store the last row into invU
   invUrelohi[rowidx] = invUrowrelohi[k]; 
   invUrehilo[rowidx] = invUrowrehilo[k]; 
   invUrelolo[rowidx] = invUrowrelolo[k]; 
   invUimhihi[rowidx] = invUrowimhihi[k];
   invUimlohi[rowidx] = invUrowimlohi[k]; 
   invUimhilo[rowidx] = invUrowimhilo[k];
   invUimlolo[rowidx] = invUrowimlolo[k]; 

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehihi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhsrelohi = 0.0;
      rhsrehilo = 0.0;
      rhsrelolo = 0.0;
      rhsimhihi = 0.0;
      rhsimlohi = 0.0;
      rhsimhilo = 0.0;
      rhsimlolo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucolrehihi[k] = Urehihi[colidx+k];
         Ucolrelohi[k] = Urelohi[colidx+k];
         Ucolrehilo[k] = Urehilo[colidx+k];
         Ucolrelolo[k] = Urelolo[colidx+k];
         Ucolimhihi[k] = Uimhihi[colidx+k];
         Ucolimlohi[k] = Uimlohi[colidx+k];
         Ucolimhilo[k] = Uimhilo[colidx+k];
         Ucolimlolo[k] = Uimlolo[colidx+k];

         rowidx = j*dim + k;                // need solution value
         invUrowrehihi[k] = invUrehihi[rowidx]; // load invU row into invUrow
         invUrowrelohi[k] = invUrelohi[rowidx];
         invUrowrehilo[k] = invUrehilo[rowidx]; 
         invUrowrelolo[k] = invUrelolo[rowidx];
         invUrowimhihi[k] = invUimhihi[rowidx];
         invUrowimlohi[k] = invUimlohi[rowidx];
         invUrowimhilo[k] = invUimhilo[rowidx];
         invUrowimlolo[k] = invUimlolo[rowidx];
         xvalrehihi = invUrowrehihi[k];
         xvalrelohi = invUrowrelohi[k];
         xvalrehilo = invUrowrehilo[k];
         xvalrelolo = invUrowrelolo[k];
         xvalimhihi = invUrowimhihi[k];
         xvalimlohi = invUrowimlohi[k];
         xvalimhilo = invUrowimhilo[k];
         xvalimlolo = invUrowimlolo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         qdg_mul(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                 xvalrehihi,   xvalrelohi,   xvalrehilo,   xvalrelolo,
                  &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdg_mul(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
                 xvalimhihi,   xvalimlohi,   xvalimhilo,   xvalimlolo,
                  &acc2hihi,    &acc2lohi,    &acc2hilo,    &acc2lolo);
         qdg_mul(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
                 xvalrehihi,   xvalrelohi,   xvalrehilo,   xvalrelolo,
                  &acc3hihi,    &acc3lohi,    &acc3hilo,    &acc3lolo);
         qdg_mul(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                 xvalimhihi,   xvalimlohi,   xvalimhilo,   xvalimlolo,
                  &acc4hihi,    &acc4lohi,    &acc4hilo,    &acc4lolo);
         qdg_dec(&rhsrehihi,&rhsrelohi,&rhsrehilo,&rhsrelolo,
                   acc1hihi,  acc1lohi,  acc1hilo,  acc1lolo);
         qdg_inc(&rhsrehihi,&rhsrelohi,&rhsrehilo,&rhsrelolo,
                   acc2hihi,  acc2lohi,  acc2hilo,  acc2lolo);
         qdg_dec(&rhsimhihi,&rhsimlohi,&rhsimhilo,&rhsimlolo,
                   acc3hihi,  acc3lohi,  acc3hilo,  acc3lolo);
         qdg_dec(&rhsimhihi,&rhsimlohi,&rhsimhilo,&rhsimlolo,
                   acc4hihi,  acc4lohi,  acc4hilo,  acc4lolo);
      }
      colidx = dim*i;                 // need column i of U
      Ucolrehihi[k] = Urehihi[colidx+k];
      Ucolrelohi[k] = Urelohi[colidx+k];
      Ucolrehilo[k] = Urehilo[colidx+k];
      Ucolrelolo[k] = Urelolo[colidx+k];
      Ucolimhihi[k] = Uimhihi[colidx+k];
      Ucolimlohi[k] = Uimlohi[colidx+k];
      Ucolimhilo[k] = Uimhilo[colidx+k];
      Ucolimlolo[k] = Uimlolo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      qdg_mul(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
              Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                &denhihi,     &denlohi,     &denhilo,     &denlolo);
      qdg_mul(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
              Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
      qdg_inc(&denhihi,&denlohi,&denhilo,&denlolo,
              acc1hihi,acc1lohi,acc1hilo,acc1lolo);
      qdg_div(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                 denhihi,      denlohi,      denhilo,      denlolo,
              &invrehihi,   &invrelohi,   &invrehilo,   &invrelolo);
      qdg_div(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
                 denhihi,      denlohi,      denhilo,      denlolo,
              &invimhihi,   &invimlohi,   &invimhilo,   &invimlolo);
      qdg_minus(&invimhihi,&invimlohi,&invimhilo,&invimlolo);
      qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
              invrehihi,invrelohi,invrehilo,invrelolo,
              &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
      qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
              invimhihi,invimlohi,invimhilo,invimlolo,
              &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
      qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
              invrehihi,invrelohi,invrehilo,invrelolo,
              &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
      qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
              invimhihi,invimlohi,invimhilo,invimlolo,
              &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
      qdg_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
               acc2hihi, acc2lohi, acc2hilo,acc2lolo);
      invUrowrehihi[k] = acc1hihi;
      invUrowrelohi[k] = acc1lohi;
      invUrowrehilo[k] = acc1hilo;
      invUrowrelolo[k] = acc1lolo;
      qdg_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
               acc4hihi, acc4lohi, acc4hilo, acc4lolo);
      invUrowimhihi[k] = acc3hihi;
      invUrowimlohi[k] = acc3lohi;
      invUrowimhilo[k] = acc3hilo;
      invUrowimlolo[k] = acc3lolo;
      invUrehihi[rowidx] = invUrowrehihi[k];
      invUrelohi[rowidx] = invUrowrelohi[k];
      invUrehilo[rowidx] = invUrowrehilo[k];
      invUrelolo[rowidx] = invUrowrelolo[k];
      invUimhihi[rowidx] = invUrowimhihi[k];
      invUimlohi[rowidx] = invUrowimlohi[k];
      invUimhilo[rowidx] = invUrowimhilo[k];
      invUimlolo[rowidx] = invUrowimlolo[k];
   }
}

__global__ void  dbl4_invert_tiles
 ( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo)
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolhihi[tabsqd_shmemsize];    // one column of U
   __shared__ double Ucollohi[tabsqd_shmemsize];
   __shared__ double Ucolhilo[tabsqd_shmemsize];
   __shared__ double Ucollolo[tabsqd_shmemsize];
   __shared__ double invUrowhihi[tabsqd_shmemsize]; // one row of invU
   __shared__ double invUrowlohi[tabsqd_shmemsize]; 
   __shared__ double invUrowhilo[tabsqd_shmemsize];
   __shared__ double invUrowlolo[tabsqd_shmemsize]; 

   double rhshihi,rhslohi,rhshilo,rhslolo;
   double xvalhihi,xvallohi,xvalhilo,xvallolo;
   double acchihi,acclohi,acchilo,acclolo;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolhihi[k] = Uhihi[colidx+k];     // load the last column
   Ucollohi[k] = Ulohi[colidx+k];
   Ucolhilo[k] = Uhilo[colidx+k]; 
   Ucollolo[k] = Ulolo[colidx+k];
   rhshihi = ((double) int(k == dim-1));  // right hand side for each thread
   rhslohi = 0.0;
   rhshilo = 0.0;
   rhslolo = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
   qdg_div(     rhshihi,        rhslohi,        rhshilo,        rhslolo,
               Ucolhihi[k],    Ucollohi[k],    Ucolhilo[k],    Ucollolo[k],
           &invUrowhihi[k],&invUrowlohi[k],&invUrowhilo[k],&invUrowlolo[k]);
   invUhihi[rowidx] = invUrowhihi[k];     // store the last row into invU
   invUlohi[rowidx] = invUrowlohi[k];
   invUhilo[rowidx] = invUrowhilo[k];
   invUlolo[rowidx] = invUrowlolo[k];

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshihi = ((double) int(k == i));   // set rhs for i-th unit vector
      rhslohi = 0.0;
      rhshilo = 0.0;
      rhslolo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;     // need column j of U
         Ucolhihi[k] = Uhihi[colidx+k];
         Ucollohi[k] = Ulohi[colidx+k];
         Ucolhilo[k] = Uhilo[colidx+k];
         Ucollolo[k] = Ulolo[colidx+k];

         rowidx = offset + j*dim + k;       // need solution value
         invUrowhihi[k] = invUhihi[rowidx]; // load invU row into invUrow
         invUrowlohi[k] = invUlohi[rowidx];
         invUrowhilo[k] = invUhilo[rowidx];
         invUrowlolo[k] = invUlolo[rowidx];
         xvalhihi = invUrowhihi[k];
         xvallohi = invUrowlohi[k];
         xvalhilo = invUrowhilo[k];
         xvallolo = invUrowlolo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         qdg_mul(Ucolhihi[i],Ucollohi[i],Ucolhilo[i],Ucollolo[i],
                 xvalhihi,   xvallohi,   xvalhilo,   xvallolo,
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdg_dec(&rhshihi,&rhslohi,&rhshilo,&rhslolo,
                  acchihi, acclohi, acchilo, acclolo);
      }
      colidx = offset + dim*i;        // need column i of U
      Ucolhihi[k] = Uhihi[colidx+k];
      Ucollohi[k] = Ulohi[colidx+k];
      Ucolhilo[k] = Uhilo[colidx+k];
      Ucollolo[k] = Ulolo[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      qdg_div(     rhshihi,        rhslohi,        rhshilo,        rhslolo,
                  Ucolhihi[i],    Ucollohi[i],    Ucolhilo[i],    Ucollolo[i],
              &invUrowhihi[k],&invUrowlohi[k],&invUrowhilo[k],&invUrowlolo[k]);
      invUhihi[rowidx] = invUrowhihi[k];
      invUlohi[rowidx] = invUrowlohi[k];
      invUhilo[rowidx] = invUrowhilo[k];
      invUlolo[rowidx] = invUrowlolo[k];
   }
}

__global__ void  cmplx4_invert_tiles
 ( int dim, 
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolrehihi[tabsqd_shmemsize];    // one column of U
   __shared__ double Ucolrelohi[tabsqd_shmemsize];
   __shared__ double Ucolrehilo[tabsqd_shmemsize];
   __shared__ double Ucolrelolo[tabsqd_shmemsize];
   __shared__ double Ucolimhihi[tabsqd_shmemsize]; 
   __shared__ double Ucolimlohi[tabsqd_shmemsize];
   __shared__ double Ucolimhilo[tabsqd_shmemsize]; 
   __shared__ double Ucolimlolo[tabsqd_shmemsize];
   __shared__ double invUrowrehihi[tabsqd_shmemsize];  // one row of invU
   __shared__ double invUrowrelohi[tabsqd_shmemsize]; 
   __shared__ double invUrowrehilo[tabsqd_shmemsize]; 
   __shared__ double invUrowrelolo[tabsqd_shmemsize]; 
   __shared__ double invUrowimhihi[tabsqd_shmemsize];
   __shared__ double invUrowimlohi[tabsqd_shmemsize]; 
   __shared__ double invUrowimhilo[tabsqd_shmemsize];
   __shared__ double invUrowimlolo[tabsqd_shmemsize]; 

   double rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo;
   double rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo;
   double xvalrehihi,xvalrelohi,xvalrehilo,xvalrelolo;
   double xvalimhihi,xvalimlohi,xvalimhilo,xvalimlolo;
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;
   double invrehihi,invrelohi,invrehilo,invrelolo;
   double invimhihi,invimlohi,invimhilo,invimlolo;
   double denhihi,denlohi,denhilo,denlolo;

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolrehihi[k] = Urehihi[colidx+k];   // load the last column
   Ucolrelohi[k] = Urelohi[colidx+k];
   Ucolrehilo[k] = Urehilo[colidx+k];
   Ucolrelolo[k] = Urelolo[colidx+k];
   Ucolimhihi[k] = Uimhihi[colidx+k];
   Ucolimlohi[k] = Uimlohi[colidx+k];
   Ucolimhilo[k] = Uimhilo[colidx+k];
   Ucolimlolo[k] = Uimlolo[colidx+k];
   rhsrehihi = ((double) int(k == dim-1)); // right hand side for each thread
   rhsrelohi = 0.0;
   rhsrehilo = 0.0;
   rhsrelolo = 0.0;
   rhsimhihi = 0.0;
   rhsimlohi = 0.0;
   rhsimhilo = 0.0;
   rhsimlolo = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
   qdg_mul(Ucolrehihi[k],Ucolrelohi[k],Ucolrehilo[k],Ucolrelolo[k],
           Ucolrehihi[k],Ucolrelohi[k],Ucolrehilo[k],Ucolrelolo[k],
             &denhihi,     &denlohi,     &denhilo,     &denlolo);
   qdg_mul(Ucolimhihi[k],Ucolimlohi[k],Ucolimhilo[k],Ucolimlolo[k],
           Ucolimhihi[k],Ucolimlohi[k],Ucolimhilo[k],Ucolimlolo[k],
            &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
   qdg_inc(&denhihi,&denlohi,&denhilo,&denlolo,
           acc1hihi,acc1lohi,acc1hilo,acc1lolo);
   qdg_div(Ucolrehihi[k],Ucolrelohi[k],Ucolrehilo[k],Ucolrelolo[k],
              denhihi,      denlohi,      denhilo,      denlolo,
           &invrehihi,   &invrelohi,   &invrehilo,   &invrelolo);
   qdg_div(Ucolimhihi[k],Ucolimlohi[k],Ucolimhilo[k],Ucolimlolo[k],
              denhihi,      denlohi,      denhilo,      denlolo,
           &invimhihi,   &invimlohi,   &invimhilo,   &invimlolo);
   qdg_minus(&invimhihi,&invimlohi,&invimhilo,&invimlolo);
   qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
           invrehihi,invrelohi,invrehilo,invrelolo,
           &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
   qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
           invimhihi,invimlohi,invimhilo,invimlolo,
           &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
   qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
           invrehihi,invrelohi,invrehilo,invrelolo,
           &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
   qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
           invimhihi,invimlohi,invimhilo,invimlolo,
           &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
   qdg_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
            acc2hihi, acc2lohi, acc2hilo, acc2lolo);
   invUrowrehihi[k] = acc1hihi;
   invUrowrelohi[k] = acc1lohi;
   invUrowrehilo[k] = acc1hilo;
   invUrowrelolo[k] = acc1lolo;
   qdg_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
            acc4hihi, acc4lohi, acc4hilo, acc4lolo);
   invUrowimhihi[k] = acc3hihi;
   invUrowimlohi[k] = acc3lohi;
   invUrowimhilo[k] = acc3hilo;
   invUrowimlolo[k] = acc3lolo;
   invUrehihi[rowidx] = invUrowrehihi[k];   // store the last row into invU
   invUrelohi[rowidx] = invUrowrelohi[k];
   invUrehilo[rowidx] = invUrowrehilo[k]; 
   invUrelolo[rowidx] = invUrowrelolo[k];
   invUimhihi[rowidx] = invUrowimhihi[k];
   invUimlohi[rowidx] = invUrowimlohi[k];
   invUimhilo[rowidx] = invUrowimhilo[k];
   invUimlolo[rowidx] = invUrowimlolo[k];

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehihi = ((double) int(k == i));   // set rhs for i-th unit vector
      rhsrelohi = 0.0;
      rhsrehilo = 0.0;
      rhsrelolo = 0.0;
      rhsimhihi = 0.0;
      rhsimlohi = 0.0;
      rhsimhilo = 0.0;
      rhsimlolo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;        // need column j of U
         Ucolrehihi[k] = Urehihi[colidx+k];
         Ucolrelohi[k] = Urelohi[colidx+k];
         Ucolrehilo[k] = Urehilo[colidx+k];
         Ucolrelolo[k] = Urelolo[colidx+k];
         Ucolimhihi[k] = Uimhihi[colidx+k];
         Ucolimlohi[k] = Uimlohi[colidx+k];
         Ucolimhilo[k] = Uimhilo[colidx+k];
         Ucolimlolo[k] = Uimlolo[colidx+k];

         rowidx = offset + j*dim + k;       // need solution value
         invUrowrehihi[k] = invUrehihi[rowidx]; // load invU row into invUrow
         invUrowrelohi[k] = invUrelohi[rowidx];
         invUrowrehilo[k] = invUrehilo[rowidx];
         invUrowrelolo[k] = invUrelolo[rowidx];
         invUrowimhihi[k] = invUimhihi[rowidx];
         invUrowimlohi[k] = invUimlohi[rowidx];
         invUrowimhilo[k] = invUimhilo[rowidx];
         invUrowimlolo[k] = invUimlolo[rowidx];
         xvalrehihi = invUrowrehihi[k];
         xvalrelohi = invUrowrelohi[k];
         xvalrehilo = invUrowrehilo[k];
         xvalrelolo = invUrowrelolo[k];
         xvalimhihi = invUrowimhihi[k];
         xvalimlohi = invUrowimlohi[k];
         xvalimhilo = invUrowimhilo[k];
         xvalimlolo = invUrowimlolo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         qdg_mul(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                 xvalrehihi,   xvalrelohi,   xvalrehilo,   xvalrelolo,
                  &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
         qdg_mul(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
                 xvalimhihi,   xvalimlohi,   xvalimhilo,   xvalimlolo,
                  &acc2hihi,    &acc2lohi,    &acc2hilo,    &acc2lolo);
         qdg_mul(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
                 xvalrehihi,   xvalrelohi,   xvalrehilo,   xvalrelolo,
                  &acc3hihi,    &acc3lohi,    &acc3hilo,    &acc3lolo);
         qdg_mul(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                 xvalimhihi,   xvalimlohi,   xvalimhilo,   xvalimlolo,
                  &acc4hihi,    &acc4lohi,    &acc4hilo,    &acc4lolo);
         qdg_dec(&rhsrehihi,&rhsrelohi,&rhsrehilo,&rhsrelolo,
                   acc1hihi,  acc1lohi,  acc1hilo,  acc1lolo);
         qdg_inc(&rhsrehihi,&rhsrelohi,&rhsrehilo,&rhsrelolo,
                   acc2hihi,  acc2lohi,  acc2hilo,  acc2lolo);
         qdg_dec(&rhsimhihi,&rhsimlohi,&rhsimhilo,&rhsimlolo,
                   acc3hihi,  acc3lohi,  acc3hilo,  acc3lolo);
         qdg_dec(&rhsimhihi,&rhsimlohi,&rhsimhilo,&rhsimlolo,
                   acc4hihi,  acc4lohi,  acc4hilo,  acc4lolo);
      }
      colidx = offset + dim*i;        // need column i of U
      Ucolrehihi[k] = Urehihi[colidx+k];
      Ucolrelohi[k] = Urelohi[colidx+k];
      Ucolrehilo[k] = Urehilo[colidx+k];
      Ucolrelolo[k] = Urelolo[colidx+k];
      Ucolimhihi[k] = Uimhihi[colidx+k];
      Ucolimlohi[k] = Uimlohi[colidx+k];
      Ucolimhilo[k] = Uimhilo[colidx+k];
      Ucolimlolo[k] = Uimlolo[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      qdg_mul(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
              Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                &denhihi,     &denlohi,     &denhilo,     &denlolo);
      qdg_mul(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
              Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
               &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
      qdg_inc(&denhihi,&denlohi,&denhilo,&denlolo,
              acc1hihi,acc1lohi,acc1hilo,acc1lolo);
      qdg_div(Ucolrehihi[i],Ucolrelohi[i],Ucolrehilo[i],Ucolrelolo[i],
                 denhihi,      denlohi,      denhilo,      denlolo,
              &invrehihi,   &invrelohi,   &invrehilo,   &invrelolo);
      qdg_div(Ucolimhihi[i],Ucolimlohi[i],Ucolimhilo[i],Ucolimlolo[i],
                 denhihi,      denlohi,      denhilo,      denlolo,
              &invimhihi,   &invimlohi,   &invimhilo,   &invimlolo);
      qdg_minus(&invimhihi,&invimlohi,&invimhilo,&invimlolo);
      qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
              invrehihi,invrelohi,invrehilo,invrelolo,
              &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
      qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
              invimhihi,invimlohi,invimhilo,invimlolo,
              &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
      qdg_mul(rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
              invrehihi,invrelohi,invrehilo,invrelolo,
              &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
      qdg_mul(rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
              invimhihi,invimlohi,invimhilo,invimlolo,
              &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
      qdg_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
               acc2hihi, acc2lohi, acc2hilo, acc2lolo);
      invUrowrehihi[k] = acc1hihi;
      invUrowrelohi[k] = acc1lohi;
      invUrowrehilo[k] = acc1hilo;
      invUrowrelolo[k] = acc1lolo;
      qdg_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
               acc4hihi, acc4lohi, acc4hilo, acc4lolo);
      invUrowimhihi[k] = acc3hihi;
      invUrowimlohi[k] = acc3lohi;
      invUrowimhilo[k] = acc3hilo;
      invUrowimlolo[k] = acc3lolo;
      invUrehihi[rowidx] = invUrowrehihi[k];
      invUrelohi[rowidx] = invUrowrelohi[k];
      invUrehilo[rowidx] = invUrowrehilo[k];
      invUrelolo[rowidx] = invUrowrelolo[k];
      invUimhihi[rowidx] = invUrowimhihi[k];
      invUimlohi[rowidx] = invUrowimlohi[k];
      invUimhilo[rowidx] = invUrowimhilo[k];
      invUimlolo[rowidx] = invUrowimlolo[k];
   }
}

__global__ void dbl4_multiply_inverse
 ( int dim, int idx,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo,
   double *whihi, double *wlohi, double *whilo, double *wlolo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workhihi[tabsqd_shmemsize];      // copy of w
   __shared__ double worklohi[tabsqd_shmemsize];
   __shared__ double workhilo[tabsqd_shmemsize];
   __shared__ double worklolo[tabsqd_shmemsize];

   workhihi[k] = whihi[rhsoff+k];
   worklohi[k] = wlohi[rhsoff+k];
   workhilo[k] = whilo[rhsoff+k];
   worklolo[k] = wlolo[rhsoff+k];

   double resulthihi = 0.0; // each thread stores its product in result
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double coeffhihi,coefflohi,coeffhilo,coefflolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhihi = invUhihi[offset+k*dim+j]; // thread k does row k
      coefflohi = invUlohi[offset+k*dim+j];
      coeffhilo = invUhilo[offset+k*dim+j];
      coefflolo = invUlolo[offset+k*dim+j];
      // result = result + coeff*work[j];
      qdg_mul(coeffhihi,  coefflohi,  coeffhilo,  coefflolo,
               workhihi[j],worklohi[j],workhilo[j],worklolo[j],
               &acchihi,   &acclohi,   &acchilo,   &acclolo);
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                  acchihi,    acclohi,    acchilo,    acclolo);
   }
   whihi[rhsoff+k] = resulthihi;
   wlohi[rhsoff+k] = resultlohi;
   whilo[rhsoff+k] = resulthilo;
   wlolo[rhsoff+k] = resultlolo;
}

__global__ void cmplx4_multiply_inverse
 ( int dim, int idx,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workrehihi[tabsqd_shmemsize];      // copy of w
   __shared__ double workrelohi[tabsqd_shmemsize]; 
   __shared__ double workrehilo[tabsqd_shmemsize];
   __shared__ double workrelolo[tabsqd_shmemsize]; 
   __shared__ double workimhihi[tabsqd_shmemsize];
   __shared__ double workimlohi[tabsqd_shmemsize];
   __shared__ double workimhilo[tabsqd_shmemsize];
   __shared__ double workimlolo[tabsqd_shmemsize];

   workrehihi[k] = wrehihi[rhsoff+k];
   workrelohi[k] = wrelohi[rhsoff+k];
   workrehilo[k] = wrehilo[rhsoff+k];
   workrelolo[k] = wrelolo[rhsoff+k];
   workimhihi[k] = wimhihi[rhsoff+k];
   workimlohi[k] = wimlohi[rhsoff+k];
   workimhilo[k] = wimhilo[rhsoff+k];
   workimlolo[k] = wimlolo[rhsoff+k];

   double resultrehihi = 0.0; // each thread stores its product in result
   double resultrelohi = 0.0;
   double resultrehilo = 0.0;
   double resultrelolo = 0.0;
   double resultimhihi = 0.0;
   double resultimlohi = 0.0;
   double resultimhilo = 0.0;
   double resultimlolo = 0.0;
   double coeffrehihi,coeffrelohi,coeffrehilo,coeffrelolo;
   double coeffimhihi,coeffimlohi,coeffimhilo,coeffimlolo;
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffrehihi = invUrehihi[offset+k*dim+j]; // thread k does row k
      coeffrelohi = invUrelohi[offset+k*dim+j];
      coeffrehilo = invUrehilo[offset+k*dim+j];
      coeffrelolo = invUrelolo[offset+k*dim+j];
      coeffimhihi = invUimhihi[offset+k*dim+j];
      coeffimlohi = invUimlohi[offset+k*dim+j];
      coeffimhilo = invUimhilo[offset+k*dim+j];
      coeffimlolo = invUimlolo[offset+k*dim+j];
      // result = result + coeff*work[j];
      qdg_mul(coeffrehihi,  coeffrelohi,  coeffrehilo,  coeffrelolo,
               workrehihi[j],workrelohi[j],workrehilo[j],workrelolo[j],
                &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
      qdg_mul(coeffimhihi,  coeffimlohi,  coeffimhilo,  coeffimlolo,
               workimhihi[j],workimlohi[j],workimhilo[j],workimlolo[j],
                &acc2hihi,    &acc2lohi,    &acc2hilo,    &acc2lolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
      qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                   acc2hihi,     acc2lohi,     acc2hilo,     acc2lolo);
      qdg_mul(coeffimhihi,  coeffimlohi,  coeffimhilo,  coeffimlolo,
               workrehihi[j],workrelohi[j],workrehilo[j],workrelolo[j],
                &acc1hihi,    &acc1lohi,    &acc1hilo,    &acc1lolo);
      qdg_mul(coeffrehihi,  coeffrelohi,  coeffrehilo,  coeffrelolo,
               workimhihi[j],workimlohi[j],workimhilo[j],workimlolo[j],
                &acc2hihi,    &acc2lohi,    &acc2hilo,    &acc2lolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                   acc2hihi,     acc2lohi,     acc2hilo,    acc2lolo);
   }
   wrehihi[rhsoff+k] = resultrehihi;
   wrelohi[rhsoff+k] = resultrelohi;
   wrehilo[rhsoff+k] = resultrehilo;
   wrelolo[rhsoff+k] = resultrelolo;
   wimhihi[rhsoff+k] = resultimhihi;
   wimlohi[rhsoff+k] = resultimlohi;
   wimhilo[rhsoff+k] = resultimhilo;
   wimlolo[rhsoff+k] = resultimlolo;
}

__global__ void dbl4_back_substitute
 ( int dim, int idx, 
   double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo, 
   double *whihi, double *wlohi, double *whilo, double *wlolo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkhihi[tabsqd_shmemsize];  // copy of w
   __shared__ double wrklohi[tabsqd_shmemsize]; 
   __shared__ double wrkhilo[tabsqd_shmemsize];
   __shared__ double wrklolo[tabsqd_shmemsize]; 
   __shared__ double solhihi[tabsqd_shmemsize];  // solution to update with
   __shared__ double sollohi[tabsqd_shmemsize];
   __shared__ double solhilo[tabsqd_shmemsize];
   __shared__ double sollolo[tabsqd_shmemsize];

   wrkhihi[k] = whihi[B*dim+k];    // block B updates B-th slice of w
   wrklohi[k] = wlohi[B*dim+k];
   wrkhilo[k] = whilo[B*dim+k];
   wrklolo[k] = wlolo[B*dim+k];
   solhihi[k] = whihi[idx*dim+k];  // solution that is back substituted
   sollohi[k] = wlohi[idx*dim+k];
   solhilo[k] = whilo[idx*dim+k];
   sollolo[k] = wlolo[idx*dim+k];

   double resulthihi = 0.0; // each thread stores its product in result
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;
   double coeffhihi,coefflohi,coeffhilo,coefflolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhihi = Uhihi[offset+k*dim+j];
      coefflohi = Ulohi[offset+k*dim+j];
      coeffhilo = Uhilo[offset+k*dim+j];
      coefflolo = Ulolo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      qdg_mul(coeffhihi, coefflohi, coeffhilo, coefflolo,
                solhihi[j],sollohi[j],solhilo[j],sollolo[j],
               &acchihi,  &acclohi,  &acchilo,  &acclolo);
      qdg_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
                  acchihi,    acclohi,    acchilo,    acclolo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   qdg_dec(  &wrkhihi[k],&wrklohi[k],&wrkhilo[k],&wrklolo[k],
           resulthihi, resultlohi, resulthilo, resultlolo);
   whihi[B*dim+k] = wrkhihi[k];
   wlohi[B*dim+k] = wrklohi[k];
   whilo[B*dim+k] = wrkhilo[k];
   wlolo[B*dim+k] = wrklolo[k];
}

__global__ void cmplx4_back_substitute
 ( int dim, int idx,
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *wrehihi, double *wrelohi, double *wrehilo, double *wrelolo,
   double *wimhihi, double *wimlohi, double *wimhilo, double *wimlolo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkrehihi[tabsqd_shmemsize]; // copy of w
   __shared__ double wrkrelohi[tabsqd_shmemsize]; 
   __shared__ double wrkrehilo[tabsqd_shmemsize];
   __shared__ double wrkrelolo[tabsqd_shmemsize]; 
   __shared__ double wrkimhihi[tabsqd_shmemsize];
   __shared__ double wrkimlohi[tabsqd_shmemsize]; 
   __shared__ double wrkimhilo[tabsqd_shmemsize];
   __shared__ double wrkimlolo[tabsqd_shmemsize]; 
   __shared__ double solrehihi[tabsqd_shmemsize]; // solution to update with
   __shared__ double solrelohi[tabsqd_shmemsize];
   __shared__ double solrehilo[tabsqd_shmemsize];
   __shared__ double solrelolo[tabsqd_shmemsize];
   __shared__ double solimhihi[tabsqd_shmemsize];
   __shared__ double solimlohi[tabsqd_shmemsize];
   __shared__ double solimhilo[tabsqd_shmemsize];
   __shared__ double solimlolo[tabsqd_shmemsize];

   wrkrehihi[k] = wrehihi[B*dim+k];    // block B updates B-th slice of w
   wrkrelohi[k] = wrelohi[B*dim+k];
   wrkrehilo[k] = wrehilo[B*dim+k];
   wrkrelolo[k] = wrelolo[B*dim+k];
   wrkimhihi[k] = wimhihi[B*dim+k];
   wrkimlohi[k] = wimlohi[B*dim+k];
   wrkimhilo[k] = wimhilo[B*dim+k];
   wrkimlolo[k] = wimlolo[B*dim+k];
   solrehihi[k] = wrehihi[idx*dim+k];  // solution that is back substituted
   solrelohi[k] = wrelohi[idx*dim+k];
   solrehilo[k] = wrehilo[idx*dim+k];
   solrelolo[k] = wrelolo[idx*dim+k];
   solimhihi[k] = wimhihi[idx*dim+k];
   solimlohi[k] = wimlohi[idx*dim+k];
   solimhilo[k] = wimhilo[idx*dim+k];
   solimlolo[k] = wimlolo[idx*dim+k];

   double resultrehihi = 0.0; // each thread stores its product in result
   double resultrelohi = 0.0;
   double resultrehilo = 0.0;
   double resultrelolo = 0.0;
   double resultimhihi = 0.0;
   double resultimlohi = 0.0;
   double resultimhilo = 0.0;
   double resultimlolo = 0.0;
   double coeffrehihi,coeffrelohi,coeffrehilo,coeffrelolo;
   double coeffimhihi,coeffimlohi,coeffimhilo,coeffimlolo;
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffrehihi = Urehihi[offset+k*dim+j];
      coeffrelohi = Urelohi[offset+k*dim+j];
      coeffrehilo = Urehilo[offset+k*dim+j];
      coeffrelolo = Urelolo[offset+k*dim+j];
      coeffimhihi = Uimhihi[offset+k*dim+j];
      coeffimlohi = Uimlohi[offset+k*dim+j];
      coeffimhilo = Uimhilo[offset+k*dim+j];
      coeffimlolo = Uimlolo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      qdg_mul(coeffrehihi, coeffrelohi, coeffrehilo, coeffrelolo,
                solrehihi[j],solrelohi[j],solrehilo[j],solrelolo[j],
                &acc1hihi,   &acc1lohi,   &acc1hilo,   &acc1lolo);
      qdg_mul(coeffimhihi, coeffimlohi, coeffimhilo, coeffimlolo,
                solimhihi[j],solimlohi[j],solimhilo[j],solimlolo[j],
                &acc2hihi,   &acc2lohi,   &acc2hilo,   &acc2lolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
      qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                   acc2hihi,     acc2lohi,     acc2hilo,     acc2lolo);
      qdg_mul(coeffimhihi, coeffimlohi, coeffimhilo, coeffimlolo,
                solrehihi[j],solrelohi[j],solrehilo[j],solrelolo[j],
                &acc1hihi,   &acc1lohi,   &acc1hilo,   &acc1lolo);
      qdg_mul(coeffrehihi, coeffrelohi, coeffrehilo, coeffrelolo,
                solimhihi[j],solimlohi[j],solimhilo[j],solimlolo[j],
                &acc2hihi,   &acc2lohi,   &acc2hilo,   &acc2lolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                   acc1hihi,     acc1lohi,     acc1hilo,     acc1lolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                   acc2hihi,     acc2lohi,     acc2hilo,     acc2lolo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   qdg_dec(  &wrkrehihi[k],&wrkrelohi[k],&wrkrehilo[k],&wrkrelolo[k],
           resultrehihi, resultrelohi, resultrehilo, resultrelolo);
   qdg_dec(  &wrkimhihi[k],&wrkimlohi[k],&wrkimhilo[k],&wrkimlolo[k],
           resultimhihi, resultimlohi, resultimhilo, resultimlolo);
   wrehihi[B*dim+k] = wrkrehihi[k];
   wrelohi[B*dim+k] = wrkrelohi[k];
   wrehilo[B*dim+k] = wrkrehilo[k];
   wrelolo[B*dim+k] = wrkrelolo[k];
   wimhihi[B*dim+k] = wrkimhihi[k];
   wimlohi[B*dim+k] = wrkimlohi[k];
   wimhilo[B*dim+k] = wrkimhilo[k];
   wimlolo[B*dim+k] = wrkimlolo[k];
}

*/

void GPU_dbl8_upper_inverse
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double **invUhihihi, double **invUlohihi,
   double **invUhilohi, double **invUlolohi,
   double **invUhihilo, double **invUlohilo,
   double **invUhilolo, double **invUlololo,
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Uhihihi_h = new double[szU];   // the columns of U
   double *Ulohihi_h = new double[szU];    
   double *Uhilohi_h = new double[szU];   
   double *Ulolohi_h = new double[szU]; 
   double *Uhihilo_h = new double[szU];
   double *Ulohilo_h = new double[szU];    
   double *Uhilolo_h = new double[szU];   
   double *Ulololo_h = new double[szU]; 
   double *Uhihihi_d;                     // the columns of U on the device
   double *Ulohihi_d;
   double *Uhilohi_d;
   double *Ulolohi_d;
   double *Uhihilo_d; 
   double *Ulohilo_d;
   double *Uhilolo_d;
   double *Ulololo_d;
   double *invUhihihi_h = new double[szU];  // columns of the inverse of U
   double *invUlohihi_h = new double[szU]; 
   double *invUhilohi_h = new double[szU]; 
   double *invUlolohi_h = new double[szU]; 
   double *invUhihilo_h = new double[szU]; 
   double *invUlohilo_h = new double[szU]; 
   double *invUhilolo_h = new double[szU]; 
   double *invUlololo_h = new double[szU]; 
   double *invUhihihi_d;                    // inverse of U on the device
   double *invUlohihi_d;
   double *invUhilohi_d;
   double *invUlolohi_d;
   double *invUhihilo_d;
   double *invUlohilo_d;
   double *invUhilolo_d;
   double *invUlololo_d;

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Uhihihi_h[ix]   = Uhihihi[i][j];
         Ulohihi_h[ix]   = Ulohihi[i][j];
         Uhilohi_h[ix]   = Uhilohi[i][j];
         Ulolohi_h[ix]   = Ulolohi[i][j];
         Uhihilo_h[ix]   = Uhihilo[i][j];
         Ulohilo_h[ix]   = Ulohilo[i][j];
         Uhilolo_h[ix]   = Uhilolo[i][j];
         Ulololo_h[ix++] = Ulololo[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Uhihihi_d,szmat);
   cudaMalloc((void**)&Ulohihi_d,szmat);
   cudaMalloc((void**)&Uhilohi_d,szmat);
   cudaMalloc((void**)&Ulolohi_d,szmat);
   cudaMalloc((void**)&Uhihilo_d,szmat);
   cudaMalloc((void**)&Ulohilo_d,szmat);
   cudaMalloc((void**)&Uhilolo_d,szmat);
   cudaMalloc((void**)&Ulololo_d,szmat);
   cudaMalloc((void**)&invUhihihi_d,szmat);
   cudaMalloc((void**)&invUlohihi_d,szmat);
   cudaMalloc((void**)&invUhilohi_d,szmat);
   cudaMalloc((void**)&invUlolohi_d,szmat);
   cudaMalloc((void**)&invUhihilo_d,szmat);
   cudaMalloc((void**)&invUlohilo_d,szmat);
   cudaMalloc((void**)&invUhilolo_d,szmat);
   cudaMalloc((void**)&invUlololo_d,szmat);
   cudaMemcpy(Uhihihi_d,Uhihihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulohihi_d,Ulohihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uhilohi_d,Uhilohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulolohi_d,Ulolohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uhihilo_d,Uhihilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulohilo_d,Ulohilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uhilolo_d,Uhilolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulololo_d,Ulololo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      dbl8_small_invert_upper<<<1,dim>>>
         (dim,   Uhihihi_d,   Ulohihi_d,   Uhilohi_d,   Ulolohi_d,
                 Uhihilo_d,   Ulohilo_d,   Uhilolo_d,   Ulololo_d,
              invUhihihi_d,invUlohihi_d,invUhilohi_d,invUlolohi_d,
              invUhihilo_d,invUlohilo_d,invUhilolo_d,invUlololo_d);
   else
      dbl8_medium_invert_upper<<<1,dim>>>
         (dim,   Uhihihi_d,   Ulohihi_d,   Uhilohi_d,   Ulolohi_d,
                 Uhihilo_d,   Ulohilo_d,   Uhilolo_d,   Ulololo_d,
              invUhihihi_d,invUlohihi_d,invUhilohi_d,invUlolohi_d,
              invUhihilo_d,invUlohilo_d,invUhilolo_d,invUlololo_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUhihihi_h,invUhihihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlohihi_h,invUlohihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUhilohi_h,invUhilohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlolohi_h,invUlolohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUhihilo_h,invUhihilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlohilo_h,invUlohilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUhilolo_h,invUhilolo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlololo_h,invUlololo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUhihihi[i][j] = invUhihihi_h[ix];
         invUlohihi[i][j] = invUlohihi_h[ix];
         invUhilohi[i][j] = invUhilohi_h[ix];
         invUlolohi[i][j] = invUlolohi_h[ix];
         invUhihilo[i][j] = invUhihilo_h[ix];
         invUlohilo[i][j] = invUlohilo_h[ix];
         invUhilolo[i][j] = invUhilolo_h[ix];
         invUlololo[i][j] = invUlololo_h[ix++];
      }

   free(Uhihihi_h); free(Ulohihi_h); free(Uhilohi_h); free(Ulolohi_h);
   free(Uhihilo_h); free(Ulohilo_h); free(Uhilolo_h); free(Ulololo_h);
   free(invUhihihi_h); free(invUlohihi_h);
   free(invUhilohi_h); free(invUlolohi_h);
   free(invUhihilo_h); free(invUlohilo_h);
   free(invUhilolo_h); free(invUlololo_h);
}

/*

void GPU_cmplx4_upper_inverse
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double **invUrehihi, double **invUrelohi,
   double **invUrehilo, double **invUrelolo,
   double **invUimhihi, double **invUimlohi, 
   double **invUimhilo, double **invUimlolo, 
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Urehihi_h = new double[szU];  // real parts of U
   double *Urelohi_h = new double[szU]; 
   double *Urehilo_h = new double[szU]; 
   double *Urelolo_h = new double[szU]; 
   double *Uimhihi_h = new double[szU];  // imaginary parts of U
   double *Uimlohi_h = new double[szU]; 
   double *Uimhilo_h = new double[szU]; 
   double *Uimlolo_h = new double[szU]; 
   double *Urehihi_d;                    // real parts on the device
   double *Urelohi_d;
   double *Urehilo_d;
   double *Urelolo_d;
   double *Uimhihi_d;                    // imaginary parts on the device
   double *Uimlohi_d;
   double *Uimhilo_d;
   double *Uimlolo_d;
   double *invUrehihi_h = new double[szU]; // real parts of inverse
   double *invUrelohi_h = new double[szU];
   double *invUrehilo_h = new double[szU];
   double *invUrelolo_h = new double[szU];
   double *invUimhihi_h = new double[szU]; // imaginary parts of inverse
   double *invUimlohi_h = new double[szU];
   double *invUimhilo_h = new double[szU];
   double *invUimlolo_h = new double[szU];
   double *invUrehihi_d;           // invUrehihi_d ~ invUrehihi_h on device
   double *invUrelohi_d;           // invUrelohi_d ~ invUrelohi_h on device
   double *invUrehilo_d;           // invUrehilo_d ~ invUrehilo_h on device
   double *invUrelolo_d;           // invUrelolo_d ~ invUrelolo_h on device
   double *invUimhihi_d;           // invUimhihi_d ~ invUimhihi_h on device
   double *invUimlohi_d;           // invUimlohi_d ~ invUimlohi_h on device
   double *invUimhilo_d;           // invUimhilo_d ~ invUimhilo_h on device
   double *invUimlolo_d;           // invUimlolo_d ~ invUimlolo_h on device

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Urehihi_h[ix]   = Urehihi[i][j];
         Urelohi_h[ix]   = Urelohi[i][j];
         Urehilo_h[ix]   = Urehilo[i][j];
         Urelolo_h[ix]   = Urelolo[i][j];
         Uimhihi_h[ix]   = Uimhihi[i][j];
         Uimlohi_h[ix]   = Uimlohi[i][j];
         Uimhilo_h[ix]   = Uimhilo[i][j];
         Uimlolo_h[ix++] = Uimlolo[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Urehihi_d,szmat);
   cudaMalloc((void**)&Urelohi_d,szmat);
   cudaMalloc((void**)&Urehilo_d,szmat);
   cudaMalloc((void**)&Urelolo_d,szmat);
   cudaMalloc((void**)&Uimhihi_d,szmat);
   cudaMalloc((void**)&Uimlohi_d,szmat);
   cudaMalloc((void**)&Uimhilo_d,szmat);
   cudaMalloc((void**)&Uimlolo_d,szmat);
   cudaMalloc((void**)&invUrehihi_d,szmat);
   cudaMalloc((void**)&invUrelohi_d,szmat);
   cudaMalloc((void**)&invUrehilo_d,szmat);
   cudaMalloc((void**)&invUrelolo_d,szmat);
   cudaMalloc((void**)&invUimhihi_d,szmat);
   cudaMalloc((void**)&invUimlohi_d,szmat);
   cudaMalloc((void**)&invUimhilo_d,szmat);
   cudaMalloc((void**)&invUimlolo_d,szmat);
   cudaMemcpy(Urehihi_d,Urehihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelohi_d,Urelohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urehilo_d,Urehilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelolo_d,Urelolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhihi_d,Uimhihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlohi_d,Uimlohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhilo_d,Uimhilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlolo_d,Uimlolo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      cmplx4_small_invert_upper<<<1,dim>>>
         (dim,   Urehihi_d,   Urelohi_d,   Urehilo_d,   Urelolo_d,
                 Uimhihi_d,   Uimlohi_d,   Uimhilo_d,   Uimlolo_d,
              invUrehihi_d,invUrelohi_d,invUrehilo_d,invUrelolo_d,
              invUimhihi_d,invUimlohi_d,invUimhilo_d,invUimlolo_d);
   else
      cmplx4_medium_invert_upper<<<1,dim>>>
         (dim,   Urehihi_d,   Urelohi_d,   Urehilo_d,   Urelolo_d,
                 Uimhihi_d,   Uimlohi_d,   Uimhilo_d,   Uimlolo_d,
              invUrehihi_d,invUrelohi_d,invUrehilo_d,invUrelolo_d,
              invUimhihi_d,invUimlohi_d,invUimhilo_d,invUimlolo_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUrehihi_h,invUrehihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelohi_h,invUrelohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrehilo_h,invUrehilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelolo_h,invUrelolo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhihi_h,invUimhihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlohi_h,invUimlohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhilo_h,invUimhilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlolo_h,invUimlolo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUrehihi[i][j] = invUrehihi_h[ix];
         invUrelohi[i][j] = invUrelohi_h[ix];
         invUrehilo[i][j] = invUrehilo_h[ix];
         invUrelolo[i][j] = invUrelolo_h[ix];
         invUimhihi[i][j] = invUimhihi_h[ix];
         invUimlohi[i][j] = invUimlohi_h[ix];
         invUimhilo[i][j] = invUimhilo_h[ix];
         invUimlolo[i][j] = invUimlolo_h[ix++];
      }

   free(Urehihi_h); free(Urelohi_h); free(Urehilo_h); free(Urelolo_h);
   free(Uimhihi_h); free(Uimlohi_h); free(Uimhilo_h); free(Uimlolo_h);
   free(invUrehihi_h); free(invUrelohi_h);
   free(invUrehilo_h); free(invUrelolo_h);
   free(invUimhihi_h); free(invUimlohi_h);
   free(invUimhilo_h); free(invUimlolo_h);
}

void GPU_dbl4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *Dhihi_h = new double[nbr];  // the diagonal tiles on the host
   double *Dlohi_h = new double[nbr];  // second highest
   double *Dhilo_h = new double[nbr];  // second lowest diagonal tiles
   double *Dlolo_h = new double[nbr];  // lowest doubles of diagonal tiles
   double *Dhihi_d;                    // diagonal tiles on the device
   double *Dlohi_d;                    // second highest diagonal tiles
   double *Dhilo_d;                    // second lowest diagonal tiles
   double *Dlolo_d;                    // lowest doubles of diagonal tiles
   double *invDhihi_h = new double[nbr]; // inverse of diagonal tiles on host 
   double *invDlohi_h = new double[nbr]; 
   double *invDhilo_h = new double[nbr]; 
   double *invDlolo_h = new double[nbr];
   double *invDhihi_d;            // invDhihi_d is invDhihi_h on device
   double *invDlohi_d;            // invDlohi_d is invDlohi_h on device
   double *invDhilo_d;            // invDhilo_d is invDhilo_h on device
   double *invDlolo_d;            // invDlolo_d is invDlolo_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Dhihi_h[ix]   = Uhihi[offset+i][offset+j];
            Dlohi_h[ix]   = Ulohi[offset+i][offset+j];
            Dhilo_h[ix]   = Uhilo[offset+i][offset+j];
            Dlolo_h[ix++] = Ulolo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Dhihi_d,sznum);
   cudaMalloc((void**)&Dlohi_d,sznum);
   cudaMalloc((void**)&Dhilo_d,sznum);
   cudaMalloc((void**)&Dlolo_d,sznum);
   cudaMalloc((void**)&invDhihi_d,sznum);
   cudaMalloc((void**)&invDlohi_d,sznum);
   cudaMalloc((void**)&invDhilo_d,sznum);
   cudaMalloc((void**)&invDlolo_d,sznum);
   cudaMemcpy(Dhihi_d,Dhihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlohi_d,Dlohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dhilo_d,Dhilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlolo_d,Dlolo_h,sznum,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *invlapms = 0.0;
   *mullapms = 0.0;
   *sublapms = 0.0;
   *totlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   dbl4_invert_tiles<<<nbt,szt>>>
      (szt,Dhihi_d,   Dlohi_d,   Dhilo_d,   Dlolo_d,
        invDhihi_d,invDlohi_d,invDhilo_d,invDlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhshihi_d;                    // right hand side on device
   double *rhslohi_d;
   double *rhshilo_d;
   double *rhslolo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhshihi_d,szrhs);
   cudaMalloc((void**)&rhslohi_d,szrhs);
   cudaMalloc((void**)&rhshilo_d,szrhs);
   cudaMalloc((void**)&rhslolo_d,szrhs);
   cudaMemcpy(rhshihi_d,bhihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslohi_d,blohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhshilo_d,bhilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslolo_d,blolo,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl4_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDhihi_d,invDlohi_d,invDhilo_d,invDlolo_d,
                  rhshihi_d, rhslohi_d, rhshilo_d, rhslolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;             // #doubles in column of U
   double *Ucolhihi_h = new double[nbrUcol];  // column of U on host
   double *Ucollohi_h = new double[nbrUcol]; 
   double *Ucolhilo_h = new double[nbrUcol];
   double *Ucollolo_h = new double[nbrUcol];
   double *Ucolhihi_d;
   double *Ucollohi_d;
   double *Ucolhilo_d;
   double *Ucollolo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolhihi_d,szUcol);
   cudaMalloc((void**)&Ucollohi_d,szUcol);
   cudaMalloc((void**)&Ucolhilo_d,szUcol);
   cudaMalloc((void**)&Ucollolo_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++)
            {
               Ucolhihi_h[ix]   = Uhihi[rowoff+i][coloff+j];
               Ucollohi_h[ix]   = Ulohi[rowoff+i][coloff+j];
               Ucolhilo_h[ix]   = Uhilo[rowoff+i][coloff+j];
               Ucollolo_h[ix++] = Ulolo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolhihi_d,Ucolhihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollohi_d,Ucollohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolhilo_d,Ucolhilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollolo_d,Ucollolo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      dbl4_back_substitute<<<k,szt>>>
         (szt,k,Ucolhihi_d,Ucollohi_d,Ucolhilo_d,Ucollolo_d,
                 rhshihi_d, rhslohi_d, rhshilo_d, rhslolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_dbl_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      dbl4_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDhihi_d,invDlohi_d,invDhilo_d,invDlolo_d,
                   rhshihi_d, rhslohi_d, rhshilo_d, rhslolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *mullapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);

      nbrUcol = nbrUcol - szt*szt; // one tile less used in update
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(xhihi,rhshihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlohi,rhslohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xhilo,rhshilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlolo,rhslolo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDhihi_h,invDhihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlohi_h,invDlohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDhilo_h,invDhilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlolo_h,invDlolo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhihi[offset+i][offset+j] = invDhihi_h[ix];
            Ulohi[offset+i][offset+j] = invDlohi_h[ix];
            Uhilo[offset+i][offset+j] = invDhilo_h[ix];
            Ulolo[offset+i][offset+j] = invDlolo_h[ix++];
         }
   }
   free(Dhihi_h); free(Dlohi_h); free(Dhilo_h); free(Dlolo_h);
   free(invDhihi_h); free(invDlohi_h);
   free(invDhilo_h); free(invDlolo_h);
   free(Ucolhihi_h); free(Ucollohi_h);
   free(Ucolhilo_h); free(Ucollolo_h);
}

void GPU_cmplx4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;       // number of doubles on diagonal tiles
   double *Drehihi_h = new double[nbr];  // real parts of diagonal tiles
   double *Drelohi_h = new double[nbr];  
   double *Drehilo_h = new double[nbr];  
   double *Drelolo_h = new double[nbr]; 
   double *Dimhihi_h = new double[nbr];  // imaginary parts of diagonal tiles
   double *Dimlohi_h = new double[nbr];  
   double *Dimhilo_h = new double[nbr];  
   double *Dimlolo_h = new double[nbr];  
   double *Drehihi_d;                    // diagonal tiles on the device
   double *Drelohi_d;
   double *Drehilo_d;
   double *Drelolo_d;
   double *Dimhihi_d;
   double *Dimlohi_d;
   double *Dimhilo_d; 
   double *Dimlolo_d;
   double *invDrehihi_h = new double[nbr]; // real parts of inverse tiles
   double *invDrelohi_h = new double[nbr];
   double *invDrehilo_h = new double[nbr]; 
   double *invDrelolo_h = new double[nbr];
   double *invDimhihi_h = new double[nbr]; // imaginary parts of inverse tiles
   double *invDimlohi_h = new double[nbr];
   double *invDimhilo_h = new double[nbr];
   double *invDimlolo_h = new double[nbr];
   double *invDrehihi_d;           // invDrehihi_d ~ invDrehihi_h on device
   double *invDrelohi_d;           // invDrelohi_d ~ invDrelohi_h on device
   double *invDrehilo_d;           // invDrehilo_d ~ invDrehilo_h on device
   double *invDrelolo_d;           // invDrelolo_d ~ invDrelolo_h on device
   double *invDimhihi_d;           // invDimhihi_d ~ invDimhihi_h on device
   double *invDimlohi_d;           // invDimlohi_d ~ invDimlohi_h on device
   double *invDimhilo_d;           // invDimhilo_d ~ invDimhilo_h on device
   double *invDimlolo_d;           // invDimlolo_d ~ invDimlolo_h on device
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Drehihi_h[ix]   = Urehihi[offset+i][offset+j];
            Drelohi_h[ix]   = Urelohi[offset+i][offset+j];
            Drehilo_h[ix]   = Urehilo[offset+i][offset+j];
            Drelolo_h[ix]   = Urelolo[offset+i][offset+j];
            Dimhihi_h[ix]   = Uimhihi[offset+i][offset+j];
            Dimlohi_h[ix]   = Uimlohi[offset+i][offset+j];
            Dimhilo_h[ix]   = Uimhilo[offset+i][offset+j];
            Dimlolo_h[ix++] = Uimlolo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Drehihi_d,sznum);
   cudaMalloc((void**)&Drelohi_d,sznum);
   cudaMalloc((void**)&Drehilo_d,sznum);
   cudaMalloc((void**)&Drelolo_d,sznum);
   cudaMalloc((void**)&Dimhihi_d,sznum);
   cudaMalloc((void**)&Dimlohi_d,sznum);
   cudaMalloc((void**)&Dimhilo_d,sznum);
   cudaMalloc((void**)&Dimlolo_d,sznum);
   cudaMalloc((void**)&invDrehihi_d,sznum);
   cudaMalloc((void**)&invDrelohi_d,sznum);
   cudaMalloc((void**)&invDrehilo_d,sznum);
   cudaMalloc((void**)&invDrelolo_d,sznum);
   cudaMalloc((void**)&invDimhihi_d,sznum);
   cudaMalloc((void**)&invDimlohi_d,sznum);
   cudaMalloc((void**)&invDimhilo_d,sznum);
   cudaMalloc((void**)&invDimlolo_d,sznum);
   cudaMemcpy(Drehihi_d,Drehihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelohi_d,Drelohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drehilo_d,Drehilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelolo_d,Drelolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhihi_d,Dimhihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlohi_d,Dimlohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhilo_d,Dimhilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlolo_d,Dimlolo_h,sznum,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *invlapms = 0.0;
   *mullapms = 0.0;
   *sublapms = 0.0;
   *totlapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);
   cmplx4_invert_tiles<<<nbt,szt>>>
      (szt, Drehihi_d,     Drelohi_d,   Drehilo_d,   Drelolo_d,
            Dimhihi_d,     Dimlohi_d,   Dimhilo_d,   Dimlolo_d,
           invDrehihi_d,invDrelohi_d,invDrehilo_d,invDrelolo_d,
           invDimhihi_d,invDimlohi_d,invDimhilo_d,invDimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhsrehihi_d;                    // right hand side on device
   double *rhsrelohi_d;
   double *rhsrehilo_d;
   double *rhsrelolo_d;
   double *rhsimhihi_d;
   double *rhsimlohi_d;
   double *rhsimhilo_d;
   double *rhsimlolo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhsrehihi_d,szrhs);
   cudaMalloc((void**)&rhsrelohi_d,szrhs);
   cudaMalloc((void**)&rhsrehilo_d,szrhs);
   cudaMalloc((void**)&rhsrelolo_d,szrhs);
   cudaMalloc((void**)&rhsimhihi_d,szrhs);
   cudaMalloc((void**)&rhsimlohi_d,szrhs);
   cudaMalloc((void**)&rhsimhilo_d,szrhs);
   cudaMalloc((void**)&rhsimlolo_d,szrhs);
   cudaMemcpy(rhsrehihi_d,brehihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelohi_d,brelohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrehilo_d,brehilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelolo_d,brelolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhihi_d,bimhihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlohi_d,bimlohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhilo_d,bimhilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlolo_d,bimlolo,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   cmplx4_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDrehihi_d,invDrelohi_d,invDrehilo_d,invDrelolo_d,
                 invDimhihi_d,invDimlohi_d,invDimhilo_d,invDimlolo_d,
                  rhsrehihi_d, rhsrelohi_d, rhsrehilo_d, rhsrelolo_d,
                  rhsimhihi_d, rhsimlohi_d, rhsimhilo_d, rhsimlolo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;               // #doubles in column of U
   double *Ucolrehihi_h = new double[nbrUcol];  // column of U on host
   double *Ucolrelohi_h = new double[nbrUcol];
   double *Ucolrehilo_h = new double[nbrUcol];
   double *Ucolrelolo_h = new double[nbrUcol];
   double *Ucolimhihi_h = new double[nbrUcol];
   double *Ucolimlohi_h = new double[nbrUcol];
   double *Ucolimhilo_h = new double[nbrUcol];
   double *Ucolimlolo_h = new double[nbrUcol];
   double *Ucolrehihi_d;
   double *Ucolrelohi_d;
   double *Ucolrehilo_d;
   double *Ucolrelolo_d;
   double *Ucolimhihi_d;
   double *Ucolimlohi_d;
   double *Ucolimhilo_d;
   double *Ucolimlolo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolrehihi_d,szUcol);
   cudaMalloc((void**)&Ucolrelohi_d,szUcol);
   cudaMalloc((void**)&Ucolrehilo_d,szUcol);
   cudaMalloc((void**)&Ucolrelolo_d,szUcol);
   cudaMalloc((void**)&Ucolimhihi_d,szUcol);
   cudaMalloc((void**)&Ucolimlohi_d,szUcol);
   cudaMalloc((void**)&Ucolimhilo_d,szUcol);
   cudaMalloc((void**)&Ucolimlolo_d,szUcol);

   int coloff,rowoff;

   for(int k=nbt-1; k>0; k--)      // update with solution tile k
   {
      coloff = k*szt;      // column offset to update with solution tile k
      ix = 0;
      for(int L=0; L<k; L++)       // copy k tiles of U
      {
         rowoff = L*szt;           // row offset for update data
         for(int i=0; i<szt; i++)
            for(int j=0; j<szt; j++)
            {
               Ucolrehihi_h[ix]   = Urehihi[rowoff+i][coloff+j];
               Ucolrelohi_h[ix]   = Urelohi[rowoff+i][coloff+j];
               Ucolrehilo_h[ix]   = Urehilo[rowoff+i][coloff+j];
               Ucolrelolo_h[ix]   = Urelolo[rowoff+i][coloff+j];
               Ucolimhihi_h[ix]   = Uimhihi[rowoff+i][coloff+j];
               Ucolimlohi_h[ix]   = Uimlohi[rowoff+i][coloff+j];
               Ucolimhilo_h[ix]   = Uimhilo[rowoff+i][coloff+j];
               Ucolimlolo_h[ix++] = Uimlolo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolrehihi_d,Ucolrehihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelohi_d,Ucolrelohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrehilo_d,Ucolrehilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelolo_d,Ucolrelolo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhihi_d,Ucolimhihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlohi_d,Ucolimlohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhilo_d,Ucolimhilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlolo_d,Ucolimlolo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      cmplx4_back_substitute<<<k,szt>>>
         (szt,k,Ucolrehihi_d,Ucolrelohi_d,Ucolrehilo_d,Ucolrelolo_d,
                Ucolimhihi_d,Ucolimlohi_d,Ucolimhilo_d,Ucolimlolo_d,
                 rhsrehihi_d, rhsrelohi_d, rhsrehilo_d, rhsrelolo_d,
                 rhsimhihi_d, rhsimlohi_d, rhsimhilo_d, rhsimlolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_cmplx_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      cmplx4_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDrehihi_d,invDrelohi_d,invDrehilo_d,invDrelolo_d,
                  invDimhihi_d,invDimlohi_d,invDimhilo_d,invDimlolo_d,
                   rhsrehihi_d, rhsrelohi_d, rhsrehilo_d, rhsrelolo_d,
                   rhsimhihi_d, rhsimlohi_d, rhsimhilo_d, rhsimlolo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *mullapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);

      nbrUcol = nbrUcol - szt*szt; // one tile less used in update
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(xrehihi,rhsrehihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelohi,rhsrelohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrehilo,rhsrehilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelolo,rhsrelolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhihi,rhsimhihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlohi,rhsimlohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhilo,rhsimhilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlolo,rhsimlolo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDrehihi_h,invDrehihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelohi_h,invDrelohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrehilo_h,invDrehilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelolo_h,invDrelolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhihi_h,invDimhihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlohi_h,invDimlohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhilo_h,invDimhilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlolo_h,invDimlolo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Urehihi[offset+i][offset+j] = invDrehihi_h[ix];
            Urelohi[offset+i][offset+j] = invDrelohi_h[ix];
            Urehilo[offset+i][offset+j] = invDrehilo_h[ix];
            Urelolo[offset+i][offset+j] = invDrelolo_h[ix];
            Uimhihi[offset+i][offset+j] = invDimhihi_h[ix];
            Uimlohi[offset+i][offset+j] = invDimlohi_h[ix];
            Uimhilo[offset+i][offset+j] = invDimhilo_h[ix];
            Uimlolo[offset+i][offset+j] = invDimlolo_h[ix++];
         }
   }
   free(Drehihi_h); free(Drelohi_h); free(Drehilo_h); free(Drelolo_h);
   free(Dimhihi_h); free(Dimlohi_h); free(Dimhilo_h); free(Dimlolo_h);
   free(invDrehihi_h); free(invDrelohi_h);
   free(invDrehilo_h); free(invDrelolo_h);
   free(Ucolrehihi_h); free(Ucolrelohi_h);
   free(Ucolrehilo_h); free(Ucolrelolo_h);
   free(invDimhihi_h); free(invDimlohi_h);
   free(invDimhilo_h); free(invDimlolo_h);
   free(Ucolimhihi_h); free(Ucolimlohi_h);
   free(Ucolimhilo_h); free(Ucolimlolo_h);
}

*/
