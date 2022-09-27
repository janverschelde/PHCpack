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
   __syncthreads();
   invUrowsrehihihi[rowidx] = acc1hihihi;
   invUrowsrelohihi[rowidx] = acc1lohihi;
   invUrowsrehilohi[rowidx] = acc1hilohi;
   invUrowsrelolohi[rowidx] = acc1lolohi;
   invUrowsrehihilo[rowidx] = acc1hihilo;
   invUrowsrelohilo[rowidx] = acc1lohilo;
   invUrowsrehilolo[rowidx] = acc1hilolo;
   invUrowsrelololo[rowidx] = acc1lololo;
   __syncthreads();
   odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
           &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
            acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
            acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
   __syncthreads();
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
         __syncthreads();
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
         __syncthreads();
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
      __syncthreads();
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
      __syncthreads();
      invUrowsrehihihi[rowidx] = acc1hihihi;
      invUrowsrelohihi[rowidx] = acc1lohihi;
      invUrowsrehilohi[rowidx] = acc1hilohi;
      invUrowsrelolohi[rowidx] = acc1lolohi;
      invUrowsrehihilo[rowidx] = acc1hihilo;
      invUrowsrelohilo[rowidx] = acc1lohilo;
      invUrowsrehilolo[rowidx] = acc1hilolo;
      invUrowsrelololo[rowidx] = acc1lololo;
      __syncthreads();
      odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
              &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
               acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
               acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
      __syncthreads();
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

__global__ void cmplx8_medium_invert_upper
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
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolrehihihi[tabsod_shmemsize];  // one column of U
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
   __shared__ double invUrowrehihihi[tabsod_shmemsize]; // one row of invU
   __shared__ double invUrowrelohihi[tabsod_shmemsize]; 
   __shared__ double invUrowrehilohi[tabsod_shmemsize];
   __shared__ double invUrowrelolohi[tabsod_shmemsize]; 
   __shared__ double invUrowrehihilo[tabsod_shmemsize];
   __shared__ double invUrowrelohilo[tabsod_shmemsize]; 
   __shared__ double invUrowrehilolo[tabsod_shmemsize];
   __shared__ double invUrowrelololo[tabsod_shmemsize]; 
   __shared__ double invUrowimhihihi[tabsod_shmemsize]; 
   __shared__ double invUrowimlohihi[tabsod_shmemsize]; 
   __shared__ double invUrowimhilohi[tabsod_shmemsize]; 
   __shared__ double invUrowimlolohi[tabsod_shmemsize]; 
   __shared__ double invUrowimhihilo[tabsod_shmemsize]; 
   __shared__ double invUrowimlohilo[tabsod_shmemsize]; 
   __shared__ double invUrowimhilolo[tabsod_shmemsize]; 
   __shared__ double invUrowimlololo[tabsod_shmemsize]; 

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

   int colidx = dim*(dim-1);           // start with the last column

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
   int rowidx = (dim - 1)*dim + k;     // the row index in the inverse

   __syncthreads();
   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
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
   odg_minus(&invimhihihi,&invimlohihi,
             &invimhilohi,&invimlolohi,
             &invimhihilo,&invimlohilo,
             &invimhilolo,&invimlololo);
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
   invUrowrehihihi[k] = acc1hihihi;
   invUrowrelohihi[k] = acc1lohihi;
   invUrowrehilohi[k] = acc1hilohi;
   invUrowrelolohi[k] = acc1lolohi;
   invUrowrehihilo[k] = acc1hihilo;
   invUrowrelohilo[k] = acc1lohilo;
   invUrowrehilolo[k] = acc1hilolo;
   invUrowrelololo[k] = acc1lololo;
   odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
           &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
            acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
            acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
   invUrowimhihihi[k] = acc3hihihi;
   invUrowimlohihi[k] = acc3lohihi;
   invUrowimhilohi[k] = acc3hilohi;
   invUrowimlolohi[k] = acc3lolohi;
   invUrowimhihilo[k] = acc3hihilo;
   invUrowimlohilo[k] = acc3lohilo;
   invUrowimhilolo[k] = acc3hilolo;
   invUrowimlololo[k] = acc3lololo;
   invUrehihihi[rowidx] = invUrowrehihihi[k];  // store last row into invU
   invUrelohihi[rowidx] = invUrowrelohihi[k]; 
   invUrehilohi[rowidx] = invUrowrehilohi[k]; 
   invUrelolohi[rowidx] = invUrowrelolohi[k]; 
   invUrehihilo[rowidx] = invUrowrehihilo[k];
   invUrelohilo[rowidx] = invUrowrelohilo[k]; 
   invUrehilolo[rowidx] = invUrowrehilolo[k]; 
   invUrelololo[rowidx] = invUrowrelololo[k]; 
   invUimhihihi[rowidx] = invUrowimhihihi[k];
   invUimlohihi[rowidx] = invUrowimlohihi[k]; 
   invUimhilohi[rowidx] = invUrowimhilohi[k];
   invUimlolohi[rowidx] = invUrowimlolohi[k]; 
   invUimhihilo[rowidx] = invUrowimhihilo[k];
   invUimlohilo[rowidx] = invUrowimlohilo[k]; 
   invUimhilolo[rowidx] = invUrowimhilolo[k];
   invUimlololo[rowidx] = invUrowimlololo[k]; 

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

         rowidx = j*dim + k;                      // need solution value
         invUrowrehihihi[k] = invUrehihihi[rowidx]; // load invU row
         invUrowrelohihi[k] = invUrelohihi[rowidx];
         invUrowrehilohi[k] = invUrehilohi[rowidx]; 
         invUrowrelolohi[k] = invUrelolohi[rowidx];
         invUrowrehihilo[k] = invUrehihilo[rowidx];
         invUrowrelohilo[k] = invUrelohilo[rowidx];
         invUrowrehilolo[k] = invUrehilolo[rowidx]; 
         invUrowrelololo[k] = invUrelololo[rowidx];
         invUrowimhihihi[k] = invUimhihihi[rowidx];
         invUrowimlohihi[k] = invUimlohihi[rowidx];
         invUrowimhilohi[k] = invUimhilohi[rowidx];
         invUrowimlolohi[k] = invUimlolohi[rowidx];
         invUrowimhihilo[k] = invUimhihilo[rowidx];
         invUrowimlohilo[k] = invUimlohilo[rowidx];
         invUrowimhilolo[k] = invUimhilolo[rowidx];
         invUrowimlololo[k] = invUimlololo[rowidx];
         xvalrehihihi = invUrowrehihihi[k];
         xvalrelohihi = invUrowrelohihi[k];
         xvalrehilohi = invUrowrehilohi[k];
         xvalrelolohi = invUrowrelolohi[k];
         xvalrehihilo = invUrowrehihilo[k];
         xvalrelohilo = invUrowrelohilo[k];
         xvalrehilolo = invUrowrehilolo[k];
         xvalrelololo = invUrowrelololo[k];
         xvalimhihihi = invUrowimhihihi[k];
         xvalimlohihi = invUrowimlohihi[k];
         xvalimhilohi = invUrowimhilohi[k];
         xvalimlolohi = invUrowimlolohi[k];
         xvalimhihilo = invUrowimhihilo[k];
         xvalimlohilo = invUrowimlohilo[k];
         xvalimhilolo = invUrowimhilolo[k];
         xvalimlololo = invUrowimlololo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                 xvalrehihihi,   xvalrelohihi,   xvalrehilohi,   xvalrelolohi,
                 xvalrehihilo,   xvalrelohilo,   xvalrehilolo,   xvalrelololo,
                  &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
                  &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
         odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
                 xvalimhihihi,   xvalimlohihi,   xvalimhilohi,   xvalimlolohi,
                 xvalimhihilo,   xvalimlohilo,   xvalimhilolo,   xvalimlololo,
                  &acc2hihihi,    &acc2lohihi,    &acc2hilohi,    &acc2lolohi,
                  &acc2hihilo,    &acc2lohilo,    &acc2hilolo,    &acc2lololo);
         odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
                 xvalrehihihi,   xvalrelohihi,   xvalrehilohi,   xvalrelolohi,
                 xvalrehihilo,   xvalrelohilo,   xvalrehilolo,   xvalrelololo,
                  &acc3hihihi,    &acc3lohihi,    &acc3hilohi,    &acc3lolohi,
                  &acc3hihilo,    &acc3lohilo,    &acc3hilolo,    &acc3lololo);
         odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                 xvalimhihihi,   xvalimlohihi,   xvalimhilohi,   xvalimlolohi,
                 xvalimhihilo,   xvalimlohilo,   xvalimhilolo,   xvalimlololo,
                  &acc4hihihi,    &acc4lohihi,    &acc4hilohi,    &acc4lolohi,
                  &acc4hihilo,    &acc4lohilo,    &acc4hilolo,    &acc4lololo);
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
      Ucolimhilohi[k] = Uimhilohi[colidx+k];
      Ucolimlolohi[k] = Uimlolohi[colidx+k];
      Ucolimhihilo[k] = Uimhihilo[colidx+k];
      Ucolimlohilo[k] = Uimlohilo[colidx+k];
      Ucolimhilolo[k] = Uimhilolo[colidx+k];
      Ucolimlololo[k] = Uimlololo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
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
      odg_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
              &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
               acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
               acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
      invUrowrehihihi[k] = acc1hihihi;
      invUrowrelohihi[k] = acc1lohihi;
      invUrowrehilohi[k] = acc1hilohi;
      invUrowrelolohi[k] = acc1lolohi;
      invUrowrehihilo[k] = acc1hihilo;
      invUrowrelohilo[k] = acc1lohilo;
      invUrowrehilolo[k] = acc1hilolo;
      invUrowrelololo[k] = acc1lololo;
      odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
              &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
               acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
               acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
      invUrowimhihihi[k] = acc3hihihi;
      invUrowimlohihi[k] = acc3lohihi;
      invUrowimhilohi[k] = acc3hilohi;
      invUrowimlolohi[k] = acc3lolohi;
      invUrowimhihilo[k] = acc3hihilo;
      invUrowimlohilo[k] = acc3lohilo;
      invUrowimhilolo[k] = acc3hilolo;
      invUrowimlololo[k] = acc3lololo;
      invUrehihihi[rowidx] = invUrowrehihihi[k];
      invUrelohihi[rowidx] = invUrowrelohihi[k];
      invUrehilohi[rowidx] = invUrowrehilohi[k];
      invUrelolohi[rowidx] = invUrowrelolohi[k];
      invUrehihilo[rowidx] = invUrowrehihilo[k];
      invUrelohilo[rowidx] = invUrowrelohilo[k];
      invUrehilolo[rowidx] = invUrowrehilolo[k];
      invUrelololo[rowidx] = invUrowrelololo[k];
      invUimhihihi[rowidx] = invUrowimhihihi[k];
      invUimlohihi[rowidx] = invUrowimlohihi[k];
      invUimhilohi[rowidx] = invUrowimhilohi[k];
      invUimlolohi[rowidx] = invUrowimlolohi[k];
      invUimhihilo[rowidx] = invUrowimhihilo[k];
      invUimlohilo[rowidx] = invUrowimlohilo[k];
      invUimhilolo[rowidx] = invUrowimhilolo[k];
      invUimlololo[rowidx] = invUrowimlololo[k];
   }
}

__global__ void  dbl8_invert_tiles
 ( int dim,
   double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi,
   double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo,
   double *invUhihihi, double *invUlohihi,
   double *invUhilohi, double *invUlolohi,
   double *invUhihilo, double *invUlohilo,
   double *invUhilolo, double *invUlololo )
{
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

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

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolhihihi[k] = Uhihihi[colidx+k];  // load the last column
   Ucollohihi[k] = Ulohihi[colidx+k];
   Ucolhilohi[k] = Uhilohi[colidx+k]; 
   Ucollolohi[k] = Ulolohi[colidx+k];
   Ucolhihilo[k] = Uhihilo[colidx+k];
   Ucollohilo[k] = Ulohilo[colidx+k];
   Ucolhilolo[k] = Uhilolo[colidx+k]; 
   Ucollololo[k] = Ulololo[colidx+k];
   rhshihihi = ((double) int(k == dim-1));  // right hand side for threads
   rhslohihi = 0.0;
   rhshilohi = 0.0;
   rhslolohi = 0.0;
   rhshihilo = 0.0;
   rhslohilo = 0.0;
   rhshilolo = 0.0;
   rhslololo = 0.0;
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
   __syncthreads();

   invUhihihi[rowidx] = 0.0; // initialize in case of zero divisor
   invUlohihi[rowidx] = 0.0;
   invUhilohi[rowidx] = 0.0;
   invUlolohi[rowidx] = 0.0;
   invUhihilo[rowidx] = 0.0;
   invUlohilo[rowidx] = 0.0;
   invUhilolo[rowidx] = 0.0;
   invUlololo[rowidx] = 0.0;

   if(1.0 + Ucolhihihi[k] != 1.0)
   {
      odg_div(     rhshihihi,        rhslohihi,    rhshilohi,    rhslolohi,
                   rhshihilo,        rhslohilo,    rhshilolo,    rhslololo,
                  Ucolhihihi[k],    Ucollohihi[k],Ucolhilohi[k],Ucollolohi[k],
                  Ucolhihilo[k],    Ucollohilo[k],Ucolhilolo[k],Ucollololo[k],
              &invUrowhihihi[k],&invUrowlohihi[k],
              &invUrowhilohi[k],&invUrowlolohi[k],
              &invUrowhihilo[k],&invUrowlohilo[k],
              &invUrowhilolo[k],&invUrowlololo[k]);
      // store the last row into invU
      invUhihihi[rowidx] = invUrowhihihi[k];
      invUlohihi[rowidx] = invUrowlohihi[k];
      invUhilohi[rowidx] = invUrowhilohi[k];
      invUlolohi[rowidx] = invUrowlolohi[k];
      invUhihilo[rowidx] = invUrowhihilo[k];
      invUlohilo[rowidx] = invUrowlohilo[k];
      invUhilolo[rowidx] = invUrowhilolo[k];
      invUlololo[rowidx] = invUrowlololo[k];
   }
   __syncthreads();
   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshihihi = ((double) int(k == i));   // set rhs for i-th unit vector
      rhslohihi = 0.0;
      rhshilohi = 0.0;
      rhslolohi = 0.0;
      rhshihilo = 0.0;
      rhslohilo = 0.0;
      rhshilolo = 0.0;
      rhslololo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = offset + dim*j;     // need column j of U
         __syncthreads();
         Ucolhihihi[k] = Uhihihi[colidx+k];
         Ucollohihi[k] = Ulohihi[colidx+k];
         Ucolhilohi[k] = Uhilohi[colidx+k];
         Ucollolohi[k] = Ulolohi[colidx+k];
         Ucolhihilo[k] = Uhihilo[colidx+k];
         Ucollohilo[k] = Ulohilo[colidx+k];
         Ucolhilolo[k] = Uhilolo[colidx+k];
         Ucollololo[k] = Ulololo[colidx+k];

         rowidx = offset + j*dim + k;       // need solution value
         __syncthreads();
         invUrowhihihi[k] = invUhihihi[rowidx]; // load invU row into invUrow
         invUrowlohihi[k] = invUlohihi[rowidx];
         invUrowhilohi[k] = invUhilohi[rowidx];
         invUrowlolohi[k] = invUlolohi[rowidx];
         invUrowhihilo[k] = invUhihilo[rowidx];
         invUrowlohilo[k] = invUlohilo[rowidx];
         invUrowhilolo[k] = invUhilolo[rowidx];
         invUrowlololo[k] = invUlololo[rowidx];
         __syncthreads();
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
      colidx = offset + dim*i;        // need column i of U
      __syncthreads();
      Ucolhihihi[k] = Uhihihi[colidx+k];
      Ucollohihi[k] = Ulohihi[colidx+k];
      Ucolhilohi[k] = Uhilohi[colidx+k];
      Ucollolohi[k] = Ulolohi[colidx+k];
      Ucolhihilo[k] = Uhihilo[colidx+k];
      Ucollohilo[k] = Ulohilo[colidx+k];
      Ucolhilolo[k] = Uhilolo[colidx+k];
      Ucollololo[k] = Ulololo[colidx+k];
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      invUhihihi[rowidx] = 0.0;    // initialize in case of zero divisor
      invUlohihi[rowidx] = 0.0;
      invUhilohi[rowidx] = 0.0;
      invUlolohi[rowidx] = 0.0;
      invUhihilo[rowidx] = 0.0;
      invUlohilo[rowidx] = 0.0;
      invUhilolo[rowidx] = 0.0;
      invUlololo[rowidx] = 0.0;

      if(1.0 + Ucolhihihi[i] != 1.0)
      {
         odg_div(  rhshihihi,        rhslohihi,    rhshilohi,    rhslolohi,
                   rhshihilo,        rhslohilo,    rhshilolo,    rhslololo,
                  Ucolhihihi[i],    Ucollohihi[i],Ucolhilohi[i],Ucollolohi[i],
                  Ucolhihilo[i],    Ucollohilo[i],Ucolhilolo[i],Ucollololo[i],
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
}

__global__ void cmplx8_invert_tiles
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
   const int B = blockIdx.x;   // block index
   const int k = threadIdx.x;  // thread k computes k-th column of inverse
   const int offset = dim*dim*B; // offset in U and invU

   __shared__ double Ucolrehihihi[tabsod_shmemsize]; // one column of U
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
   __shared__ double invUrowrehihihi[tabsod_shmemsize]; // one row of invU
   __shared__ double invUrowrelohihi[tabsod_shmemsize]; 
   __shared__ double invUrowrehilohi[tabsod_shmemsize]; 
   __shared__ double invUrowrelolohi[tabsod_shmemsize]; 
   __shared__ double invUrowrehihilo[tabsod_shmemsize];
   __shared__ double invUrowrelohilo[tabsod_shmemsize]; 
   __shared__ double invUrowrehilolo[tabsod_shmemsize]; 
   __shared__ double invUrowrelololo[tabsod_shmemsize]; 
   __shared__ double invUrowimhihihi[tabsod_shmemsize];
   __shared__ double invUrowimlohihi[tabsod_shmemsize]; 
   __shared__ double invUrowimhilohi[tabsod_shmemsize];
   __shared__ double invUrowimlolohi[tabsod_shmemsize]; 
   __shared__ double invUrowimhihilo[tabsod_shmemsize];
   __shared__ double invUrowimlohilo[tabsod_shmemsize]; 
   __shared__ double invUrowimhilolo[tabsod_shmemsize];
   __shared__ double invUrowimlololo[tabsod_shmemsize]; 

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

   int colidx = offset + dim*(dim-1); // start with the last column

   Ucolrehihihi[k] = Urehihihi[colidx+k];   // load the last column
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
   int rowidx = offset + (dim - 1)*dim + k; // row index in the inverse

   // invUrow[k] = rhs/Ucol[k];       // last row of the inverse
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

   invUrehihihi[rowidx] = 0.0;   // initialize in case of zero denominator
   invUrelohihi[rowidx] = 0.0;
   invUrehilohi[rowidx] = 0.0; 
   invUrelolohi[rowidx] = 0.0;
   invUrehihilo[rowidx] = 0.0;
   invUrelohilo[rowidx] = 0.0;
   invUrehilolo[rowidx] = 0.0; 
   invUrelololo[rowidx] = 0.0;
   invUimhihihi[rowidx] = 0.0;
   invUimlohihi[rowidx] = 0.0;
   invUimhilohi[rowidx] = 0.0;
   invUimlolohi[rowidx] = 0.0; 
   invUimhihilo[rowidx] = 0.0;
   invUimlohilo[rowidx] = 0.0;
   invUimhilolo[rowidx] = 0.0;
   invUimlololo[rowidx] = 0.0;

   if(1.0 + denhihihi != 1.0)
   {
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
      // __syncthreads();
      invUrowrehihihi[k] = acc1hihihi;
      invUrowrelohihi[k] = acc1lohihi;
      invUrowrehilohi[k] = acc1hilohi;
      invUrowrelolohi[k] = acc1lolohi;
      invUrowrehihilo[k] = acc1hihilo;
      invUrowrelohilo[k] = acc1lohilo;
      invUrowrehilolo[k] = acc1hilolo;
      invUrowrelololo[k] = acc1lololo;
      // __syncthreads();
      odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
              &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
               acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
               acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
      // __syncthreads();
      invUrowimhihihi[k] = acc3hihihi;
      invUrowimlohihi[k] = acc3lohihi;
      invUrowimhilohi[k] = acc3hilohi;
      invUrowimlolohi[k] = acc3lolohi;
      invUrowimhihilo[k] = acc3hihilo;
      invUrowimlohilo[k] = acc3lohilo;
      invUrowimhilolo[k] = acc3hilolo;
      invUrowimlololo[k] = acc3lololo;
      // __syncthreads();
      invUrehihihi[rowidx] = invUrowrehihihi[k];   // store the last row
      invUrelohihi[rowidx] = invUrowrelohihi[k];
      invUrehilohi[rowidx] = invUrowrehilohi[k]; 
      invUrelolohi[rowidx] = invUrowrelolohi[k];
      invUrehihilo[rowidx] = invUrowrehihilo[k];
      invUrelohilo[rowidx] = invUrowrelohilo[k];
      invUrehilolo[rowidx] = invUrowrehilolo[k]; 
      invUrelololo[rowidx] = invUrowrelololo[k];
      invUimhihihi[rowidx] = invUrowimhihihi[k];
      invUimlohihi[rowidx] = invUrowimlohihi[k];
      invUimhilohi[rowidx] = invUrowimhilohi[k];
      invUimlolohi[rowidx] = invUrowimlolohi[k];
      invUimhihilo[rowidx] = invUrowimhihilo[k];
      invUimlohilo[rowidx] = invUrowimlohilo[k];
      invUimhilolo[rowidx] = invUrowimhilolo[k];
      invUimlololo[rowidx] = invUrowimlololo[k];
   }
   __syncthreads();
   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhsrehihihi = ((double) int(k == i));   // set rhs for i-th unit vector
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
         colidx = offset + dim*j;        // need column j of U
         __syncthreads();
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

         rowidx = offset + j*dim + k;       // need solution value
         __syncthreads();
         invUrowrehihihi[k] = invUrehihihi[rowidx]; // load invU row
         invUrowrelohihi[k] = invUrelohihi[rowidx];
         invUrowrehilohi[k] = invUrehilohi[rowidx];
         invUrowrelolohi[k] = invUrelolohi[rowidx];
         invUrowrehihilo[k] = invUrehihilo[rowidx];
         invUrowrelohilo[k] = invUrelohilo[rowidx];
         invUrowrehilolo[k] = invUrehilolo[rowidx];
         invUrowrelololo[k] = invUrelololo[rowidx];
         invUrowimhihihi[k] = invUimhihihi[rowidx];
         invUrowimlohihi[k] = invUimlohihi[rowidx];
         invUrowimhilohi[k] = invUimhilohi[rowidx];
         invUrowimlolohi[k] = invUimlolohi[rowidx];
         invUrowimhihilo[k] = invUimhihilo[rowidx];
         invUrowimlohilo[k] = invUimlohilo[rowidx];
         invUrowimhilolo[k] = invUimhilolo[rowidx];
         invUrowimlololo[k] = invUimlololo[rowidx];
         __syncthreads();
         xvalrehihihi = invUrowrehihihi[k];
         xvalrelohihi = invUrowrelohihi[k];
         xvalrehilohi = invUrowrehilohi[k];
         xvalrelolohi = invUrowrelolohi[k];
         xvalrehihilo = invUrowrehihilo[k];
         xvalrelohilo = invUrowrelohilo[k];
         xvalrehilolo = invUrowrehilolo[k];
         xvalrelololo = invUrowrelololo[k];
         xvalimhihihi = invUrowimhihihi[k];
         xvalimlohihi = invUrowimlohihi[k];
         xvalimhilohi = invUrowimhilohi[k];
         xvalimlolohi = invUrowimlolohi[k];
         xvalimhihilo = invUrowimhihilo[k];
         xvalimlohilo = invUrowimlohilo[k];
         xvalimhilolo = invUrowimhilolo[k];
         xvalimlololo = invUrowimlololo[k];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval;    // update right hand side
         odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                 xvalrehihihi,   xvalrelohihi,   xvalrehilohi, xvalrelolohi,
                 xvalrehihilo,   xvalrelohilo,   xvalrehilolo, xvalrelololo,
                  &acc1hihihi,    &acc1lohihi,    &acc1hilohi,  &acc1lolohi,
                  &acc1hihilo,    &acc1lohilo,    &acc1hilolo,  &acc1lololo);
         odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
                 xvalimhihihi,   xvalimlohihi,   xvalimhilohi, xvalimlolohi,
                 xvalimhihilo,   xvalimlohilo,   xvalimhilolo, xvalimlololo,
                  &acc2hihihi,    &acc2lohihi,    &acc2hilohi,  &acc2lolohi,
                  &acc2hihilo,    &acc2lohilo,    &acc2hilolo,  &acc2lololo);
         odg_mul(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
                 xvalrehihihi,   xvalrelohihi,   xvalrehilohi, xvalrelolohi,
                 xvalrehihilo,   xvalrelohilo,   xvalrehilolo, xvalrelololo,
                  &acc3hihihi,    &acc3lohihi,    &acc3hilohi,  &acc3lolohi,
                  &acc3hihilo,    &acc3lohilo,    &acc3hilolo,  &acc3lololo);
         odg_mul(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                 xvalimhihihi,   xvalimlohihi,   xvalimhilohi, xvalimlolohi,
                 xvalimhihilo,   xvalimlohilo,   xvalimhilolo, xvalimlololo,
                  &acc4hihihi,    &acc4lohihi,    &acc4hilohi,  &acc4lolohi,
                  &acc4hihilo,    &acc4lohilo,    &acc4hilolo,  &acc4lololo);
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
      colidx = offset + dim*i;        // need column i of U
      __syncthreads();
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
      rowidx = offset + i*dim + k;    // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
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
      odg_inc(&denhihihi,&denlohihi,&denhilohi,&denlolohi,
              &denhihilo,&denlohilo,&denhilolo,&denlololo,
              acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
              acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);

      invUrehihihi[rowidx] = 0.0; // initialize in case
      invUrelohihi[rowidx] = 0.0; // of zero denominator
      invUrehilohi[rowidx] = 0.0;
      invUrelolohi[rowidx] = 0.0;
      invUrehihilo[rowidx] = 0.0;
      invUrelohilo[rowidx] = 0.0;
      invUrehilolo[rowidx] = 0.0;
      invUrelololo[rowidx] = 0.0;
      invUimhihihi[rowidx] = 0.0;
      invUimlohihi[rowidx] = 0.0;
      invUimhilohi[rowidx] = 0.0;
      invUimlolohi[rowidx] = 0.0;
      invUimhihilo[rowidx] = 0.0;
      invUimlohilo[rowidx] = 0.0;
      invUimhilolo[rowidx] = 0.0;
      invUimlololo[rowidx] = 0.0;

      if(1.0 + denhihihi != 1.0)
      {
         odg_div(Ucolrehihihi[i],Ucolrelohihi[i],
                 Ucolrehilohi[i],Ucolrelolohi[i],
                 Ucolrehihilo[i],Ucolrelohilo[i],
                 Ucolrehilolo[i],Ucolrelololo[i],
                    denhihihi,      denlohihi,      denhilohi,      denlolohi,
                    denhihilo,      denlohilo,      denhilolo,      denlololo,
                 &invrehihihi,   &invrelohihi,   &invrehilohi,   &invrelolohi,
                 &invrehihilo,   &invrelohilo,   &invrehilolo,   &invrelololo);
         odg_div(Ucolimhihihi[i],Ucolimlohihi[i],
                 Ucolimhilohi[i],Ucolimlolohi[i],
                 Ucolimhihilo[i],Ucolimlohilo[i],
                 Ucolimhilolo[i],Ucolimlololo[i],
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
      }
      __syncthreads(); // moved synchronization outside if statements
      if(1.0 + denhihihi != 1.0)
      {
         invUrowrehihihi[k] = acc1hihihi;
         invUrowrelohihi[k] = acc1lohihi;
         invUrowrehilohi[k] = acc1hilohi;
         invUrowrelolohi[k] = acc1lolohi;
         invUrowrehihilo[k] = acc1hihilo;
         invUrowrelohilo[k] = acc1lohilo;
         invUrowrehilolo[k] = acc1hilolo;
         invUrowrelololo[k] = acc1lololo;
      }
      __syncthreads();
      if(1.0 + denhihihi != 1.0)
      {
         odg_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                 &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
                  acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
                  acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
      }
      __syncthreads();
      if(1.0 + denhihihi != 1.0)
      {
         invUrowimhihihi[k] = acc3hihihi;
         invUrowimlohihi[k] = acc3lohihi;
         invUrowimhilohi[k] = acc3hilohi;
         invUrowimlolohi[k] = acc3lolohi;
         invUrowimhihilo[k] = acc3hihilo;
         invUrowimlohilo[k] = acc3lohilo;
         invUrowimhilolo[k] = acc3hilolo;
         invUrowimlololo[k] = acc3lololo;
      }
      __syncthreads();
      if(1.0 + denhihihi != 1.0)
      {
         invUrehihihi[rowidx] = invUrowrehihihi[k];
         invUrelohihi[rowidx] = invUrowrelohihi[k];
         invUrehilohi[rowidx] = invUrowrehilohi[k];
         invUrelolohi[rowidx] = invUrowrelolohi[k];
         invUrehihilo[rowidx] = invUrowrehihilo[k];
         invUrelohilo[rowidx] = invUrowrelohilo[k];
         invUrehilolo[rowidx] = invUrowrehilolo[k];
         invUrelololo[rowidx] = invUrowrelololo[k];
         invUimhihihi[rowidx] = invUrowimhihihi[k];
         invUimlohihi[rowidx] = invUrowimlohihi[k];
         invUimhilohi[rowidx] = invUrowimhilohi[k];
         invUimlolohi[rowidx] = invUrowimlolohi[k];
         invUimhihilo[rowidx] = invUrowimhihilo[k];
         invUimlohilo[rowidx] = invUrowimlohilo[k];
         invUimhilolo[rowidx] = invUrowimhilolo[k];
         invUimlololo[rowidx] = invUrowimlololo[k];
      }
   }
}

__global__ void dbl8_multiply_inverse
 ( int dim, int idx,
   double *invUhihihi, double *invUlohihi,
   double *invUhilohi, double *invUlolohi,
   double *invUhihilo, double *invUlohilo,
   double *invUhilolo, double *invUlololo,
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workhihihi[tabsod_shmemsize];    // copy of w
   __shared__ double worklohihi[tabsod_shmemsize];
   __shared__ double workhilohi[tabsod_shmemsize];
   __shared__ double worklolohi[tabsod_shmemsize];
   __shared__ double workhihilo[tabsod_shmemsize];
   __shared__ double worklohilo[tabsod_shmemsize];
   __shared__ double workhilolo[tabsod_shmemsize];
   __shared__ double worklololo[tabsod_shmemsize];

   double resulthihihi = 0.0; // each thread stores its product in result
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;
   double coeffhihihi,coefflohihi,coeffhilohi,coefflolohi;
   double coeffhihilo,coefflohilo,coeffhilolo,coefflololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   workhihihi[k] = whihihi[rhsoff+k];
   worklohihi[k] = wlohihi[rhsoff+k];
   workhilohi[k] = whilohi[rhsoff+k];
   worklolohi[k] = wlolohi[rhsoff+k];
   workhihilo[k] = whihilo[rhsoff+k];
   worklohilo[k] = wlohilo[rhsoff+k];
   workhilolo[k] = whilolo[rhsoff+k];
   worklololo[k] = wlololo[rhsoff+k];

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      __syncthreads();
      coeffhihihi = invUhihihi[offset+k*dim+j]; // thread k does row k
      coefflohihi = invUlohihi[offset+k*dim+j];
      coeffhilohi = invUhilohi[offset+k*dim+j];
      coefflolohi = invUlolohi[offset+k*dim+j];
      coeffhihilo = invUhihilo[offset+k*dim+j];
      coefflohilo = invUlohilo[offset+k*dim+j];
      coeffhilolo = invUhilolo[offset+k*dim+j];
      coefflololo = invUlololo[offset+k*dim+j];
      // result = result + coeff*work[j];
      __syncthreads();
      odg_mul(coeffhihihi,  coefflohihi,  coeffhilohi,  coefflolohi,
              coeffhihilo,  coefflohilo,  coeffhilolo,  coefflololo,
               workhihihi[j],worklohihi[j],workhilohi[j],worklolohi[j],
               workhihilo[j],worklohilo[j],workhilolo[j],worklololo[j],
               &acchihihi,   &acclohihi,   &acchilohi,   &acclolohi,
               &acchihilo,   &acclohilo,   &acchilolo,   &acclololo);
      __syncthreads();
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo);
   }
   __syncthreads();
   whihihi[rhsoff+k] = resulthihihi;
   wlohihi[rhsoff+k] = resultlohihi;
   whilohi[rhsoff+k] = resulthilohi;
   wlolohi[rhsoff+k] = resultlolohi;
   whihilo[rhsoff+k] = resulthihilo;
   wlohilo[rhsoff+k] = resultlohilo;
   whilolo[rhsoff+k] = resulthilolo;
   wlololo[rhsoff+k] = resultlololo;
}

__global__ void cmplx8_multiply_inverse
 ( int dim, int idx,
   double *invUrehihihi, double *invUrelohihi,
   double *invUrehilohi, double *invUrelolohi,
   double *invUrehihilo, double *invUrelohilo,
   double *invUrehilolo, double *invUrelololo,
   double *invUimhihihi, double *invUimlohihi,
   double *invUimhilohi, double *invUimlolohi,
   double *invUimhihilo, double *invUimlohilo,
   double *invUimhilolo, double *invUimlololo,
   double *wrehihihi, double *wrelohihi,
   double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo,
   double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
   const int k = threadIdx.x;     // thread k computes k-th product
   const int rhsoff = dim*idx;    // offset for the right hand size
   const int offset = dim*rhsoff; // offset for diagonal tile

   __shared__ double workrehihihi[tabsod_shmemsize];  // copy of w
   __shared__ double workrelohihi[tabsod_shmemsize]; 
   __shared__ double workrehilohi[tabsod_shmemsize];
   __shared__ double workrelolohi[tabsod_shmemsize]; 
   __shared__ double workrehihilo[tabsod_shmemsize];
   __shared__ double workrelohilo[tabsod_shmemsize]; 
   __shared__ double workrehilolo[tabsod_shmemsize];
   __shared__ double workrelololo[tabsod_shmemsize]; 
   __shared__ double workimhihihi[tabsod_shmemsize];
   __shared__ double workimlohihi[tabsod_shmemsize];
   __shared__ double workimhilohi[tabsod_shmemsize];
   __shared__ double workimlolohi[tabsod_shmemsize];
   __shared__ double workimhihilo[tabsod_shmemsize];
   __shared__ double workimlohilo[tabsod_shmemsize];
   __shared__ double workimhilolo[tabsod_shmemsize];
   __shared__ double workimlololo[tabsod_shmemsize];

   double resultrehihihi = 0.0; // each thread stores its product in result
   double resultrelohihi = 0.0;
   double resultrehilohi = 0.0;
   double resultrelolohi = 0.0;
   double resultrehihilo = 0.0;
   double resultrelohilo = 0.0;
   double resultrehilolo = 0.0;
   double resultrelololo = 0.0;
   double resultimhihihi = 0.0;
   double resultimlohihi = 0.0;
   double resultimhilohi = 0.0;
   double resultimlolohi = 0.0;
   double resultimhihilo = 0.0;
   double resultimlohilo = 0.0;
   double resultimhilolo = 0.0;
   double resultimlololo = 0.0;
   double coeffrehihihi,coeffrelohihi,coeffrehilohi,coeffrelolohi;
   double coeffrehihilo,coeffrelohilo,coeffrehilolo,coeffrelololo;
   double coeffimhihihi,coeffimlohihi,coeffimhilohi,coeffimlolohi;
   double coeffimhihilo,coeffimlohilo,coeffimhilolo,coeffimlololo;
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;

   workrehihihi[k] = wrehihihi[rhsoff+k];
   workrelohihi[k] = wrelohihi[rhsoff+k];
   workrehilohi[k] = wrehilohi[rhsoff+k];
   workrelolohi[k] = wrelolohi[rhsoff+k];
   workrehihilo[k] = wrehihilo[rhsoff+k];
   workrelohilo[k] = wrelohilo[rhsoff+k];
   workrehilolo[k] = wrehilolo[rhsoff+k];
   workrelololo[k] = wrelololo[rhsoff+k];
   workimhihihi[k] = wimhihihi[rhsoff+k];
   workimlohihi[k] = wimlohihi[rhsoff+k];
   workimhilohi[k] = wimhilohi[rhsoff+k];
   workimlolohi[k] = wimlolohi[rhsoff+k];
   workimhihilo[k] = wimhihilo[rhsoff+k];
   workimlohilo[k] = wimlohilo[rhsoff+k];
   workimhilolo[k] = wimhilolo[rhsoff+k];
   workimlololo[k] = wimlololo[rhsoff+k];

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffrehihihi = invUrehihihi[offset+k*dim+j]; // thread k does row k
      coeffrelohihi = invUrelohihi[offset+k*dim+j];
      coeffrehilohi = invUrehilohi[offset+k*dim+j];
      coeffrelolohi = invUrelolohi[offset+k*dim+j];
      coeffrehihilo = invUrehihilo[offset+k*dim+j];
      coeffrelohilo = invUrelohilo[offset+k*dim+j];
      coeffrehilolo = invUrehilolo[offset+k*dim+j];
      coeffrelololo = invUrelololo[offset+k*dim+j];
      coeffimhihihi = invUimhihihi[offset+k*dim+j];
      coeffimlohihi = invUimlohihi[offset+k*dim+j];
      coeffimhilohi = invUimhilohi[offset+k*dim+j];
      coeffimlolohi = invUimlolohi[offset+k*dim+j];
      coeffimhihilo = invUimhihilo[offset+k*dim+j];
      coeffimlohilo = invUimlohilo[offset+k*dim+j];
      coeffimhilolo = invUimhilolo[offset+k*dim+j];
      coeffimlololo = invUimlololo[offset+k*dim+j];
      // result = result + coeff*work[j];
      odg_mul(coeffrehihihi,  coeffrelohihi,  coeffrehilohi,  coeffrelolohi,
              coeffrehihilo,  coeffrelohilo,  coeffrehilolo,  coeffrelololo,
               workrehihihi[j],workrelohihi[j],workrehilohi[j],workrelolohi[j],
               workrehihilo[j],workrelohilo[j],workrehilolo[j],workrelololo[j],
                &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
                &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
      odg_mul(coeffimhihihi,  coeffimlohihi,  coeffimhilohi,  coeffimlolohi,
              coeffimhihilo,  coeffimlohilo,  coeffimhilolo,  coeffimlololo,
               workimhihihi[j],workimlohihi[j],workimhilohi[j],workimlolohi[j],
               workimhihilo[j],workimlohilo[j],workimhilolo[j],workimlololo[j],
                &acc2hihihi,    &acc2lohihi,    &acc2hilohi,    &acc2lolohi,
                &acc2hihilo,    &acc2lohilo,    &acc2hilolo,    &acc2lololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                   acc1hihihi,     acc1lohihi,     acc1hilohi,     acc1lolohi,
                   acc1hihilo,     acc1lohilo,     acc1hilolo,     acc1lololo);
      odg_dec(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                   acc2hihihi,     acc2lohihi,     acc2hilohi,     acc2lolohi,
                   acc2hihilo,     acc2lohilo,     acc2hilolo,     acc2lololo);
      odg_mul(coeffimhihihi,  coeffimlohihi,  coeffimhilohi,  coeffimlolohi,
              coeffimhihilo,  coeffimlohilo,  coeffimhilolo,  coeffimlololo,
               workrehihihi[j],workrelohihi[j],workrehilohi[j],workrelolohi[j],
               workrehihilo[j],workrelohilo[j],workrehilolo[j],workrelololo[j],
                &acc1hihihi,    &acc1lohihi,    &acc1hilohi,    &acc1lolohi,
                &acc1hihilo,    &acc1lohilo,    &acc1hilolo,    &acc1lololo);
      odg_mul(coeffrehihihi,  coeffrelohihi,  coeffrehilohi,  coeffrelolohi,
              coeffrehihilo,  coeffrelohilo,  coeffrehilolo,  coeffrelololo,
               workimhihihi[j],workimlohihi[j],workimhilohi[j],workimlolohi[j],
               workimhihilo[j],workimlohilo[j],workimhilolo[j],workimlololo[j],
                &acc2hihihi,    &acc2lohihi,    &acc2hilohi,    &acc2lolohi,
                &acc2hihilo,    &acc2lohilo,    &acc2hilolo,    &acc2lololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                   acc1hihihi,     acc1lohihi,     acc1hilohi,     acc1lolohi,
                   acc1hihilo,     acc1lohilo,     acc1hilolo,     acc1lololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                   acc2hihihi,     acc2lohihi,     acc2hilohi,    acc2lolohi,
                   acc2hihilo,     acc2lohilo,     acc2hilolo,    acc2lololo);
   }
   wrehihihi[rhsoff+k] = resultrehihihi;
   wrelohihi[rhsoff+k] = resultrelohihi;
   wrehilohi[rhsoff+k] = resultrehilohi;
   wrelolohi[rhsoff+k] = resultrelolohi;
   wrehihilo[rhsoff+k] = resultrehihilo;
   wrelohilo[rhsoff+k] = resultrelohilo;
   wrehilolo[rhsoff+k] = resultrehilolo;
   wrelololo[rhsoff+k] = resultrelololo;
   wimhihihi[rhsoff+k] = resultimhihihi;
   wimlohihi[rhsoff+k] = resultimlohihi;
   wimhilohi[rhsoff+k] = resultimhilohi;
   wimlolohi[rhsoff+k] = resultimlolohi;
   wimhihilo[rhsoff+k] = resultimhihilo;
   wimlohilo[rhsoff+k] = resultimlohilo;
   wimhilolo[rhsoff+k] = resultimhilolo;
   wimlololo[rhsoff+k] = resultimlololo;
}

__global__ void dbl8_back_substitute
 ( int dim, int idx, 
   double *Uhihihi, double *Ulohihi, double *Uhilohi, double *Ulolohi, 
   double *Uhihilo, double *Ulohilo, double *Uhilolo, double *Ulololo, 
   double *whihihi, double *wlohihi, double *whilohi, double *wlolohi,
   double *whihilo, double *wlohilo, double *whilolo, double *wlololo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkhihihi[tabsod_shmemsize]; // copy of w
   __shared__ double wrklohihi[tabsod_shmemsize]; 
   __shared__ double wrkhilohi[tabsod_shmemsize];
   __shared__ double wrklolohi[tabsod_shmemsize]; 
   __shared__ double wrkhihilo[tabsod_shmemsize];
   __shared__ double wrklohilo[tabsod_shmemsize]; 
   __shared__ double wrkhilolo[tabsod_shmemsize];
   __shared__ double wrklololo[tabsod_shmemsize]; 
   __shared__ double solhihihi[tabsod_shmemsize]; // solution to update with
   __shared__ double sollohihi[tabsod_shmemsize];
   __shared__ double solhilohi[tabsod_shmemsize];
   __shared__ double sollolohi[tabsod_shmemsize];
   __shared__ double solhihilo[tabsod_shmemsize];
   __shared__ double sollohilo[tabsod_shmemsize];
   __shared__ double solhilolo[tabsod_shmemsize];
   __shared__ double sollololo[tabsod_shmemsize];

   double resulthihihi = 0.0; // each thread stores its product in result
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;
   double coeffhihihi,coefflohihi,coeffhilohi,coefflolohi;
   double coeffhihilo,coefflohilo,coeffhilolo,coefflololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   wrkhihihi[k] = whihihi[B*dim+k];   // block B updates B-th slice of w
   wrklohihi[k] = wlohihi[B*dim+k];
   wrkhilohi[k] = whilohi[B*dim+k];
   wrklolohi[k] = wlolohi[B*dim+k];
   wrkhihilo[k] = whihilo[B*dim+k]; 
   wrklohilo[k] = wlohilo[B*dim+k];
   wrkhilolo[k] = whilolo[B*dim+k];
   wrklololo[k] = wlololo[B*dim+k];
   solhihihi[k] = whihihi[idx*dim+k]; // solution that is back substituted
   sollohihi[k] = wlohihi[idx*dim+k];
   solhilohi[k] = whilohi[idx*dim+k];
   sollolohi[k] = wlolohi[idx*dim+k];
   solhihilo[k] = whihilo[idx*dim+k];
   sollohilo[k] = wlohilo[idx*dim+k];
   solhilolo[k] = whilolo[idx*dim+k];
   sollololo[k] = wlololo[idx*dim+k];

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      __syncthreads();
      coeffhihihi = Uhihihi[offset+k*dim+j];
      coefflohihi = Ulohihi[offset+k*dim+j];
      coeffhilohi = Uhilohi[offset+k*dim+j];
      coefflolohi = Ulolohi[offset+k*dim+j];
      coeffhihilo = Uhihilo[offset+k*dim+j];
      coefflohilo = Ulohilo[offset+k*dim+j];
      coeffhilolo = Uhilolo[offset+k*dim+j];
      coefflololo = Ulololo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      __syncthreads();
      odg_mul(coeffhihihi, coefflohihi, coeffhilohi, coefflolohi,
              coeffhihilo, coefflohilo, coeffhilolo, coefflololo,
                solhihihi[j],sollohihi[j],solhilohi[j],sollolohi[j],
                solhihilo[j],sollohilo[j],solhilolo[j],sollololo[j],
               &acchihihi,  &acclohihi,  &acchilohi,  &acclolohi,
               &acchihilo,  &acclohilo,  &acchilolo,  &acclololo);
      odg_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
                  acchihihi,    acclohihi,    acchilohi,    acclolohi,
                  acchihilo,    acclohilo,    acchilolo,    acclololo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   __syncthreads();
   odg_dec(  &wrkhihihi[k],&wrklohihi[k],&wrkhilohi[k],&wrklolohi[k],
             &wrkhihilo[k],&wrklohilo[k],&wrkhilolo[k],&wrklololo[k],
           resulthihihi, resultlohihi, resulthilohi, resultlolohi,
           resulthihilo, resultlohilo, resulthilolo, resultlololo);
   __syncthreads();
   whihihi[B*dim+k] = wrkhihihi[k];
   wlohihi[B*dim+k] = wrklohihi[k];
   whilohi[B*dim+k] = wrkhilohi[k];
   wlolohi[B*dim+k] = wrklolohi[k];
   whihilo[B*dim+k] = wrkhihilo[k];
   wlohilo[B*dim+k] = wrklohilo[k];
   whilolo[B*dim+k] = wrkhilolo[k];
   wlololo[B*dim+k] = wrklololo[k];
}

__global__ void cmplx8_back_substitute
 ( int dim, int idx,
   double *Urehihihi, double *Urelohihi,
   double *Urehilohi, double *Urelolohi,
   double *Urehihilo, double *Urelohilo,
   double *Urehilolo, double *Urelololo,
   double *Uimhihihi, double *Uimlohihi,
   double *Uimhilohi, double *Uimlolohi,
   double *Uimhihilo, double *Uimlohilo,
   double *Uimhilolo, double *Uimlololo,
   double *wrehihihi, double *wrelohihi,
   double *wrehilohi, double *wrelolohi,
   double *wrehihilo, double *wrelohilo,
   double *wrehilolo, double *wrelololo,
   double *wimhihihi, double *wimlohihi,
   double *wimhilohi, double *wimlolohi,
   double *wimhihilo, double *wimlohilo,
   double *wimhilolo, double *wimlololo )
{
   const int B = blockIdx.x;     // block index
   const int k = threadIdx.x;    // thread k computes k-th product
   const int offset = B*dim*dim; // numbers to skip

   __shared__ double wrkrehihihi[tabsod_shmemsize]; // copy of w
   __shared__ double wrkrelohihi[tabsod_shmemsize]; 
   __shared__ double wrkrehilohi[tabsod_shmemsize];
   __shared__ double wrkrelolohi[tabsod_shmemsize]; 
   __shared__ double wrkrehihilo[tabsod_shmemsize];
   __shared__ double wrkrelohilo[tabsod_shmemsize]; 
   __shared__ double wrkrehilolo[tabsod_shmemsize];
   __shared__ double wrkrelololo[tabsod_shmemsize]; 
   __shared__ double wrkimhihihi[tabsod_shmemsize];
   __shared__ double wrkimlohihi[tabsod_shmemsize]; 
   __shared__ double wrkimhilohi[tabsod_shmemsize];
   __shared__ double wrkimlolohi[tabsod_shmemsize]; 
   __shared__ double wrkimhihilo[tabsod_shmemsize];
   __shared__ double wrkimlohilo[tabsod_shmemsize]; 
   __shared__ double wrkimhilolo[tabsod_shmemsize];
   __shared__ double wrkimlololo[tabsod_shmemsize]; 
   __shared__ double solrehihihi[tabsod_shmemsize]; // solution in update
   __shared__ double solrelohihi[tabsod_shmemsize];
   __shared__ double solrehilohi[tabsod_shmemsize];
   __shared__ double solrelolohi[tabsod_shmemsize];
   __shared__ double solrehihilo[tabsod_shmemsize];
   __shared__ double solrelohilo[tabsod_shmemsize];
   __shared__ double solrehilolo[tabsod_shmemsize];
   __shared__ double solrelololo[tabsod_shmemsize];
   __shared__ double solimhihihi[tabsod_shmemsize];
   __shared__ double solimlohihi[tabsod_shmemsize];
   __shared__ double solimhilohi[tabsod_shmemsize];
   __shared__ double solimlolohi[tabsod_shmemsize];
   __shared__ double solimhihilo[tabsod_shmemsize];
   __shared__ double solimlohilo[tabsod_shmemsize];
   __shared__ double solimhilolo[tabsod_shmemsize];
   __shared__ double solimlololo[tabsod_shmemsize];

   double resultrehihihi = 0.0; // each thread stores its product in result
   double resultrelohihi = 0.0;
   double resultrehilohi = 0.0;
   double resultrelolohi = 0.0;
   double resultrehihilo = 0.0;
   double resultrelohilo = 0.0;
   double resultrehilolo = 0.0;
   double resultrelololo = 0.0;
   double resultimhihihi = 0.0;
   double resultimlohihi = 0.0;
   double resultimhilohi = 0.0;
   double resultimlolohi = 0.0;
   double resultimhihilo = 0.0;
   double resultimlohilo = 0.0;
   double resultimhilolo = 0.0;
   double resultimlololo = 0.0;
   double coeffhihihi,coefflohihi,coeffhilohi,coefflolohi;
   double coeffhihilo,coefflohilo,coeffhilolo,coefflololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   wrkrehihihi[k] = wrehihihi[B*dim+k];  // block B updates B-th slice of w
   wrkrelohihi[k] = wrelohihi[B*dim+k];
   wrkrehilohi[k] = wrehilohi[B*dim+k];
   wrkrelolohi[k] = wrelolohi[B*dim+k];
   wrkrehihilo[k] = wrehihilo[B*dim+k];
   wrkrelohilo[k] = wrelohilo[B*dim+k];
   wrkrehilolo[k] = wrehilolo[B*dim+k];
   wrkrelololo[k] = wrelololo[B*dim+k];
   wrkimhihihi[k] = wimhihihi[B*dim+k];
   wrkimlohihi[k] = wimlohihi[B*dim+k];
   wrkimhilohi[k] = wimhilohi[B*dim+k];
   wrkimlolohi[k] = wimlolohi[B*dim+k];
   wrkimhihilo[k] = wimhihilo[B*dim+k];
   wrkimlohilo[k] = wimlohilo[B*dim+k];
   wrkimhilolo[k] = wimhilolo[B*dim+k];
   wrkimlololo[k] = wimlololo[B*dim+k];
   solrehihihi[k] = wrehihihi[idx*dim+k];  // solution back substituted
   solrelohihi[k] = wrelohihi[idx*dim+k];
   solrehilohi[k] = wrehilohi[idx*dim+k];
   solrelolohi[k] = wrelolohi[idx*dim+k];
   solrehihilo[k] = wrehihilo[idx*dim+k];
   solrelohilo[k] = wrelohilo[idx*dim+k];
   solrehilolo[k] = wrehilolo[idx*dim+k];
   solrelololo[k] = wrelololo[idx*dim+k];
   solimhihihi[k] = wimhihihi[idx*dim+k];
   solimlohihi[k] = wimlohihi[idx*dim+k];
   solimhilohi[k] = wimhilohi[idx*dim+k];
   solimlolohi[k] = wimlolohi[idx*dim+k];
   solimhihilo[k] = wimhihilo[idx*dim+k];
   solimlohilo[k] = wimlohilo[idx*dim+k];
   solimhilolo[k] = wimhilolo[idx*dim+k];
   solimlololo[k] = wimlololo[idx*dim+k];

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhihihi = Urehihihi[offset+k*dim+j];
      coefflohihi = Urelohihi[offset+k*dim+j];
      coeffhilohi = Urehilohi[offset+k*dim+j];
      coefflolohi = Urelolohi[offset+k*dim+j];
      coeffhihilo = Urehihilo[offset+k*dim+j];
      coefflohilo = Urelohilo[offset+k*dim+j];
      coeffhilolo = Urehilolo[offset+k*dim+j];
      coefflololo = Urelololo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      odg_mul(coeffhihihi, coefflohihi, coeffhilohi, coefflolohi,
              coeffhihilo, coefflohilo, coeffhilolo, coefflololo,
                solrehihihi[j],solrelohihi[j],solrehilohi[j],solrelolohi[j],
                solrehihilo[j],solrelohilo[j],solrehilolo[j],solrelololo[j],
                 &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                 &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      odg_mul(coeffhihihi, coefflohihi, coeffhilohi, coefflolohi,
              coeffhihilo, coefflohilo, coeffhilolo, coefflololo,
                solimhihihi[j],solimlohihi[j],solimhilohi[j],solimlolohi[j],
                solimhihilo[j],solimlohilo[j],solimhilolo[j],solimlololo[j],
                 &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                 &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      coeffhihihi = Uimhihihi[offset+k*dim+j];
      coefflohihi = Uimlohihi[offset+k*dim+j];
      coeffhilohi = Uimhilohi[offset+k*dim+j];
      coefflolohi = Uimlolohi[offset+k*dim+j];
      coeffhihilo = Uimhihilo[offset+k*dim+j];
      coefflohilo = Uimlohilo[offset+k*dim+j];
      coeffhilolo = Uimhilolo[offset+k*dim+j];
      coefflololo = Uimlololo[offset+k*dim+j];
      odg_mul(coeffhihihi, coefflohihi, coeffhilohi, coefflolohi,
              coeffhihilo, coefflohilo, coeffhilolo, coefflololo,
                solimhihihi[j],solimlohihi[j],solimhilohi[j],solimlolohi[j],
                solimhihilo[j],solimlohilo[j],solimhilolo[j],solimlololo[j],
                 &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                 &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_dec(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
      odg_mul(coeffhihihi, coefflohihi, coeffhilohi, coefflolohi,
              coeffhihilo, coefflohilo, coeffhilolo, coefflololo,
                solrehihihi[j],solrelohihi[j],solrehilohi[j],solrelolohi[j],
                solrehihilo[j],solrelohilo[j],solrehilolo[j],solrelololo[j],
                 &acchihihi,    &acclohihi,    &acchilohi,    &acclolohi,
                 &acchihilo,    &acclohilo,    &acchilolo,    &acclololo);
      odg_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
                    acchihihi,      acclohihi,      acchilohi,      acclolohi,
                    acchihilo,      acclohilo,      acchilolo,      acclololo);
   }
   // wrk[k] = wrk[k] - result; // subtract product
   odg_dec(  &wrkrehihihi[k],&wrkrelohihi[k],&wrkrehilohi[k],&wrkrelolohi[k],
             &wrkrehihilo[k],&wrkrelohilo[k],&wrkrehilolo[k],&wrkrelololo[k],
           resultrehihihi, resultrelohihi, resultrehilohi, resultrelolohi,
           resultrehihilo, resultrelohilo, resultrehilolo, resultrelololo);
   odg_dec(  &wrkimhihihi[k],&wrkimlohihi[k],&wrkimhilohi[k],&wrkimlolohi[k],
             &wrkimhihilo[k],&wrkimlohilo[k],&wrkimhilolo[k],&wrkimlololo[k],
           resultimhihihi, resultimlohihi, resultimhilohi, resultimlolohi,
           resultimhihilo, resultimlohilo, resultimhilolo, resultimlololo);
   wrehihihi[B*dim+k] = wrkrehihihi[k];
   wrelohihi[B*dim+k] = wrkrelohihi[k];
   wrehilohi[B*dim+k] = wrkrehilohi[k];
   wrelolohi[B*dim+k] = wrkrelolohi[k];
   wrehihilo[B*dim+k] = wrkrehihilo[k];
   wrelohilo[B*dim+k] = wrkrelohilo[k];
   wrehilolo[B*dim+k] = wrkrehilolo[k];
   wrelololo[B*dim+k] = wrkrelololo[k];
   wimhihihi[B*dim+k] = wrkimhihihi[k];
   wimlohihi[B*dim+k] = wrkimlohihi[k];
   wimhilohi[B*dim+k] = wrkimhilohi[k];
   wimlolohi[B*dim+k] = wrkimlolohi[k];
   wimhihilo[B*dim+k] = wrkimhihilo[k];
   wimlohilo[B*dim+k] = wrkimlohilo[k];
   wimhilolo[B*dim+k] = wrkimhilolo[k];
   wimlololo[B*dim+k] = wrkimlololo[k];
}

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

   if(dim <= 13)
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

void GPU_cmplx8_upper_inverse
 ( int dim,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double **invUrehihihi, double **invUrelohihi,
   double **invUrehilohi, double **invUrelolohi,
   double **invUrehihilo, double **invUrelohilo,
   double **invUrehilolo, double **invUrelololo,
   double **invUimhihihi, double **invUimlohihi, 
   double **invUimhilohi, double **invUimlolohi, 
   double **invUimhihilo, double **invUimlohilo, 
   double **invUimhilolo, double **invUimlololo, 
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Urehihihi_h = new double[szU];  // real parts of U
   double *Urelohihi_h = new double[szU]; 
   double *Urehilohi_h = new double[szU]; 
   double *Urelolohi_h = new double[szU]; 
   double *Urehihilo_h = new double[szU];
   double *Urelohilo_h = new double[szU]; 
   double *Urehilolo_h = new double[szU]; 
   double *Urelololo_h = new double[szU]; 
   double *Uimhihihi_h = new double[szU];  // imaginary parts of U
   double *Uimlohihi_h = new double[szU]; 
   double *Uimhilohi_h = new double[szU]; 
   double *Uimlolohi_h = new double[szU]; 
   double *Uimhihilo_h = new double[szU];
   double *Uimlohilo_h = new double[szU]; 
   double *Uimhilolo_h = new double[szU]; 
   double *Uimlololo_h = new double[szU]; 
   double *Urehihihi_d;                    // real parts on the device
   double *Urelohihi_d;
   double *Urehilohi_d;
   double *Urelolohi_d;
   double *Urehihilo_d;
   double *Urelohilo_d;
   double *Urehilolo_d;
   double *Urelololo_d;
   double *Uimhihihi_d;                    // imaginary parts on the device
   double *Uimlohihi_d;
   double *Uimhilohi_d;
   double *Uimlolohi_d;
   double *Uimhihilo_d;
   double *Uimlohilo_d;
   double *Uimhilolo_d;
   double *Uimlololo_d;
   double *invUrehihihi_h = new double[szU]; // real parts of inverse
   double *invUrelohihi_h = new double[szU];
   double *invUrehilohi_h = new double[szU];
   double *invUrelolohi_h = new double[szU];
   double *invUrehihilo_h = new double[szU];
   double *invUrelohilo_h = new double[szU];
   double *invUrehilolo_h = new double[szU];
   double *invUrelololo_h = new double[szU];
   double *invUimhihihi_h = new double[szU]; // imaginary parts of inverse
   double *invUimlohihi_h = new double[szU];
   double *invUimhilohi_h = new double[szU];
   double *invUimlolohi_h = new double[szU];
   double *invUimhihilo_h = new double[szU];
   double *invUimlohilo_h = new double[szU];
   double *invUimhilolo_h = new double[szU];
   double *invUimlololo_h = new double[szU];
   double *invUrehihihi_d;           // inverse on device
   double *invUrelohihi_d;
   double *invUrehilohi_d;
   double *invUrelolohi_d;
   double *invUrehihilo_d;
   double *invUrelohilo_d;
   double *invUrehilolo_d;
   double *invUrelololo_d;
   double *invUimhihihi_d; 
   double *invUimlohihi_d;
   double *invUimhilohi_d;
   double *invUimlolohi_d;
   double *invUimhihilo_d; 
   double *invUimlohilo_d;
   double *invUimhilolo_d;
   double *invUimlololo_d;

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Urehihihi_h[ix]   = Urehihihi[i][j];
         Urelohihi_h[ix]   = Urelohihi[i][j];
         Urehilohi_h[ix]   = Urehilohi[i][j];
         Urelolohi_h[ix]   = Urelolohi[i][j];
         Urehihilo_h[ix]   = Urehihilo[i][j];
         Urelohilo_h[ix]   = Urelohilo[i][j];
         Urehilolo_h[ix]   = Urehilolo[i][j];
         Urelololo_h[ix]   = Urelololo[i][j];
         Uimhihihi_h[ix]   = Uimhihihi[i][j];
         Uimlohihi_h[ix]   = Uimlohihi[i][j];
         Uimhilohi_h[ix]   = Uimhilohi[i][j];
         Uimlolohi_h[ix]   = Uimlolohi[i][j];
         Uimhihilo_h[ix]   = Uimhihilo[i][j];
         Uimlohilo_h[ix]   = Uimlohilo[i][j];
         Uimhilolo_h[ix]   = Uimhilolo[i][j];
         Uimlololo_h[ix++] = Uimlololo[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Urehihihi_d,szmat);
   cudaMalloc((void**)&Urelohihi_d,szmat);
   cudaMalloc((void**)&Urehilohi_d,szmat);
   cudaMalloc((void**)&Urelolohi_d,szmat);
   cudaMalloc((void**)&Urehihilo_d,szmat);
   cudaMalloc((void**)&Urelohilo_d,szmat);
   cudaMalloc((void**)&Urehilolo_d,szmat);
   cudaMalloc((void**)&Urelololo_d,szmat);
   cudaMalloc((void**)&Uimhihihi_d,szmat);
   cudaMalloc((void**)&Uimlohihi_d,szmat);
   cudaMalloc((void**)&Uimhilohi_d,szmat);
   cudaMalloc((void**)&Uimlolohi_d,szmat);
   cudaMalloc((void**)&Uimhihilo_d,szmat);
   cudaMalloc((void**)&Uimlohilo_d,szmat);
   cudaMalloc((void**)&Uimhilolo_d,szmat);
   cudaMalloc((void**)&Uimlololo_d,szmat);
   cudaMalloc((void**)&invUrehihihi_d,szmat);
   cudaMalloc((void**)&invUrelohihi_d,szmat);
   cudaMalloc((void**)&invUrehilohi_d,szmat);
   cudaMalloc((void**)&invUrelolohi_d,szmat);
   cudaMalloc((void**)&invUrehihilo_d,szmat);
   cudaMalloc((void**)&invUrelohilo_d,szmat);
   cudaMalloc((void**)&invUrehilolo_d,szmat);
   cudaMalloc((void**)&invUrelololo_d,szmat);
   cudaMalloc((void**)&invUimhihihi_d,szmat);
   cudaMalloc((void**)&invUimlohihi_d,szmat);
   cudaMalloc((void**)&invUimhilohi_d,szmat);
   cudaMalloc((void**)&invUimlolohi_d,szmat);
   cudaMalloc((void**)&invUimhihilo_d,szmat);
   cudaMalloc((void**)&invUimlohilo_d,szmat);
   cudaMalloc((void**)&invUimhilolo_d,szmat);
   cudaMalloc((void**)&invUimlololo_d,szmat);
   cudaMemcpy(Urehihihi_d,Urehihihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelohihi_d,Urelohihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urehilohi_d,Urehilohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelolohi_d,Urelolohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urehihilo_d,Urehihilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelohilo_d,Urelohilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urehilolo_d,Urehilolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Urelololo_d,Urelololo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhihihi_d,Uimhihihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlohihi_d,Uimlohihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhilohi_d,Uimhilohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlolohi_d,Uimlolohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhihilo_d,Uimhihilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlohilo_d,Uimlohilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimhilolo_d,Uimhilolo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uimlololo_d,Uimlololo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 13)
      cmplx8_small_invert_upper<<<1,dim>>>
         (dim,   Urehihihi_d,   Urelohihi_d,   Urehilohi_d,   Urelolohi_d,
                 Urehihilo_d,   Urelohilo_d,   Urehilolo_d,   Urelololo_d,
                 Uimhihihi_d,   Uimlohihi_d,   Uimhilohi_d,   Uimlolohi_d,
                 Uimhihilo_d,   Uimlohilo_d,   Uimhilolo_d,   Uimlololo_d,
              invUrehihihi_d,invUrelohihi_d,invUrehilohi_d,invUrelolohi_d,
              invUrehihilo_d,invUrelohilo_d,invUrehilolo_d,invUrelololo_d,
              invUimhihihi_d,invUimlohihi_d,invUimhilohi_d,invUimlolohi_d,
              invUimhihilo_d,invUimlohilo_d,invUimhilolo_d,invUimlololo_d);
   else
      cmplx8_medium_invert_upper<<<1,dim>>>
         (dim,   Urehihihi_d,   Urelohihi_d,   Urehilohi_d,   Urelolohi_d,
                 Urehihilo_d,   Urelohilo_d,   Urehilolo_d,   Urelololo_d,
                 Uimhihihi_d,   Uimlohihi_d,   Uimhilohi_d,   Uimlolohi_d,
                 Uimhihilo_d,   Uimlohilo_d,   Uimhilolo_d,   Uimlololo_d,
              invUrehihihi_d,invUrelohihi_d,invUrehilohi_d,invUrelolohi_d,
              invUrehihilo_d,invUrelohilo_d,invUrehilolo_d,invUrelololo_d,
              invUimhihihi_d,invUimlohihi_d,invUimhilohi_d,invUimlolohi_d,
              invUimhihilo_d,invUimlohilo_d,invUimhilolo_d,invUimlololo_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUrehihihi_h,invUrehihihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelohihi_h,invUrelohihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrehilohi_h,invUrehilohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelolohi_h,invUrelolohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrehihilo_h,invUrehihilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelohilo_h,invUrelohilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrehilolo_h,invUrehilolo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUrelololo_h,invUrelololo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhihihi_h,invUimhihihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlohihi_h,invUimlohihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhilohi_h,invUimhilohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlolohi_h,invUimlolohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhihilo_h,invUimhihilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlohilo_h,invUimlohilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimhilolo_h,invUimhilolo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUimlololo_h,invUimlololo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUrehihihi[i][j] = invUrehihihi_h[ix];
         invUrelohihi[i][j] = invUrelohihi_h[ix];
         invUrehilohi[i][j] = invUrehilohi_h[ix];
         invUrelolohi[i][j] = invUrelolohi_h[ix];
         invUrehihilo[i][j] = invUrehihilo_h[ix];
         invUrelohilo[i][j] = invUrelohilo_h[ix];
         invUrehilolo[i][j] = invUrehilolo_h[ix];
         invUrelololo[i][j] = invUrelololo_h[ix];
         invUimhihihi[i][j] = invUimhihihi_h[ix];
         invUimlohihi[i][j] = invUimlohihi_h[ix];
         invUimhilohi[i][j] = invUimhilohi_h[ix];
         invUimlolohi[i][j] = invUimlolohi_h[ix];
         invUimhihilo[i][j] = invUimhihilo_h[ix];
         invUimlohilo[i][j] = invUimlohilo_h[ix];
         invUimhilolo[i][j] = invUimhilolo_h[ix];
         invUimlololo[i][j] = invUimlololo_h[ix++];
      }

   free(Urehihihi_h); free(Urelohihi_h);
   free(Urehilohi_h); free(Urelolohi_h);
   free(Urehihilo_h); free(Urelohilo_h);
   free(Urehilolo_h); free(Urelololo_h);
   free(Uimhihihi_h); free(Uimlohihi_h);
   free(Uimhilohi_h); free(Uimlolohi_h);
   free(Uimhihilo_h); free(Uimlohilo_h);
   free(Uimhilolo_h); free(Uimlololo_h);
   free(invUrehihihi_h); free(invUrelohihi_h);
   free(invUrehilohi_h); free(invUrelolohi_h);
   free(invUrehihilo_h); free(invUrelohilo_h);
   free(invUrehilolo_h); free(invUrelololo_h);
   free(invUimhihihi_h); free(invUimlohihi_h);
   free(invUimhilohi_h); free(invUimlolohi_h);
   free(invUimhihilo_h); free(invUimlohilo_h);
   free(invUimhilolo_h); free(invUimlololo_h);
}

void GPU_dbl8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;   // number of doubles on diagonal tiles
   double *Dhihihi_h = new double[nbr]; // diagonal tiles on the host
   double *Dlohihi_h = new double[nbr]; 
   double *Dhilohi_h = new double[nbr]; 
   double *Dlolohi_h = new double[nbr];
   double *Dhihilo_h = new double[nbr];
   double *Dlohilo_h = new double[nbr]; 
   double *Dhilolo_h = new double[nbr]; 
   double *Dlololo_h = new double[nbr];
   double *Dhihihi_d;                   // diagonal tiles on the device
   double *Dlohihi_d; 
   double *Dhilohi_d;
   double *Dlolohi_d; 
   double *Dhihilo_d;
   double *Dlohilo_d; 
   double *Dhilolo_d;
   double *Dlololo_d; 
   double *invDhihihi_h = new double[nbr]; // inverse of tiles on host 
   double *invDlohihi_h = new double[nbr]; 
   double *invDhilohi_h = new double[nbr]; 
   double *invDlolohi_h = new double[nbr];
   double *invDhihilo_h = new double[nbr];
   double *invDlohilo_h = new double[nbr]; 
   double *invDhilolo_h = new double[nbr]; 
   double *invDlololo_h = new double[nbr];
   double *invDhihihi_d;            // inverse diagonal tiles on device
   double *invDlohihi_d;
   double *invDhilohi_d;
   double *invDlolohi_d;
   double *invDhihilo_d;
   double *invDlohilo_d;
   double *invDhilolo_d;
   double *invDlololo_d;
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Dhihihi_h[ix]   = Uhihihi[offset+i][offset+j];
            Dlohihi_h[ix]   = Ulohihi[offset+i][offset+j];
            Dhilohi_h[ix]   = Uhilohi[offset+i][offset+j];
            Dlolohi_h[ix]   = Ulolohi[offset+i][offset+j];
            Dhihilo_h[ix]   = Uhihilo[offset+i][offset+j];
            Dlohilo_h[ix]   = Ulohilo[offset+i][offset+j];
            Dhilolo_h[ix]   = Uhilolo[offset+i][offset+j];
            Dlololo_h[ix++] = Ulololo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Dhihihi_d,sznum);
   cudaMalloc((void**)&Dlohihi_d,sznum);
   cudaMalloc((void**)&Dhilohi_d,sznum);
   cudaMalloc((void**)&Dlolohi_d,sznum);
   cudaMalloc((void**)&Dhihilo_d,sznum);
   cudaMalloc((void**)&Dlohilo_d,sznum);
   cudaMalloc((void**)&Dhilolo_d,sznum);
   cudaMalloc((void**)&Dlololo_d,sznum);
   cudaMalloc((void**)&invDhihihi_d,sznum);
   cudaMalloc((void**)&invDlohihi_d,sznum);
   cudaMalloc((void**)&invDhilohi_d,sznum);
   cudaMalloc((void**)&invDlolohi_d,sznum);
   cudaMalloc((void**)&invDhihilo_d,sznum);
   cudaMalloc((void**)&invDlohilo_d,sznum);
   cudaMalloc((void**)&invDhilolo_d,sznum);
   cudaMalloc((void**)&invDlololo_d,sznum);
   cudaMemcpy(Dhihihi_d,Dhihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlohihi_d,Dlohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dhilohi_d,Dhilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlolohi_d,Dlolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dhihilo_d,Dhihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlohilo_d,Dlohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dhilolo_d,Dhilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dlololo_d,Dlololo_h,sznum,cudaMemcpyHostToDevice);

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
   dbl8_invert_tiles<<<nbt,szt>>>
      (szt,Dhihihi_d,   Dlohihi_d,   Dhilohi_d,   Dlolohi_d,
           Dhihilo_d,   Dlohilo_d,   Dhilolo_d,   Dlololo_d,
        invDhihihi_d,invDlohihi_d,invDhilohi_d,invDlolohi_d,
        invDhihilo_d,invDlohilo_d,invDhilolo_d,invDlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhshihihi_d;                    // right hand side on device
   double *rhslohihi_d;
   double *rhshilohi_d;
   double *rhslolohi_d;
   double *rhshihilo_d;
   double *rhslohilo_d;
   double *rhshilolo_d;
   double *rhslololo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhshihihi_d,szrhs);
   cudaMalloc((void**)&rhslohihi_d,szrhs);
   cudaMalloc((void**)&rhshilohi_d,szrhs);
   cudaMalloc((void**)&rhslolohi_d,szrhs);
   cudaMalloc((void**)&rhshihilo_d,szrhs);
   cudaMalloc((void**)&rhslohilo_d,szrhs);
   cudaMalloc((void**)&rhshilolo_d,szrhs);
   cudaMalloc((void**)&rhslololo_d,szrhs);
   cudaMemcpy(rhshihihi_d,bhihihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslohihi_d,blohihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhshilohi_d,bhilohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslolohi_d,blolohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhshihilo_d,bhihilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslohilo_d,blohilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhshilolo_d,bhilolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhslololo_d,blololo,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   dbl8_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDhihihi_d,invDlohihi_d,invDhilohi_d,invDlolohi_d,
                 invDhihilo_d,invDlohilo_d,invDhilolo_d,invDlololo_d,
                  rhshihihi_d, rhslohihi_d, rhshilohi_d, rhslolohi_d,
                  rhshihilo_d, rhslohilo_d, rhshilolo_d, rhslololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;             // #doubles in column of U
   double *Ucolhihihi_h = new double[nbrUcol];  // column of U on host
   double *Ucollohihi_h = new double[nbrUcol]; 
   double *Ucolhilohi_h = new double[nbrUcol];
   double *Ucollolohi_h = new double[nbrUcol];
   double *Ucolhihilo_h = new double[nbrUcol];
   double *Ucollohilo_h = new double[nbrUcol]; 
   double *Ucolhilolo_h = new double[nbrUcol];
   double *Ucollololo_h = new double[nbrUcol];
   double *Ucolhihihi_d;
   double *Ucollohihi_d;
   double *Ucolhilohi_d;
   double *Ucollolohi_d;
   double *Ucolhihilo_d;
   double *Ucollohilo_d;
   double *Ucolhilolo_d;
   double *Ucollololo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolhihihi_d,szUcol);
   cudaMalloc((void**)&Ucollohihi_d,szUcol);
   cudaMalloc((void**)&Ucolhilohi_d,szUcol);
   cudaMalloc((void**)&Ucollolohi_d,szUcol);
   cudaMalloc((void**)&Ucolhihilo_d,szUcol);
   cudaMalloc((void**)&Ucollohilo_d,szUcol);
   cudaMalloc((void**)&Ucolhilolo_d,szUcol);
   cudaMalloc((void**)&Ucollololo_d,szUcol);

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
               Ucolhihihi_h[ix]   = Uhihihi[rowoff+i][coloff+j];
               Ucollohihi_h[ix]   = Ulohihi[rowoff+i][coloff+j];
               Ucolhilohi_h[ix]   = Uhilohi[rowoff+i][coloff+j];
               Ucollolohi_h[ix]   = Ulolohi[rowoff+i][coloff+j];
               Ucolhihilo_h[ix]   = Uhihilo[rowoff+i][coloff+j];
               Ucollohilo_h[ix]   = Ulohilo[rowoff+i][coloff+j];
               Ucolhilolo_h[ix]   = Uhilolo[rowoff+i][coloff+j];
               Ucollololo_h[ix++] = Ulololo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolhihihi_d,Ucolhihihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollohihi_d,Ucollohihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolhilohi_d,Ucolhilohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollolohi_d,Ucollolohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolhihilo_d,Ucolhihilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollohilo_d,Ucollohilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolhilolo_d,Ucolhilolo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucollololo_d,Ucollololo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      dbl8_back_substitute<<<k,szt>>>
         (szt,k,Ucolhihihi_d,Ucollohihi_d,Ucolhilohi_d,Ucollolohi_d,
                Ucolhihilo_d,Ucollohilo_d,Ucolhilolo_d,Ucollololo_d,
                 rhshihihi_d, rhslohihi_d, rhshilohi_d, rhslolohi_d,
                 rhshihilo_d, rhslohilo_d, rhshilolo_d, rhslololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_dbl_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      dbl8_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDhihihi_d,invDlohihi_d,invDhilohi_d,invDlolohi_d,
                  invDhihilo_d,invDlohilo_d,invDhilolo_d,invDlololo_d,
                   rhshihihi_d, rhslohihi_d, rhshilohi_d, rhslolohi_d,
                   rhshihilo_d, rhslohilo_d, rhshilolo_d, rhslololo_d);
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

   cudaMemcpy(xhihihi,rhshihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlohihi,rhslohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xhilohi,rhshilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlolohi,rhslolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xhihilo,rhshihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlohilo,rhslohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xhilolo,rhshilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xlololo,rhslololo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDhihihi_h,invDhihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlohihi_h,invDlohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDhilohi_h,invDhilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlolohi_h,invDlolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDhihilo_h,invDhihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlohilo_h,invDlohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDhilolo_h,invDhilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDlololo_h,invDlololo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhihihi[offset+i][offset+j] = invDhihihi_h[ix];
            Ulohihi[offset+i][offset+j] = invDlohihi_h[ix];
            Uhilohi[offset+i][offset+j] = invDhilohi_h[ix];
            Ulolohi[offset+i][offset+j] = invDlolohi_h[ix];
            Uhihilo[offset+i][offset+j] = invDhihilo_h[ix];
            Ulohilo[offset+i][offset+j] = invDlohilo_h[ix];
            Uhilolo[offset+i][offset+j] = invDhilolo_h[ix];
            Ulololo[offset+i][offset+j] = invDlololo_h[ix++];
         }
   }
   free(Dhihihi_h); free(Dlohihi_h); free(Dhilohi_h); free(Dlolohi_h);
   free(Dhihilo_h); free(Dlohilo_h); free(Dhilolo_h); free(Dlololo_h);
   free(invDhihihi_h); free(invDlohihi_h);
   free(invDhilohi_h); free(invDlolohi_h);
   free(invDhihilo_h); free(invDlohilo_h);
   free(invDhilolo_h); free(invDlololo_h);
   free(Ucolhihihi_h); free(Ucollohihi_h);
   free(Ucolhilohi_h); free(Ucollolohi_h);
   free(Ucolhihilo_h); free(Ucollohilo_h);
   free(Ucolhilolo_h); free(Ucollololo_h);
}

void GPU_cmplx8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihihi, double **Urelohihi,
   double **Urehilohi, double **Urelolohi,
   double **Urehihilo, double **Urelohilo,
   double **Urehilolo, double **Urelololo,
   double **Uimhihihi, double **Uimlohihi,
   double **Uimhilohi, double **Uimlolohi,
   double **Uimhihilo, double **Uimlohilo,
   double **Uimhilolo, double **Uimlololo,
   double *brehihihi, double *brelohihi, double *brehilohi, double *brelolohi,
   double *brehihilo, double *brelohilo, double *brehilolo, double *brelololo,
   double *bimhihihi, double *bimlohihi, double *bimhilohi, double *bimlolohi,
   double *bimhihilo, double *bimlohilo, double *bimhilolo, double *bimlololo,
   double *xrehihihi, double *xrelohihi, double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo, double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi, double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo,
   double *invlapms, double *mullapms, double *sublapms, double *totlapms,
   double *walltimesec,
   long long int *addcnt, long long int *mulcnt, long long int *divcnt )
{
   const int nbr = nbt*szt*szt;       // number of doubles on diagonal tiles
   double *Drehihihi_h = new double[nbr]; // real parts of diagonal tiles
   double *Drelohihi_h = new double[nbr];  
   double *Drehilohi_h = new double[nbr];  
   double *Drelolohi_h = new double[nbr]; 
   double *Drehihilo_h = new double[nbr]; 
   double *Drelohilo_h = new double[nbr];  
   double *Drehilolo_h = new double[nbr];  
   double *Drelololo_h = new double[nbr]; 
   double *Dimhihihi_h = new double[nbr]; // imag parts of diagonal tiles
   double *Dimlohihi_h = new double[nbr];  
   double *Dimhilohi_h = new double[nbr];  
   double *Dimlolohi_h = new double[nbr];  
   double *Dimhihilo_h = new double[nbr];
   double *Dimlohilo_h = new double[nbr];  
   double *Dimhilolo_h = new double[nbr];  
   double *Dimlololo_h = new double[nbr];  
   double *Drehihihi_d;                    // diagonal tiles on the device
   double *Drelohihi_d;
   double *Drehilohi_d;
   double *Drelolohi_d;
   double *Drehihilo_d;
   double *Drelohilo_d;
   double *Drehilolo_d;
   double *Drelololo_d;
   double *Dimhihihi_d;
   double *Dimlohihi_d;
   double *Dimhilohi_d; 
   double *Dimlolohi_d;
   double *Dimhihilo_d;
   double *Dimlohilo_d;
   double *Dimhilolo_d; 
   double *Dimlololo_d;
   double *invDrehihihi_h = new double[nbr]; // real parts of inverse tiles
   double *invDrelohihi_h = new double[nbr];
   double *invDrehilohi_h = new double[nbr]; 
   double *invDrelolohi_h = new double[nbr];
   double *invDrehihilo_h = new double[nbr];
   double *invDrelohilo_h = new double[nbr];
   double *invDrehilolo_h = new double[nbr]; 
   double *invDrelololo_h = new double[nbr];
   double *invDimhihihi_h = new double[nbr]; // imag parts of inverse tiles
   double *invDimlohihi_h = new double[nbr];
   double *invDimhilohi_h = new double[nbr];
   double *invDimlolohi_h = new double[nbr];
   double *invDimhihilo_h = new double[nbr];
   double *invDimlohilo_h = new double[nbr];
   double *invDimhilolo_h = new double[nbr];
   double *invDimlololo_h = new double[nbr];
   double *invDrehihihi_d;                   // inverse tiles on device
   double *invDrelohihi_d; 
   double *invDrehilohi_d;
   double *invDrelolohi_d;
   double *invDrehihilo_d;
   double *invDrelohilo_d; 
   double *invDrehilolo_d;
   double *invDrelololo_d;
   double *invDimhihihi_d; 
   double *invDimlohihi_d;
   double *invDimhilohi_d;
   double *invDimlolohi_d;
   double *invDimhihilo_d; 
   double *invDimlohilo_d;
   double *invDimhilolo_d;
   double *invDimlololo_d;
   int offset;
   int ix = 0;

   for(int k=0; k<nbt; k++) // copy columns of the k-th tile
   {
      offset = k*szt;
      for(int j=0; j<szt; j++)
         for(int i=0; i<szt; i++)
         {
            Drehihihi_h[ix]   = Urehihihi[offset+i][offset+j];
            Drelohihi_h[ix]   = Urelohihi[offset+i][offset+j];
            Drehilohi_h[ix]   = Urehilohi[offset+i][offset+j];
            Drelolohi_h[ix]   = Urelolohi[offset+i][offset+j];
            Drehihilo_h[ix]   = Urehihilo[offset+i][offset+j];
            Drelohilo_h[ix]   = Urelohilo[offset+i][offset+j];
            Drehilolo_h[ix]   = Urehilolo[offset+i][offset+j];
            Drelololo_h[ix]   = Urelololo[offset+i][offset+j];
            Dimhihihi_h[ix]   = Uimhihihi[offset+i][offset+j];
            Dimlohihi_h[ix]   = Uimlohihi[offset+i][offset+j];
            Dimhilohi_h[ix]   = Uimhilohi[offset+i][offset+j];
            Dimlolohi_h[ix]   = Uimlolohi[offset+i][offset+j];
            Dimhihilo_h[ix]   = Uimhihilo[offset+i][offset+j];
            Dimlohilo_h[ix]   = Uimlohilo[offset+i][offset+j];
            Dimhilolo_h[ix]   = Uimhilolo[offset+i][offset+j];
            Dimlololo_h[ix++] = Uimlololo[offset+i][offset+j];
         }
   }
   const size_t sznum = nbr*sizeof(double);
   cudaMalloc((void**)&Drehihihi_d,sznum);
   cudaMalloc((void**)&Drelohihi_d,sznum);
   cudaMalloc((void**)&Drehilohi_d,sznum);
   cudaMalloc((void**)&Drelolohi_d,sznum);
   cudaMalloc((void**)&Drehihilo_d,sznum);
   cudaMalloc((void**)&Drelohilo_d,sznum);
   cudaMalloc((void**)&Drehilolo_d,sznum);
   cudaMalloc((void**)&Drelololo_d,sznum);
   cudaMalloc((void**)&Dimhihihi_d,sznum);
   cudaMalloc((void**)&Dimlohihi_d,sznum);
   cudaMalloc((void**)&Dimhilohi_d,sznum);
   cudaMalloc((void**)&Dimlolohi_d,sznum);
   cudaMalloc((void**)&Dimhihilo_d,sznum);
   cudaMalloc((void**)&Dimlohilo_d,sznum);
   cudaMalloc((void**)&Dimhilolo_d,sznum);
   cudaMalloc((void**)&Dimlololo_d,sznum);
   cudaMalloc((void**)&invDrehihihi_d,sznum);
   cudaMalloc((void**)&invDrelohihi_d,sznum);
   cudaMalloc((void**)&invDrehilohi_d,sznum);
   cudaMalloc((void**)&invDrelolohi_d,sznum);
   cudaMalloc((void**)&invDrehihilo_d,sznum);
   cudaMalloc((void**)&invDrelohilo_d,sznum);
   cudaMalloc((void**)&invDrehilolo_d,sznum);
   cudaMalloc((void**)&invDrelololo_d,sznum);
   cudaMalloc((void**)&invDimhihihi_d,sznum);
   cudaMalloc((void**)&invDimlohihi_d,sznum);
   cudaMalloc((void**)&invDimhilohi_d,sznum);
   cudaMalloc((void**)&invDimlolohi_d,sznum);
   cudaMalloc((void**)&invDimhihilo_d,sznum);
   cudaMalloc((void**)&invDimlohilo_d,sznum);
   cudaMalloc((void**)&invDimhilolo_d,sznum);
   cudaMalloc((void**)&invDimlololo_d,sznum);
   cudaMemcpy(Drehihihi_d,Drehihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelohihi_d,Drelohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drehilohi_d,Drehilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelolohi_d,Drelolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drehihilo_d,Drehihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelohilo_d,Drelohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drehilolo_d,Drehilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Drelololo_d,Drelololo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhihihi_d,Dimhihihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlohihi_d,Dimlohihi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhilohi_d,Dimhilohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlolohi_d,Dimlolohi_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhihilo_d,Dimhihilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlohilo_d,Dimlohilo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimhilolo_d,Dimhilolo_h,sznum,cudaMemcpyHostToDevice);
   cudaMemcpy(Dimlololo_d,Dimlololo_h,sznum,cudaMemcpyHostToDevice);

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
   cmplx8_invert_tiles<<<nbt,szt>>>
      (szt, Drehihihi_d,   Drelohihi_d,   Drehilohi_d,   Drelolohi_d,
            Drehihilo_d,   Drelohilo_d,   Drehilolo_d,   Drelololo_d,
            Dimhihihi_d,   Dimlohihi_d,   Dimhilohi_d,   Dimlolohi_d,
            Dimhihilo_d,   Dimlohilo_d,   Dimhilolo_d,   Dimlololo_d,
         invDrehihihi_d,invDrelohihi_d,invDrehilohi_d,invDrelolohi_d,
         invDrehihilo_d,invDrelohilo_d,invDrehilolo_d,invDrelololo_d,
         invDimhihihi_d,invDimlohihi_d,invDimhilohi_d,invDimlolohi_d,
         invDimhihilo_d,invDimlohilo_d,invDimhilolo_d,invDimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *invlapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);

   double *rhsrehihihi_d;                  // right hand side on device
   double *rhsrelohihi_d;
   double *rhsrehilohi_d;
   double *rhsrelolohi_d;
   double *rhsrehihilo_d; 
   double *rhsrelohilo_d;
   double *rhsrehilolo_d;
   double *rhsrelololo_d;
   double *rhsimhihihi_d;
   double *rhsimlohihi_d;
   double *rhsimhilohi_d;
   double *rhsimlolohi_d;
   double *rhsimhihilo_d;
   double *rhsimlohilo_d;
   double *rhsimhilolo_d;
   double *rhsimlololo_d;
   const size_t szrhs = dim*sizeof(double);
   cudaMalloc((void**)&rhsrehihihi_d,szrhs);
   cudaMalloc((void**)&rhsrelohihi_d,szrhs);
   cudaMalloc((void**)&rhsrehilohi_d,szrhs);
   cudaMalloc((void**)&rhsrelolohi_d,szrhs);
   cudaMalloc((void**)&rhsrehihilo_d,szrhs);
   cudaMalloc((void**)&rhsrelohilo_d,szrhs);
   cudaMalloc((void**)&rhsrehilolo_d,szrhs);
   cudaMalloc((void**)&rhsrelololo_d,szrhs);
   cudaMalloc((void**)&rhsimhihihi_d,szrhs);
   cudaMalloc((void**)&rhsimlohihi_d,szrhs);
   cudaMalloc((void**)&rhsimhilohi_d,szrhs);
   cudaMalloc((void**)&rhsimlolohi_d,szrhs);
   cudaMalloc((void**)&rhsimhihilo_d,szrhs);
   cudaMalloc((void**)&rhsimlohilo_d,szrhs);
   cudaMalloc((void**)&rhsimhilolo_d,szrhs);
   cudaMalloc((void**)&rhsimlololo_d,szrhs);
   cudaMemcpy(rhsrehihihi_d,brehihihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelohihi_d,brelohihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrehilohi_d,brehilohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelolohi_d,brelolohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrehihilo_d,brehihilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelohilo_d,brelohilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrehilolo_d,brehilolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsrelololo_d,brelololo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhihihi_d,bimhihihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlohihi_d,bimlohihi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhilohi_d,bimhilohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlolohi_d,bimlolohi,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhihilo_d,bimhihilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlohilo_d,bimlohilo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimhilolo_d,bimhilolo,szrhs,cudaMemcpyHostToDevice);
   cudaMemcpy(rhsimlololo_d,bimlololo,szrhs,cudaMemcpyHostToDevice);

   cudaEventRecord(start);
   cmplx8_multiply_inverse<<<1,szt>>>
      (szt,nbt-1,invDrehihihi_d,invDrelohihi_d,invDrehilohi_d,invDrelolohi_d,
                 invDrehihilo_d,invDrelohilo_d,invDrehilolo_d,invDrelololo_d,
                 invDimhihihi_d,invDimlohihi_d,invDimhilohi_d,invDimlolohi_d,
                 invDimhihilo_d,invDimlohilo_d,invDimhilolo_d,invDimlololo_d,
                  rhsrehihihi_d, rhsrelohihi_d, rhsrehilohi_d, rhsrelolohi_d,
                  rhsrehihilo_d, rhsrelohilo_d, rhsrehilolo_d, rhsrelololo_d,
                  rhsimhihihi_d, rhsimlohihi_d, rhsimhilohi_d, rhsimlolohi_d,
                  rhsimhihilo_d, rhsimlohilo_d, rhsimhilolo_d, rhsimlololo_d);
   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *mullapms += milliseconds;
   *totlapms += milliseconds;
   flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);

   int nbrUcol = (nbt-1)*szt*szt;               // #doubles in column of U
   double *Ucolrehihihi_h = new double[nbrUcol];  // column of U on host
   double *Ucolrelohihi_h = new double[nbrUcol];
   double *Ucolrehilohi_h = new double[nbrUcol];
   double *Ucolrelolohi_h = new double[nbrUcol];
   double *Ucolrehihilo_h = new double[nbrUcol];
   double *Ucolrelohilo_h = new double[nbrUcol];
   double *Ucolrehilolo_h = new double[nbrUcol];
   double *Ucolrelololo_h = new double[nbrUcol];
   double *Ucolimhihihi_h = new double[nbrUcol];
   double *Ucolimlohihi_h = new double[nbrUcol];
   double *Ucolimhilohi_h = new double[nbrUcol];
   double *Ucolimlolohi_h = new double[nbrUcol];
   double *Ucolimhihilo_h = new double[nbrUcol];
   double *Ucolimlohilo_h = new double[nbrUcol];
   double *Ucolimhilolo_h = new double[nbrUcol];
   double *Ucolimlololo_h = new double[nbrUcol];
   double *Ucolrehihihi_d;
   double *Ucolrelohihi_d;
   double *Ucolrehilohi_d;
   double *Ucolrelolohi_d;
   double *Ucolrehihilo_d;
   double *Ucolrelohilo_d;
   double *Ucolrehilolo_d;
   double *Ucolrelololo_d;
   double *Ucolimhihihi_d;
   double *Ucolimlohihi_d;
   double *Ucolimhilohi_d;
   double *Ucolimlolohi_d;
   double *Ucolimhihilo_d;
   double *Ucolimlohilo_d;
   double *Ucolimhilolo_d;
   double *Ucolimlololo_d;
   const size_t szUcol = nbrUcol*sizeof(double);
   cudaMalloc((void**)&Ucolrehihihi_d,szUcol);
   cudaMalloc((void**)&Ucolrelohihi_d,szUcol);
   cudaMalloc((void**)&Ucolrehilohi_d,szUcol);
   cudaMalloc((void**)&Ucolrelolohi_d,szUcol);
   cudaMalloc((void**)&Ucolrehihilo_d,szUcol);
   cudaMalloc((void**)&Ucolrelohilo_d,szUcol);
   cudaMalloc((void**)&Ucolrehilolo_d,szUcol);
   cudaMalloc((void**)&Ucolrelololo_d,szUcol);
   cudaMalloc((void**)&Ucolimhihihi_d,szUcol);
   cudaMalloc((void**)&Ucolimlohihi_d,szUcol);
   cudaMalloc((void**)&Ucolimhilohi_d,szUcol);
   cudaMalloc((void**)&Ucolimlolohi_d,szUcol);
   cudaMalloc((void**)&Ucolimhihilo_d,szUcol);
   cudaMalloc((void**)&Ucolimlohilo_d,szUcol);
   cudaMalloc((void**)&Ucolimhilolo_d,szUcol);
   cudaMalloc((void**)&Ucolimlololo_d,szUcol);

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
               Ucolrehihihi_h[ix]   = Urehihihi[rowoff+i][coloff+j];
               Ucolrelohihi_h[ix]   = Urelohihi[rowoff+i][coloff+j];
               Ucolrehilohi_h[ix]   = Urehilohi[rowoff+i][coloff+j];
               Ucolrelolohi_h[ix]   = Urelolohi[rowoff+i][coloff+j];
               Ucolrehihilo_h[ix]   = Urehihilo[rowoff+i][coloff+j];
               Ucolrelohilo_h[ix]   = Urelohilo[rowoff+i][coloff+j];
               Ucolrehilolo_h[ix]   = Urehilolo[rowoff+i][coloff+j];
               Ucolrelololo_h[ix]   = Urelololo[rowoff+i][coloff+j];
               Ucolimhihihi_h[ix]   = Uimhihihi[rowoff+i][coloff+j];
               Ucolimlohihi_h[ix]   = Uimlohihi[rowoff+i][coloff+j];
               Ucolimhilohi_h[ix]   = Uimhilohi[rowoff+i][coloff+j];
               Ucolimlolohi_h[ix]   = Uimlolohi[rowoff+i][coloff+j];
               Ucolimhihilo_h[ix]   = Uimhihilo[rowoff+i][coloff+j];
               Ucolimlohilo_h[ix]   = Uimlohilo[rowoff+i][coloff+j];
               Ucolimhilolo_h[ix]   = Uimhilolo[rowoff+i][coloff+j];
               Ucolimlololo_h[ix++] = Uimlololo[rowoff+i][coloff+j];
            }
      }
      cudaMemcpy(Ucolrehihihi_d,Ucolrehihihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelohihi_d,Ucolrelohihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrehilohi_d,Ucolrehilohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelolohi_d,Ucolrelolohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrehihilo_d,Ucolrehihilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelohilo_d,Ucolrelohilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrehilolo_d,Ucolrehilolo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolrelololo_d,Ucolrelololo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhihihi_d,Ucolimhihihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlohihi_d,Ucolimlohihi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhilohi_d,Ucolimhilohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlolohi_d,Ucolimlolohi_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhihilo_d,Ucolimhihilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlohilo_d,Ucolimlohilo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimhilolo_d,Ucolimhilolo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);
      cudaMemcpy(Ucolimlololo_d,Ucolimlololo_h,nbrUcol*sizeof(double),
                 cudaMemcpyHostToDevice);

      cudaEventRecord(start);
      cmplx8_back_substitute<<<k,szt>>>
         (szt,k,Ucolrehihihi_d,Ucolrelohihi_d,Ucolrehilohi_d,Ucolrelolohi_d,
                Ucolrehihilo_d,Ucolrelohilo_d,Ucolrehilolo_d,Ucolrelololo_d,
                Ucolimhihihi_d,Ucolimlohihi_d,Ucolimhilohi_d,Ucolimlolohi_d,
                Ucolimhihilo_d,Ucolimlohilo_d,Ucolimhilolo_d,Ucolimlololo_d,
                 rhsrehihihi_d, rhsrelohihi_d, rhsrehilohi_d, rhsrelolohi_d,
                 rhsrehihilo_d, rhsrelohilo_d, rhsrehilolo_d, rhsrelololo_d,
                 rhsimhihihi_d, rhsimlohihi_d, rhsimhilohi_d, rhsimlolohi_d,
                 rhsimhihilo_d, rhsimlohilo_d, rhsimhilolo_d, rhsimlololo_d);
      cudaEventRecord(stop);
      cudaEventSynchronize(stop);
      cudaEventElapsedTime(&milliseconds,start,stop);
      *sublapms += milliseconds;
      *totlapms += milliseconds;
      flopcount_cmplx_back_substitute(k,szt,addcnt,mulcnt);

      // (k-1)-th solution tile is ready for inverse multiplication
      cudaEventRecord(start);
      cmplx8_multiply_inverse<<<1,szt>>>
         (szt,k-1,invDrehihihi_d,invDrelohihi_d,invDrehilohi_d,invDrelolohi_d,
                  invDrehihilo_d,invDrelohilo_d,invDrehilolo_d,invDrelololo_d,
                  invDimhihihi_d,invDimlohihi_d,invDimhilohi_d,invDimlolohi_d,
                  invDimhihilo_d,invDimlohilo_d,invDimhilolo_d,invDimlololo_d,
                   rhsrehihihi_d, rhsrelohihi_d, rhsrehilohi_d, rhsrelolohi_d,
                   rhsrehihilo_d, rhsrelohilo_d, rhsrehilolo_d, rhsrelololo_d,
                   rhsimhihihi_d, rhsimlohihi_d, rhsimhilohi_d, rhsimlolohi_d,
                   rhsimhihilo_d, rhsimlohilo_d, rhsimhilolo_d, rhsimlololo_d);
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

   cudaMemcpy(xrehihihi,rhsrehihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelohihi,rhsrelohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrehilohi,rhsrehilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelolohi,rhsrelolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrehihilo,rhsrehihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelohilo,rhsrelohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrehilolo,rhsrehilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(xrelololo,rhsrelololo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhihihi,rhsimhihihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlohihi,rhsimlohihi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhilohi,rhsimhilohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlolohi,rhsimlolohi_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhihilo,rhsimhihilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlohilo,rhsimlohilo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximhilolo,rhsimhilolo_d,szrhs,cudaMemcpyDeviceToHost);
   cudaMemcpy(ximlololo,rhsimlololo_d,szrhs,cudaMemcpyDeviceToHost);

   // copy of invD_d is needed only for testing purposes
   cudaMemcpy(invDrehihihi_h,invDrehihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelohihi_h,invDrelohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrehilohi_h,invDrehilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelolohi_h,invDrelolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrehihilo_h,invDrehihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelohilo_h,invDrelohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrehilolo_h,invDrehilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDrelololo_h,invDrelololo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhihihi_h,invDimhihihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlohihi_h,invDimlohihi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhilohi_h,invDimhilohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlolohi_h,invDimlolohi_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhihilo_h,invDimhihilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlohilo_h,invDimlohilo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimhilolo_h,invDimhilolo_d,sznum,cudaMemcpyDeviceToHost);
   cudaMemcpy(invDimlololo_h,invDimlololo_d,sznum,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int k=0; k<nbt; k++) // copy rows of the inverse of the k-th tile
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Urehihihi[offset+i][offset+j] = invDrehihihi_h[ix];
            Urelohihi[offset+i][offset+j] = invDrelohihi_h[ix];
            Urehilohi[offset+i][offset+j] = invDrehilohi_h[ix];
            Urelolohi[offset+i][offset+j] = invDrelolohi_h[ix];
            Urehihilo[offset+i][offset+j] = invDrehihilo_h[ix];
            Urelohilo[offset+i][offset+j] = invDrelohilo_h[ix];
            Urehilolo[offset+i][offset+j] = invDrehilolo_h[ix];
            Urelololo[offset+i][offset+j] = invDrelololo_h[ix];
            Uimhihihi[offset+i][offset+j] = invDimhihihi_h[ix];
            Uimlohihi[offset+i][offset+j] = invDimlohihi_h[ix];
            Uimhilohi[offset+i][offset+j] = invDimhilohi_h[ix];
            Uimlolohi[offset+i][offset+j] = invDimlolohi_h[ix];
            Uimhihilo[offset+i][offset+j] = invDimhihilo_h[ix];
            Uimlohilo[offset+i][offset+j] = invDimlohilo_h[ix];
            Uimhilolo[offset+i][offset+j] = invDimhilolo_h[ix];
            Uimlololo[offset+i][offset+j] = invDimlololo_h[ix++];
         }
   }
   free(Drehihihi_h); free(Drelohihi_h); free(Drehilohi_h); free(Drelolohi_h);
   free(Drehihilo_h); free(Drelohilo_h); free(Drehilolo_h); free(Drelololo_h);
   free(Dimhihihi_h); free(Dimlohihi_h); free(Dimhilohi_h); free(Dimlolohi_h);
   free(Dimhihilo_h); free(Dimlohilo_h); free(Dimhilolo_h); free(Dimlololo_h);
   free(invDrehihihi_h); free(invDrelohihi_h);
   free(invDrehilohi_h); free(invDrelolohi_h);
   free(invDrehihilo_h); free(invDrelohilo_h);
   free(invDrehilolo_h); free(invDrelololo_h);
   free(Ucolrehihihi_h); free(Ucolrelohihi_h);
   free(Ucolrehilohi_h); free(Ucolrelolohi_h);
   free(Ucolrehihilo_h); free(Ucolrelohilo_h);
   free(Ucolrehilolo_h); free(Ucolrelololo_h);
   free(invDimhihihi_h); free(invDimlohihi_h);
   free(invDimhilohi_h); free(invDimlolohi_h);
   free(invDimhihilo_h); free(invDimlohilo_h);
   free(invDimhilolo_h); free(invDimlololo_h);
   free(Ucolimhihihi_h); free(Ucolimlohihi_h);
   free(Ucolimhilohi_h); free(Ucolimlolohi_h);
   free(Ucolimhihilo_h); free(Ucolimlohilo_h);
   free(Ucolimhilolo_h); free(Ucolimlololo_h);
}
