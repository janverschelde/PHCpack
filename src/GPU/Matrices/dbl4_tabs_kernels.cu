/* The file dbl4_tabs_kernels.cu defines the functions specified in
 * the file dbl4_tabs_kernels.h. */

#include <iostream>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#ifdef gpufun
#include "double_double_gpufun.cu"
#include "quad_double_gpufun.cu"
#endif
#include "dbl4_tabs_kernels.h"
#include "dbl_tabs_flopcounts.h"

using namespace std;

__global__ void dbl4_small_invert_upper 
( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
  double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolhihi[tabsqd_shmemsize];
   __shared__ double Ucollohi[tabsqd_shmemsize];
   __shared__ double Ucolhilo[tabsqd_shmemsize];
   __shared__ double Ucollolo[tabsqd_shmemsize];
   __shared__ double invUrowshihi[tabsqd_shmemsize];
   __shared__ double invUrowslohi[tabsqd_shmemsize];
   __shared__ double invUrowshilo[tabsqd_shmemsize];
   __shared__ double invUrowslolo[tabsqd_shmemsize];

   double rhshihi,rhslohi,rhshilo,rhslolo;
   double xvalhihi,xvallohi,xvalhilo,xvallolo;
   double acchihi,acclohi,acchilo,acclolo;

   int colidx = dim*(dim-1);          // start with the last column

   Ucolhihi[k] = Uhihi[colidx+k];     // load the last column
   Ucollohi[k] = Ulohi[colidx+k];
   Ucolhilo[k] = Uhilo[colidx+k];
   Ucollolo[k] = Ulolo[colidx+k];

   rhshihi = ((double) int(k == dim-1)); // right hand side for each thread
   rhslohi = 0.0;
   rhshilo = 0.0;
   rhslolo = 0.0;
   int rowidx = (dim - 1)*dim + k;      // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k]; // last row of the inverse
   qdg_div(rhshihi,rhslohi,rhshilo,rhslolo,
           Ucolhihi[k],Ucollohi[k],Ucolhilo[k],Ucollolo[k],
           &invUrowshihi[rowidx],&invUrowslohi[rowidx],
           &invUrowshilo[rowidx],&invUrowslolo[rowidx]);

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshihi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslohi = 0.0;
      rhshilo = 0.0;
      rhslolo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U

         Ucolhihi[k] = Uhihi[colidx+k];
         Ucollohi[k] = Ulohi[colidx+k];
         Ucolhilo[k] = Uhilo[colidx+k];
         Ucollolo[k] = Ulolo[colidx+k];

         rowidx = j*dim + k;          // need solution value

         xvalhihi = invUrowshihi[rowidx];
         xvallohi = invUrowslohi[rowidx];
         xvalhilo = invUrowshilo[rowidx];
         xvallolo = invUrowslolo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
         qdg_mul(Ucolhihi[i],Ucollohi[i],Ucolhilo[i],Ucollolo[i],
                 xvalhihi,   xvallohi,   xvalhilo,   xvallolo,
                 &acchihi,   &acclohi,   &acchilo,   &acclolo);
         qdg_dec(&rhshihi,&rhslohi,&rhshilo,&rhslolo,
                  acchihi, acclohi, acchilo, acclolo);
      }
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolhihi[k] = Uhihi[colidx+k];
      Ucolhilo[k] = Uhilo[colidx+k];
      Ucollohi[k] = Ulohi[colidx+k];
      Ucollolo[k] = Ulolo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
      qdg_div(rhshihi,rhslohi,rhshilo,rhslolo,
              Ucolhihi[i],Ucollohi[i],Ucolhilo[i],Ucollolo[i],
              &invUrowshihi[rowidx],&invUrowslohi[rowidx],
              &invUrowshilo[rowidx],&invUrowslolo[rowidx]);
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUhihi[rowidx+k] = invUrowshihi[rowidx+k];
      invUlohi[rowidx+k] = invUrowslohi[rowidx+k];
      invUhilo[rowidx+k] = invUrowshilo[rowidx+k];
      invUlolo[rowidx+k] = invUrowslolo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void cmplx4_small_invert_upper
 ( int dim,
   double *Urehihi, double *Urelohi, double *Urehilo, double *Urelolo,
   double *Uimhihi, double *Uimlohi, double *Uimhilo, double *Uimlolo,
   double *invUrehihi, double *invUrelohi,
   double *invUrehilo, double *invUrelolo,
   double *invUimhihi, double *invUimlohi,
   double *invUimhilo, double *invUimlolo )
{
   const int k = threadIdx.x; // thread k computes k-th column of inverse

   __shared__ double Ucolrehihi[tabsqd_shmemsize];
   __shared__ double Ucolrelohi[tabsqd_shmemsize];
   __shared__ double Ucolrehilo[tabsqd_shmemsize];
   __shared__ double Ucolrelolo[tabsqd_shmemsize];
   __shared__ double Ucolimhihi[tabsqd_shmemsize];
   __shared__ double Ucolimlohi[tabsqd_shmemsize];
   __shared__ double Ucolimhilo[tabsqd_shmemsize];
   __shared__ double Ucolimlolo[tabsqd_shmemsize];
   __shared__ double invUrowsrehihi[tabsqd_shmemsize];
   __shared__ double invUrowsrelohi[tabsqd_shmemsize];
   __shared__ double invUrowsrehilo[tabsqd_shmemsize];
   __shared__ double invUrowsrelolo[tabsqd_shmemsize];
   __shared__ double invUrowsimhihi[tabsqd_shmemsize];
   __shared__ double invUrowsimlohi[tabsqd_shmemsize];
   __shared__ double invUrowsimhilo[tabsqd_shmemsize];
   __shared__ double invUrowsimlolo[tabsqd_shmemsize];

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

   int colidx = dim*(dim-1);             // start with the last column

   Ucolrehihi[k] = Urehihi[colidx+k];    // load the last column
   Ucolrelohi[k] = Urelohi[colidx+k];
   Ucolrehilo[k] = Urehilo[colidx+k];    // load the last column
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
   int rowidx = (dim - 1)*dim + k;       // the row index in the inverse

   __syncthreads();
   // invUrows[rowidx] = rhs/Ucol[k];    // last row of the inverse
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
   invUrowsrehihi[rowidx] = acc1hihi;
   invUrowsrelohi[rowidx] = acc1lohi;
   invUrowsrehilo[rowidx] = acc1hilo;
   invUrowsrelolo[rowidx] = acc1lolo;
   qdg_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
            acc4hihi, acc4lohi, acc4hilo, acc4lolo);
   invUrowsimhihi[rowidx] = acc3hihi;
   invUrowsimlohi[rowidx] = acc3lohi;
   invUrowsimhilo[rowidx] = acc3hilo;
   invUrowsimlolo[rowidx] = acc3lolo;

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

         rowidx = j*dim + k;          // need solution value

         xvalrehihi = invUrowsrehihi[rowidx];
         xvalrelohi = invUrowsrelohi[rowidx];
         xvalrehilo = invUrowsrehilo[rowidx];
         xvalrelolo = invUrowsrelolo[rowidx];
         xvalimhihi = invUrowsimhihi[rowidx];
         xvalimlohi = invUrowsimlohi[rowidx];
         xvalimhilo = invUrowsimhilo[rowidx];
         xvalimlolo = invUrowsimlolo[rowidx];

         __syncthreads();
         // rhs = rhs - Ucol[i]*xval; // update right hand side
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
      rowidx = i*dim + k;             // save in i-th row of inverse

      colidx = dim*i;                 // need column i of U
      Ucolrehihi[k] = Urehihi[colidx+k];
      Ucolrelohi[k] = Urelohi[colidx+k];
      Ucolrehilo[k] = Urehilo[colidx+k];
      Ucolrelolo[k] = Urelolo[colidx+k];
      Ucolimhihi[k] = Uimhihi[colidx+k];
      Ucolimlohi[k] = Uimlohi[colidx+k];
      Ucolimhilo[k] = Uimhilo[colidx+k];
      Ucolimlolo[k] = Uimlolo[colidx+k];

      __syncthreads();
      // invUrows[rowidx] = rhs/Ucol[i];
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
               acc2hihi,acc2lohi,acc2hilo,acc2lolo);
      invUrowsrehihi[rowidx] = acc1hihi;
      invUrowsrelohi[rowidx] = acc1lohi;
      invUrowsrehilo[rowidx] = acc1hilo;
      invUrowsrelolo[rowidx] = acc1lolo;
      qdg_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
               acc4hihi, acc4lohi, acc4hilo, acc4lolo);
      invUrowsimhihi[rowidx] = acc3hihi;
      invUrowsimlohi[rowidx] = acc3lohi;
      invUrowsimhilo[rowidx] = acc3hilo;
      invUrowsimlolo[rowidx] = acc3lolo;
   }
   rowidx = 0;
   for(int i=0; i<dim; i++)
   {
      __syncthreads();
      invUrehihi[rowidx+k] = invUrowsrehihi[rowidx+k];
      invUrelohi[rowidx+k] = invUrowsrelohi[rowidx+k];
      invUrehilo[rowidx+k] = invUrowsrehilo[rowidx+k];
      invUrelolo[rowidx+k] = invUrowsrelolo[rowidx+k];
      invUimhihi[rowidx+k] = invUrowsimhihi[rowidx+k];
      invUimlohi[rowidx+k] = invUrowsimlohi[rowidx+k];
      invUimhilo[rowidx+k] = invUrowsimhilo[rowidx+k];
      invUimlolo[rowidx+k] = invUrowsimlolo[rowidx+k];
      rowidx = rowidx + dim;
   }
}

__global__ void dbl4_medium_invert_upper
 ( int dim, double *Uhihi, double *Ulohi, double *Uhilo, double *Ulolo,
   double *invUhihi, double *invUlohi, double *invUhilo, double *invUlolo)
{
   const int k = threadIdx.x;  // thread k computes k-th column of inverse

   __shared__ double Ucolhihi[tabsqd_shmemsize];      // one column of U
   __shared__ double Ucollohi[tabsqd_shmemsize];      // one column of U
   __shared__ double Ucolhilo[tabsqd_shmemsize];      // one column of U
   __shared__ double Ucollolo[tabsqd_shmemsize];      // one column of U
   __shared__ double invUrowhihi[tabsqd_shmemsize];   // one row of invU
   __shared__ double invUrowlohi[tabsqd_shmemsize];   // one row of invU
   __shared__ double invUrowhilo[tabsqd_shmemsize];   // one row of invU
   __shared__ double invUrowlolo[tabsqd_shmemsize];   // one row of invU

   double rhshihi,rhslohi,rhshilo,rhslolo;
   double xvalhihi,xvallohi,xvalhilo,xvallolo;
   double acchihi,acclohi,acchilo,acclolo;

   int colidx = dim*(dim-1);           // start with the last column

   Ucolhihi[k] = Uhihi[colidx+k];      // load the last column
   Ucollohi[k] = Ulohi[colidx+k];
   Ucolhilo[k] = Uhilo[colidx+k];      // load the last column
   Ucollolo[k] = Ulolo[colidx+k];
   rhshihi = ((double) int(k == dim-1)); // right hand side for each thread
   rhslohi = 0.0;
   rhshilo = 0.0;
   rhslolo = 0.0;
   int rowidx = (dim - 1)*dim + k;       // the row index in the inverse

   // invUrow[k] = rhs/Ucol[k];          // last row of the inverse
   qdg_div( rhshihi,        rhslohi,        rhshilo,        rhslolo,
           Ucolhihi[k],    Ucollohi[k],    Ucolhilo[k],    Ucollolo[k],
       &invUrowhihi[k],&invUrowlohi[k],&invUrowhilo[k],&invUrowlolo[k]);
   invUhihi[rowidx] = invUrowhihi[k];    // store the last row into invU
   invUlohi[rowidx] = invUrowlohi[k]; 
   invUhilo[rowidx] = invUrowhilo[k];
   invUlolo[rowidx] = invUrowlolo[k]; 

   for(int i=dim-2; i>=0; i--)        // compute row with index i
   {
      rhshihi = ((double) int(k == i)); // set rhs for i-th unit vector
      rhslohi = 0.0;
      rhshilo = 0.0;
      rhslolo = 0.0;

      for(int j=i+1; j<dim; j++)
      {
         colidx = dim*j;              // need column j of U
         Ucolhihi[k] = Uhihi[colidx+k];
         Ucollohi[k] = Ulohi[colidx+k];
         Ucolhilo[k] = Uhilo[colidx+k];
         Ucollolo[k] = Ulolo[colidx+k];

         rowidx = j*dim + k;                // need solution value
         invUrowhihi[k] = invUhihi[rowidx]; // load invU row into invUrow
         invUrowlohi[k] = invUlohi[rowidx];
         invUrowhilo[k] = invUhilo[rowidx]; // load invU row into invUrow
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
      colidx = dim*i;                 // need column i of U
      Ucolhihi[k] = Uhihi[colidx+k];
      Ucollohi[k] = Ulohi[colidx+k];
      Ucolhilo[k] = Uhilo[colidx+k];
      Ucollolo[k] = Ulolo[colidx+k];
      rowidx = i*dim + k;             // save in i-th row of inverse

      __syncthreads();
      // invUrow[k] = rhs/Ucol[i];
      qdg_div(rhshihi,rhslohi,rhshilo,rhslolo,
              Ucolhihi[i],Ucollohi[i],Ucolhilo[i],Ucollolo[i],
              &invUrowhihi[k],&invUrowlohi[k],
              &invUrowhilo[k],&invUrowlolo[k]);
      invUhihi[rowidx] = invUrowhihi[k];
      invUlohi[rowidx] = invUrowlohi[k];
      invUhilo[rowidx] = invUrowhilo[k];
      invUlolo[rowidx] = invUrowlolo[k];
   }
}

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
   invUhihi[rowidx] = 0.0;
   invUlohi[rowidx] = 0.0;     // initialize in case of zero divisor
   invUhilo[rowidx] = 0.0;
   invUlolo[rowidx] = 0.0;
   if(1.0 + Ucolhihi[k] != 1.0)
   {
      qdg_div(     rhshihi,        rhslohi,        rhshilo,        rhslolo,
                  Ucolhihi[k],    Ucollohi[k],    Ucolhilo[k],    Ucollolo[k],
              &invUrowhihi[k],&invUrowlohi[k],&invUrowhilo[k],&invUrowlolo[k]);
      invUhihi[rowidx] = invUrowhihi[k];     // store the last row into invU
      invUlohi[rowidx] = invUrowlohi[k];
      invUhilo[rowidx] = invUrowhilo[k];
      invUlolo[rowidx] = invUrowlolo[k];
   }
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
      invUhihi[rowidx] = 0.0;
      invUlohi[rowidx] = 0.0;    // initialize in case of zero divisor
      invUhilo[rowidx] = 0.0;
      invUlolo[rowidx] = 0.0;

      if(1.0 + Ucolhihi[i] != 1.0)
      {
         qdg_div(     rhshihi,        rhslohi,        rhshilo,        rhslolo,
                  Ucolhihi[i],    Ucollohi[i],    Ucolhilo[i],    Ucollolo[i],
              &invUrowhihi[k],&invUrowlohi[k],&invUrowhilo[k],&invUrowlolo[k]);
         invUhihi[rowidx] = invUrowhihi[k];
         invUlohi[rowidx] = invUrowlohi[k];
         invUhilo[rowidx] = invUrowhilo[k];
         invUlolo[rowidx] = invUrowlolo[k];
      }
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

   invUrehihi[rowidx] = 0.0;  // initialize in case of zero denominator
   invUrelohi[rowidx] = 0.0;
   invUrehilo[rowidx] = 0.0; 
   invUrelolo[rowidx] = 0.0;
   invUimhihi[rowidx] = 0.0;
   invUimlohi[rowidx] = 0.0;
   invUimhilo[rowidx] = 0.0;
   invUimlolo[rowidx] = 0.0;

   if(1.0 + denhihi != 1.0)
   {
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
      invUrehihi[rowidx] = invUrowrehihi[k];  // store the last row into invU
      invUrelohi[rowidx] = invUrowrelohi[k];
      invUrehilo[rowidx] = invUrowrehilo[k]; 
      invUrelolo[rowidx] = invUrowrelolo[k];
      invUimhihi[rowidx] = invUrowimhihi[k];
      invUimlohi[rowidx] = invUrowimlohi[k];
      invUimhilo[rowidx] = invUrowimhilo[k];
      invUimlolo[rowidx] = invUrowimlolo[k];
   }
   __syncthreads();
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

      invUrehihi[rowidx] = 0.0; // initialize in case of zero denominator
      invUrelohi[rowidx] = 0.0;
      invUrehilo[rowidx] = 0.0;
      invUrelolo[rowidx] = 0.0;
      invUimhihi[rowidx] = 0.0;
      invUimlohi[rowidx] = 0.0;
      invUimhilo[rowidx] = 0.0;
      invUimlolo[rowidx] = 0.0;

      if(1.0 + denhihi != 1.0)
      {
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
   double coeffhihi,coefflohi,coeffhilo,coefflolo;
   double acchihi,acclohi,acchilo,acclolo;

   for(int j=0; j<dim; j++)  // column j of the inverse diagonal tile
   {
      coeffhihi = Urehihi[offset+k*dim+j];
      coefflohi = Urelohi[offset+k*dim+j];
      coeffhilo = Urehilo[offset+k*dim+j];
      coefflolo = Urelolo[offset+k*dim+j];
      qdg_mul(coeffhihi, coefflohi, coeffhilo, coefflolo,
                solrehihi[j],solrelohi[j],solrehilo[j],solrelolo[j],
                &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      qdg_mul(coeffhihi, coefflohi, coeffhilo, coefflolo,
                solimhihi[j],solimlohi[j],solimhilo[j],solimlolo[j],
                 &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      coeffhihi = Uimhihi[offset+k*dim+j];
      coefflohi = Uimlohi[offset+k*dim+j];
      coeffhilo = Uimhilo[offset+k*dim+j];
      coefflolo = Uimlolo[offset+k*dim+j];
      // result = result + coeff*sol[j];
      qdg_mul(coeffhihi, coefflohi, coeffhilo, coefflolo,
                solimhihi[j],solimlohi[j],solimhilo[j],solimlolo[j],
                 &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_dec(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
      qdg_mul(coeffhihi, coefflohi, coeffhilo, coefflolo,
                solrehihi[j],solrelohi[j],solrehilo[j],solrelolo[j],
                 &acchihi,    &acclohi,    &acchilo,    &acclolo);
      qdg_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
                    acchihi,      acclohi,      acchilo,      acclolo);
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

void GPU_dbl4_upper_inverse
 ( int dim, double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double **invUhihi, double **invUlohi, double **invUhilo, double **invUlolo,
   double *lapms, double *walltimesec )
{
   const int szU = dim*dim;

   double *Uhihi_h = new double[szU];   // the columns of U
   double *Ulohi_h = new double[szU];    
   double *Uhilo_h = new double[szU];   
   double *Ulolo_h = new double[szU]; 
   double *Uhihi_d;                     // the columns of U on the device
   double *Ulohi_d;
   double *Uhilo_d;
   double *Ulolo_d;
   double *invUhihi_h = new double[szU];  // columns of the inverse of U
   double *invUlohi_h = new double[szU]; 
   double *invUhilo_h = new double[szU]; 
   double *invUlolo_h = new double[szU]; 
   double *invUhihi_d;                    // inverse of U on the device
   double *invUlohi_d;
   double *invUhilo_d;
   double *invUlolo_d;

   int ix = 0;
   for(int j=0; j<dim; j++)
      for(int i=0; i<dim; i++)
      {
         Uhihi_h[ix]   = Uhihi[i][j];
         Ulohi_h[ix]   = Ulohi[i][j];
         Uhilo_h[ix]   = Uhilo[i][j];
         Ulolo_h[ix++] = Ulolo[i][j];
      }

   size_t szmat = szU*sizeof(double);
   cudaMalloc((void**)&Uhihi_d,szmat);
   cudaMalloc((void**)&Ulohi_d,szmat);
   cudaMalloc((void**)&Uhilo_d,szmat);
   cudaMalloc((void**)&Ulolo_d,szmat);
   cudaMalloc((void**)&invUhihi_d,szmat);
   cudaMalloc((void**)&invUlohi_d,szmat);
   cudaMalloc((void**)&invUhilo_d,szmat);
   cudaMalloc((void**)&invUlolo_d,szmat);
   cudaMemcpy(Uhihi_d,Uhihi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulohi_d,Ulohi_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Uhilo_d,Uhilo_h,szmat,cudaMemcpyHostToDevice);
   cudaMemcpy(Ulolo_d,Ulolo_h,szmat,cudaMemcpyHostToDevice);

   cudaEvent_t start,stop;           // to measure time spent by kernels 
   cudaEventCreate(&start);
   cudaEventCreate(&stop);
   *lapms = 0.0;
   float milliseconds;
   struct timeval begintime,endtime; // wall clock time of computations

   gettimeofday(&begintime,0);

   cudaEventRecord(start);

   if(dim <= 16)
      dbl4_small_invert_upper<<<1,dim>>>
         (dim,   Uhihi_d,   Ulohi_d,   Uhilo_d,   Ulolo_d,
              invUhihi_d,invUlohi_d,invUhilo_d,invUlolo_d);
   else
      dbl4_medium_invert_upper<<<1,dim>>>
         (dim,   Uhihi_d,   Ulohi_d,   Uhilo_d,   Ulolo_d,
              invUhihi_d,invUlohi_d,invUhilo_d,invUlolo_d);

   cudaEventRecord(stop);
   cudaEventSynchronize(stop);
   cudaEventElapsedTime(&milliseconds,start,stop);
   *lapms += milliseconds;

   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   *walltimesec = seconds + microseconds*1.0e-6;

   cudaMemcpy(invUhihi_h,invUhihi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlohi_h,invUlohi_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUhilo_h,invUhilo_d,szmat,cudaMemcpyDeviceToHost);
   cudaMemcpy(invUlolo_h,invUlolo_d,szmat,cudaMemcpyDeviceToHost);

   ix = 0;
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         invUhihi[i][j] = invUhihi_h[ix];
         invUlohi[i][j] = invUlohi_h[ix];
         invUhilo[i][j] = invUhilo_h[ix];
         invUlolo[i][j] = invUlolo_h[ix++];
      }

   free(Uhihi_h); free(Ulohi_h); free(Uhilo_h); free(Ulolo_h);
   free(invUhihi_h); free(invUlohi_h); free(invUhilo_h); free(invUlolo_h);
}

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
   long long int *addcnt, double *addover,
   long long int *mulcnt, double *mulover,
   long long int *divcnt, double *divover )
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
   *addcnt = 0; *addover = 0.0;
   *mulcnt = 0; *mulover = 0.0;
   *divcnt = 0; *divover = 0.0;

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
   // flopcount_dbl_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);
   overflopcount_dbl_invert_tiles
      (nbt,szt,addcnt,addover,mulcnt,mulover,divcnt,divover);

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
   // flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);
   overflopcount_dbl_multiply_inverse(szt,addcnt,addover,mulcnt,mulover);

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
      // flopcount_dbl_back_substitute(k,szt,addcnt,mulcnt);
      overflopcount_dbl_back_substitute(k,szt,addcnt,addover,mulcnt,mulover);

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
      // flopcount_dbl_multiply_inverse(szt,addcnt,mulcnt);
      overflopcount_dbl_multiply_inverse(szt,addcnt,addover,mulcnt,mulover);

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
   long long int *addcnt, double *addover,
   long long int *mulcnt, double *mulover,
   long long int *divcnt, double *divover )
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
   *addcnt = 0; *addover = 0.0;
   *mulcnt = 0; *mulover = 0.0;
   *divcnt = 0; *divover = 0.0;

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
   // flopcount_cmplx_invert_tiles(nbt,szt,addcnt,mulcnt,divcnt);
   overflopcount_cmplx_invert_tiles
      (nbt,szt,addcnt,addover,mulcnt,mulover,divcnt,divover);

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
   // flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);
   overflopcount_cmplx_multiply_inverse(szt,addcnt,addover,mulcnt,mulover);

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
      // flopcount_cmplx_back_substitute(k,szt,addcnt,mulcnt);
      overflopcount_cmplx_back_substitute(k,szt,addcnt,addover,mulcnt,mulover);

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
      // flopcount_cmplx_multiply_inverse(szt,addcnt,mulcnt);
      overflopcount_cmplx_multiply_inverse(szt,addcnt,addover,mulcnt,mulover);

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
