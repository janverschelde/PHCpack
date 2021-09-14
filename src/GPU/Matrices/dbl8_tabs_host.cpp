/* The file dbl8_tabs_host.cpp defines functions specified in
 * the file dbl8_tabs_host.h. */

#include <cstdlib>
#include <ctime>
#include "octo_double_functions.h"
#include "dbl8_factorizations.h"
#include "dbl8_tabs_host.h"

using namespace std;

void CPU_dbl8_backsubs
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo )
{
   CPU_dbl8_factors_backward
      (dim,Uhihihi,Ulohihi,Uhilohi,Ulolohi,
           Uhihilo,Ulohilo,Uhilolo,Ulololo,
           bhihihi,blohihi,bhilohi,blolohi,
           bhihilo,blohilo,bhilolo,blololo,
           xhihihi,xlohihi,xhilohi,xlolohi,
           xhihilo,xlohilo,xhilolo,xlololo);
}

void CPU_cmplx8_backsubs
 ( int dim,
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
   double *ximhihilo, double *ximlohilo, double *ximhilolo, double *ximlololo )
{
   CPU_cmplx8_factors_backward
      (dim,Urehihihi,Urelohihi,Urehilohi,Urelolohi,
           Urehihilo,Urelohilo,Urehilolo,Urelololo,
           Uimhihihi,Uimlohihi,Uimhilohi,Uimlolohi,
           Uimhihilo,Uimlohilo,Uimhilolo,Uimlololo,
           brehihihi,brelohihi,brehilohi,brelolohi,
           brehihilo,brelohilo,brehilolo,brelololo,
           bimhihihi,bimlohihi,bimhilohi,bimlolohi,
           bimhihilo,bimlohilo,bimhilolo,bimlololo,
           xrehihihi,xrelohihi,xrehilohi,xrelolohi,
           xrehihilo,xrelohilo,xrehilolo,xrelololo,
           ximhihihi,ximlohihi,ximhilohi,ximlolohi,
           ximhihilo,ximlohilo,ximhilolo,ximlololo);
}

void CPU_dbl8_upper_inverse
 ( int dim,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double **invUhihihi, double **invUlohihi,
   double **invUhilohi, double **invUlolohi,
   double **invUhihilo, double **invUlohilo,
   double **invUhilolo, double **invUlololo, double *lapsec )
{
   double *colhihihi = new double[dim];
   double *collohihi = new double[dim];
   double *colhilohi = new double[dim];
   double *collolohi = new double[dim];
   double *colhihilo = new double[dim];
   double *collohilo = new double[dim];
   double *colhilolo = new double[dim];
   double *collololo = new double[dim];
   double *rhshihihi = new double[dim];
   double *rhslohihi = new double[dim];
   double *rhshilohi = new double[dim];
   double *rhslolohi = new double[dim];
   double *rhshihilo = new double[dim];
   double *rhslohilo = new double[dim];
   double *rhshilolo = new double[dim];
   double *rhslololo = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhshihihi[i] = 0.0; rhslohihi[i] = 0.0;
      rhshilohi[i] = 0.0; rhslolohi[i] = 0.0;
      rhshihilo[i] = 0.0; rhslohilo[i] = 0.0;
      rhshilolo[i] = 0.0; rhslololo[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhshihihi[j] = 1.0;

      CPU_dbl8_backsubs
         (dim,Uhihihi,  Ulohihi,  Uhilohi,  Ulolohi,
              Uhihilo,  Ulohilo,  Uhilolo,  Ulololo,
            rhshihihi,rhslohihi,rhshilohi,rhslolohi,
            rhshihilo,rhslohilo,rhshilolo,rhslololo,
            colhihihi,collohihi,colhilohi,collolohi,
            colhihilo,collohilo,colhilolo,collololo);

      for(int i=0; i<dim; i++)
      {
         invUhihihi[i][j] = colhihihi[i]; invUlohihi[i][j] = collohihi[i];
         invUhilohi[i][j] = colhilohi[i]; invUlolohi[i][j] = collolohi[i];
         invUhihilo[i][j] = colhihilo[i]; invUlohilo[i][j] = collohilo[i];
         invUhilolo[i][j] = colhilolo[i]; invUlololo[i][j] = collololo[i];
      }
      rhshihihi[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhshihihi); free(rhslohihi); free(rhshilohi); free(rhslolohi);
   free(rhshihilo); free(rhslohilo); free(rhshilolo); free(rhslololo);
   free(colhihihi); free(collohihi); free(colhilohi); free(collolohi);
   free(colhihilo); free(collohilo); free(colhilolo); free(collololo);
}

void CPU_cmplx8_upper_inverse
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
   double **invUimhilolo, double **invUimlololo, double *lapsec )
{
   double *colrehihihi = new double[dim];
   double *colrelohihi = new double[dim];
   double *colrehilohi = new double[dim];
   double *colrelolohi = new double[dim];
   double *colrehihilo = new double[dim];
   double *colrelohilo = new double[dim];
   double *colrehilolo = new double[dim];
   double *colrelololo = new double[dim];
   double *colimhihihi = new double[dim];
   double *colimlohihi = new double[dim];
   double *colimhilohi = new double[dim];
   double *colimlolohi = new double[dim];
   double *colimhihilo = new double[dim];
   double *colimlohilo = new double[dim];
   double *colimhilolo = new double[dim];
   double *colimlololo = new double[dim];
   double *rhsrehihihi = new double[dim];
   double *rhsrelohihi = new double[dim];
   double *rhsrehilohi = new double[dim];
   double *rhsrelolohi = new double[dim];
   double *rhsrehihilo = new double[dim];
   double *rhsrelohilo = new double[dim];
   double *rhsrehilolo = new double[dim];
   double *rhsrelololo = new double[dim];
   double *rhsimhihihi = new double[dim];
   double *rhsimlohihi = new double[dim];
   double *rhsimhilohi = new double[dim];
   double *rhsimlolohi = new double[dim];
   double *rhsimhihilo = new double[dim];
   double *rhsimlohilo = new double[dim];
   double *rhsimhilolo = new double[dim];
   double *rhsimlololo = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhsrehihihi[i] = 0.0; rhsrelohihi[i] = 0.0;
      rhsrehilohi[i] = 0.0; rhsrelolohi[i] = 0.0;
      rhsrehihilo[i] = 0.0; rhsrelohilo[i] = 0.0;
      rhsrehilolo[i] = 0.0; rhsrelololo[i] = 0.0;
      rhsimhihihi[i] = 0.0; rhsimlohihi[i] = 0.0;
      rhsimhilohi[i] = 0.0; rhsimlolohi[i] = 0.0;
      rhsimhihilo[i] = 0.0; rhsimlohilo[i] = 0.0;
      rhsimhilolo[i] = 0.0; rhsimlololo[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhsrehihihi[j] = 1.0;

      CPU_cmplx8_backsubs
         (dim, Urehihihi,  Urelohihi,  Urehilohi,  Urelolohi,
               Urehihilo,  Urelohilo,  Urehilolo,  Urelololo,
               Uimhihihi,  Uimlohihi,  Uimhilohi,  Uimlolohi,
               Uimhihilo,  Uimlohilo,  Uimhilolo,  Uimlololo,
             rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
             rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
             rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
             rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo,
             colrehihihi,colrelohihi,colrehilohi,colrelolohi,
             colrehihilo,colrelohilo,colrehilolo,colrelololo,
             colimhihihi,colimlohihi,colimhilohi,colimlolohi,
             colimhihilo,colimlohilo,colimhilolo,colimlololo);

      for(int i=0; i<dim; i++)
      {
         invUrehihihi[i][j] = colrehihihi[i];
         invUrelohihi[i][j] = colrelohihi[i];
         invUrehilohi[i][j] = colrehilohi[i];
         invUrelolohi[i][j] = colrelolohi[i];
         invUrehihilo[i][j] = colrehihilo[i];
         invUrelohilo[i][j] = colrelohilo[i];
         invUrehilolo[i][j] = colrehilolo[i];
         invUrelololo[i][j] = colrelololo[i];
         invUimhihihi[i][j] = colimhihihi[i];
         invUimlohihi[i][j] = colimlohihi[i];
         invUimhilohi[i][j] = colimhilohi[i];
         invUimlolohi[i][j] = colimlolohi[i];
         invUimhihilo[i][j] = colimhihilo[i];
         invUimlohilo[i][j] = colimlohilo[i];
         invUimhilolo[i][j] = colimhilolo[i];
         invUimlololo[i][j] = colimlololo[i];
      }
      rhsrehihihi[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhsrehihihi); free(rhsrelohihi); free(rhsrehilohi); free(rhsrelolohi);
   free(rhsrehihilo); free(rhsrelohilo); free(rhsrehilolo); free(rhsrelololo);
   free(colrehihihi); free(colrelohihi); free(colrehilohi); free(colrelolohi);
   free(colrehihilo); free(colrelohilo); free(colrehilolo); free(colrelololo);
   free(rhsimhihihi); free(rhsimlohihi); free(rhsimhilohi); free(rhsimlolohi);
   free(rhsimhihilo); free(rhsimlohilo); free(rhsimhilolo); free(rhsimlololo);
   free(colimhihihi); free(colimlohihi); free(colimhilohi); free(colimlolohi);
   free(colimhihilo); free(colimlohilo); free(colimhilolo); free(colimlololo);
}

void CPU_dbl8_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihihi, double **Ulohihi, double **Uhilohi, double **Ulolohi,
   double **Uhihilo, double **Ulohilo, double **Uhilolo, double **Ulololo,
   double *bhihihi, double *blohihi, double *bhilohi, double *blolohi,
   double *bhihilo, double *blohilo, double *bhilolo, double *blololo,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *lapsec )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;
   double prodhihihi,prodlohihi,prodhilohi,prodlolohi;
   double prodhihilo,prodlohilo,prodhilolo,prodlololo;
   double **Thihihi = new double*[szt];
   double **Tlohihi = new double*[szt];
   double **Thilohi = new double*[szt];
   double **Tlolohi = new double*[szt];
   double **Thihilo = new double*[szt];
   double **Tlohilo = new double*[szt];
   double **Thilolo = new double*[szt];
   double **Tlololo = new double*[szt];
   double **invThihihi = new double*[szt];
   double **invTlohihi = new double*[szt];
   double **invThilohi = new double*[szt];
   double **invTlolohi = new double*[szt];
   double **invThihilo = new double*[szt];
   double **invTlohilo = new double*[szt];
   double **invThilolo = new double*[szt];
   double **invTlololo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Thihihi[i] = new double[szt];
      Tlohihi[i] = new double[szt];
      Thilohi[i] = new double[szt];
      Tlolohi[i] = new double[szt];
      Thihilo[i] = new double[szt];
      Tlohilo[i] = new double[szt];
      Thilolo[i] = new double[szt];
      Tlololo[i] = new double[szt];
      invThihihi[i] = new double[szt];
      invTlohihi[i] = new double[szt];
      invThilohi[i] = new double[szt];
      invTlolohi[i] = new double[szt];
      invThihilo[i] = new double[szt];
      invTlohilo[i] = new double[szt];
      invThilolo[i] = new double[szt];
      invTlololo[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Thihihi[i][j] = Uhihihi[idx+i][idx+j];
         Tlohihi[i][j] = Ulohihi[idx+i][idx+j];
         Thilohi[i][j] = Uhilohi[idx+i][idx+j];
         Tlolohi[i][j] = Ulolohi[idx+i][idx+j];
         Thihilo[i][j] = Uhihilo[idx+i][idx+j];
         Tlohilo[i][j] = Ulohilo[idx+i][idx+j];
         Thilolo[i][j] = Uhilolo[idx+i][idx+j];
         Tlololo[i][j] = Ulololo[idx+i][idx+j];
      }

   CPU_dbl8_upper_inverse
      (szt,Thihihi,   Tlohihi,   Thilohi,   Tlolohi,
           Thihilo,   Tlohilo,   Thilolo,   Tlololo,
        invThihihi,invTlohihi,invThilohi,invTlolohi,
        invThihilo,invTlohilo,invThilolo,invTlololo,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Uhihihi[idx+i][idx+j] = invThihihi[i][j];
         Ulohihi[idx+i][idx+j] = invTlohihi[i][j];
         Uhilohi[idx+i][idx+j] = invThilohi[i][j];
         Ulolohi[idx+i][idx+j] = invTlolohi[i][j];
         Uhihilo[idx+i][idx+j] = invThihilo[i][j];
         Ulohilo[idx+i][idx+j] = invTlohilo[i][j];
         Uhilolo[idx+i][idx+j] = invThilolo[i][j];
         Ulololo[idx+i][idx+j] = invTlololo[i][j];
      }

   for(int i=0; i<szt; i++)
   {
      xhihihi[idx+i] = 0.0; xlohihi[idx+i] = 0.0;
      xhilohi[idx+i] = 0.0; xlolohi[idx+i] = 0.0;
      xhihilo[idx+i] = 0.0; xlohilo[idx+i] = 0.0;
      xhilolo[idx+i] = 0.0; xlololo[idx+i] = 0.0;

      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         odf_mul(invThihihi[i][j],invTlohihi[i][j],
                 invThilohi[i][j],invTlolohi[i][j],
                 invThihilo[i][j],invTlohilo[i][j],
                 invThilolo[i][j],invTlololo[i][j],
                    bhihihi[idx+j],  blohihi[idx+j],
                    bhilohi[idx+j],  blolohi[idx+j],
                    bhihilo[idx+j],  blohilo[idx+j],
                    bhilolo[idx+j],  blololo[idx+j],
                 &acchihihi,      &acclohihi,    &acchilohi,    &acclolohi,
                 &acchihilo,      &acclohilo,    &acchilolo,    &acclololo);
         odf_inc(&xhihihi[idx+i],&xlohihi[idx+i],
                 &xhilohi[idx+i],&xlolohi[idx+i],
                 &xhihilo[idx+i],&xlohilo[idx+i],
                 &xhilolo[idx+i],&xlololo[idx+i],
                 acchihihi,      acclohihi,      acchilohi,      acclolohi,
                 acchihilo,      acclohilo,      acchilolo,      acclololo);
      }
   }
   double *wbhihihi = new double[szt];    // work space for b
   double *wblohihi = new double[szt];
   double *wbhilohi = new double[szt];
   double *wblolohi = new double[szt];
   double *wbhihilo = new double[szt];
   double *wblohilo = new double[szt];
   double *wbhilolo = new double[szt];
   double *wblololo = new double[szt];
   double **wThihihi = new double*[szt];  // work space for a tile
   double **wTlohihi = new double*[szt];
   double **wThilohi = new double*[szt];
   double **wTlolohi = new double*[szt];
   double **wThihilo = new double*[szt];
   double **wTlohilo = new double*[szt];
   double **wThilolo = new double*[szt];
   double **wTlololo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      wThihihi[i] = new double[szt];
      wTlohihi[i] = new double[szt];
      wThilohi[i] = new double[szt];
      wTlolohi[i] = new double[szt];
      wThihilo[i] = new double[szt];
      wTlohilo[i] = new double[szt];
      wThilolo[i] = new double[szt];
      wTlololo[i] = new double[szt];
   }
   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt;  // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Thihihi[i][j] = Uhihihi[idx+i][idx+j];
            Tlohihi[i][j] = Ulohihi[idx+i][idx+j];
            Thilohi[i][j] = Uhilohi[idx+i][idx+j];
            Tlolohi[i][j] = Ulolohi[idx+i][idx+j];
            Thihilo[i][j] = Uhihilo[idx+i][idx+j];
            Tlohilo[i][j] = Ulohilo[idx+i][idx+j];
            Thilolo[i][j] = Uhilolo[idx+i][idx+j];
            Tlololo[i][j] = Ulololo[idx+i][idx+j];
         }

      CPU_dbl8_upper_inverse
         (szt,Thihihi,   Tlohihi,   Thilohi,   Tlolohi,
              Thihilo,   Tlohilo,   Thilolo,   Tlololo,
           invThihihi,invTlohihi,invThilohi,invTlolohi,
           invThihilo,invTlohilo,invThilolo,invTlololo,&timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhihihi[idx+i][idx+j] = invThihihi[i][j];
            Ulohihi[idx+i][idx+j] = invTlohihi[i][j];
            Uhilohi[idx+i][idx+j] = invThilohi[i][j];
            Ulolohi[idx+i][idx+j] = invTlolohi[i][j];
            Uhihilo[idx+i][idx+j] = invThihilo[i][j];
            Ulohilo[idx+i][idx+j] = invTlohilo[i][j];
            Uhilolo[idx+i][idx+j] = invThilolo[i][j];
            Ulololo[idx+i][idx+j] = invTlololo[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++)  // load the work space
         {
            wbhihihi[i] = bhihihi[rowidx+i];
            wblohihi[i] = blohihi[rowidx+i];
            wbhilohi[i] = bhilohi[rowidx+i];
            wblolohi[i] = blolohi[rowidx+i];
            wbhihilo[i] = bhihilo[rowidx+i];
            wblohilo[i] = blohilo[rowidx+i];
            wbhilolo[i] = bhilolo[rowidx+i];
            wblololo[i] = blololo[rowidx+i];

            for(int j=0; j<szt; j++)
            {
               wThihihi[i][j] = Uhihihi[rowidx+i][idx+szt+j];
               wTlohihi[i][j] = Ulohihi[rowidx+i][idx+szt+j];
               wThilohi[i][j] = Uhilohi[rowidx+i][idx+szt+j];
               wTlolohi[i][j] = Ulolohi[rowidx+i][idx+szt+j];
               wThihilo[i][j] = Uhihilo[rowidx+i][idx+szt+j];
               wTlohilo[i][j] = Ulohilo[rowidx+i][idx+szt+j];
               wThilolo[i][j] = Uhilolo[rowidx+i][idx+szt+j];
               wTlololo[i][j] = Ulololo[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodhihihi = 0.0; prodlohihi = 0.0;
            prodhilohi = 0.0; prodlolohi = 0.0;
            prodhihilo = 0.0; prodlohilo = 0.0;
            prodhilolo = 0.0; prodlololo = 0.0;

            for(int j=0; j<szt; j++)
            {  // prod = prod + wT[i][j]*x[idx+szt+j];
               odf_mul(wThihihi[i][j],wTlohihi[i][j],
                       wThilohi[i][j],wTlolohi[i][j],
                       wThihilo[i][j],wTlohilo[i][j],
                       wThilolo[i][j],wTlololo[i][j],
                        xhihihi[idx+szt+j],xlohihi[idx+szt+j],
                        xhilohi[idx+szt+j],xlolohi[idx+szt+j],
                        xhihilo[idx+szt+j],xlohilo[idx+szt+j],
                        xhilolo[idx+szt+j],xlololo[idx+szt+j],
                     &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                     &acchihilo,&acclohilo,&acchilolo,&acclololo);
               odf_inc(&prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
                       &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo,
                         acchihihi,  acclohihi,  acchilohi,  acclolohi,
                         acchihilo,  acclohilo,  acchilolo,  acclololo);
            }
            // wb[i] = wb[i] - prod;
            odf_dec(&wbhihihi[i],&wblohihi[i],&wbhilohi[i],&wblolohi[i],
                    &wbhihilo[i],&wblohilo[i],&wbhilolo[i],&wblololo[i],
                   prodhihihi,  prodlohihi,  prodhilohi,  prodlolohi,
                   prodhihilo,  prodlohilo,  prodhilolo,  prodlololo);
         }
         for(int i=0; i<szt; i++) // store wb into b for next update
         {
            bhihihi[rowidx+i] = wbhihihi[i];
            blohihi[rowidx+i] = wblohihi[i];
            bhilohi[rowidx+i] = wbhilohi[i];
            blolohi[rowidx+i] = wblolohi[i];
            bhihilo[rowidx+i] = wbhihilo[i];
            blohilo[rowidx+i] = wblohilo[i];
            bhilolo[rowidx+i] = wbhilolo[i];
            blololo[rowidx+i] = wblololo[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodhihihi = 0.0; prodlohihi = 0.0;
         prodhilohi = 0.0; prodlolohi = 0.0;
         prodhihilo = 0.0; prodlohilo = 0.0;
         prodhilolo = 0.0; prodlololo = 0.0;

         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            odf_mul(invThihihi[i][j],invTlohihi[i][j],
                    invThilohi[i][j],invTlolohi[i][j],
                    invThihilo[i][j],invTlohilo[i][j],
                    invThilolo[i][j],invTlololo[i][j],
                       bhihihi[idx+j],  blohihi[idx+j],
                       bhilohi[idx+j],  blolohi[idx+j],
                       bhihilo[idx+j],  blohilo[idx+j],
                       bhilolo[idx+j],  blololo[idx+j],
                    &acchihihi,      &acclohihi,    &acchilohi,   &acclolohi,
                    &acchihilo,      &acclohilo,    &acchilolo,   &acclololo);
            odf_inc(&prodhihihi,&prodlohihi,&prodhilohi,&prodlolohi,
                    &prodhihilo,&prodlohilo,&prodhilolo,&prodlololo,
                      acchihihi,  acclohihi,  acchilohi,  acclolohi,
                      acchihilo,  acclohilo,  acchilolo,  acclololo);
         }
         wbhihihi[i] = prodhihihi; wblohihi[i] = prodlohihi;
         wbhilohi[i] = prodhilohi; wblolohi[i] = prodlolohi;
         wbhihilo[i] = prodhihilo; wblohilo[i] = prodlohilo;
         wbhilolo[i] = prodhilolo; wblololo[i] = prodlololo;
      }
      for(int i=0; i<szt; i++)
      {
         xhihihi[idx+i] = wbhihihi[i]; xlohihi[idx+i] = wblohihi[i];
         xhilohi[idx+i] = wbhilohi[i]; xlolohi[idx+i] = wblolohi[i];
         xhihilo[idx+i] = wbhihilo[i]; xlohilo[idx+i] = wblohilo[i];
         xhilolo[idx+i] = wbhilolo[i]; xlololo[idx+i] = wblololo[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Thihihi[i]); free(Tlohihi[i]); free(Thilohi[i]); free(Tlolohi[i]);
      free(Thihilo[i]); free(Tlohilo[i]); free(Thilolo[i]); free(Tlololo[i]);
      free(invThihihi[i]); free(invTlohihi[i]);
      free(invThilohi[i]); free(invTlolohi[i]);
      free(invThihilo[i]); free(invTlohilo[i]);
      free(invThilolo[i]); free(invTlololo[i]);
      free(wThihihi[i]); free(wTlohihi[i]);
      free(wThilohi[i]); free(wTlolohi[i]);
      free(wThihilo[i]); free(wTlohilo[i]);
      free(wThilolo[i]); free(wTlololo[i]);
   }
   free(Thihihi); free(Tlohihi); free(Thilohi); free(Tlolohi);
   free(Thihilo); free(Tlohilo); free(Thilolo); free(Tlololo);
   free(invThihihi); free(invTlohihi); free(invThilohi); free(invTlolohi);
   free(invThihilo); free(invTlohilo); free(invThilolo); free(invTlololo);
   free(wbhihihi); free(wblohihi); free(wbhilohi); free(wblolohi);
   free(wbhihilo); free(wblohilo); free(wbhilolo); free(wblololo);
   free(wThihihi); free(wTlohihi); free(wThilohi); free(wTlolohi);
   free(wThihilo); free(wTlohilo); free(wThilolo); free(wTlololo);
}

void CPU_cmplx8_upper_tiled_solver
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
   double *lapsec )
{
   double acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi;
   double acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo;
   double acc2hihihi,acc2lohihi,acc2hilohi,acc2lolohi;
   double acc2hihilo,acc2lohilo,acc2hilolo,acc2lololo;
   double acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi;
   double acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo;
   double acc4hihihi,acc4lohihi,acc4hilohi,acc4lolohi;
   double acc4hihilo,acc4lohilo,acc4hilolo,acc4lololo;
   double prodrehihihi,prodrelohihi,prodrehilohi,prodrelolohi;
   double prodrehihilo,prodrelohilo,prodrehilolo,prodrelololo;
   double prodimhihihi,prodimlohihi,prodimhilohi,prodimlolohi;
   double prodimhihilo,prodimlohilo,prodimhilolo,prodimlololo;
   double **Trehihihi = new double*[szt];
   double **Trelohihi = new double*[szt];
   double **Trehilohi = new double*[szt];
   double **Trelolohi = new double*[szt];
   double **Trehihilo = new double*[szt];
   double **Trelohilo = new double*[szt];
   double **Trehilolo = new double*[szt];
   double **Trelololo = new double*[szt];
   double **Timhihihi = new double*[szt];
   double **Timlohihi = new double*[szt];
   double **Timhilohi = new double*[szt];
   double **Timlolohi = new double*[szt];
   double **Timhihilo = new double*[szt];
   double **Timlohilo = new double*[szt];
   double **Timhilolo = new double*[szt];
   double **Timlololo = new double*[szt];
   double **invTrehihihi = new double*[szt];
   double **invTrelohihi = new double*[szt];
   double **invTrehilohi = new double*[szt];
   double **invTrelolohi = new double*[szt];
   double **invTrehihilo = new double*[szt];
   double **invTrelohilo = new double*[szt];
   double **invTrehilolo = new double*[szt];
   double **invTrelololo = new double*[szt];
   double **invTimhihihi = new double*[szt];
   double **invTimlohihi = new double*[szt];
   double **invTimhilohi = new double*[szt];
   double **invTimlolohi = new double*[szt];
   double **invTimhihilo = new double*[szt];
   double **invTimlohilo = new double*[szt];
   double **invTimhilolo = new double*[szt];
   double **invTimlololo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Trehihihi[i] = new double[szt];
      Trelohihi[i] = new double[szt];
      Trehilohi[i] = new double[szt];
      Trelolohi[i] = new double[szt];
      Trehihilo[i] = new double[szt];
      Trelohilo[i] = new double[szt];
      Trehilolo[i] = new double[szt];
      Trelololo[i] = new double[szt];
      Timhihihi[i] = new double[szt];
      Timlohihi[i] = new double[szt];
      Timhilohi[i] = new double[szt];
      Timlolohi[i] = new double[szt];
      Timhihilo[i] = new double[szt];
      Timlohilo[i] = new double[szt];
      Timhilolo[i] = new double[szt];
      Timlololo[i] = new double[szt];
      invTrehihihi[i] = new double[szt];
      invTrelohihi[i] = new double[szt];
      invTrehilohi[i] = new double[szt];
      invTrelolohi[i] = new double[szt];
      invTrehihilo[i] = new double[szt];
      invTrelohilo[i] = new double[szt];
      invTrehilolo[i] = new double[szt];
      invTrelololo[i] = new double[szt];
      invTimhihihi[i] = new double[szt];
      invTimlohihi[i] = new double[szt];
      invTimhilohi[i] = new double[szt];
      invTimlolohi[i] = new double[szt];
      invTimhihilo[i] = new double[szt];
      invTimlohilo[i] = new double[szt];
      invTimhilolo[i] = new double[szt];
      invTimlololo[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Trehihihi[i][j] = Urehihihi[idx+i][idx+j];
         Trelohihi[i][j] = Urelohihi[idx+i][idx+j];
         Trehilohi[i][j] = Urehilohi[idx+i][idx+j];
         Trelolohi[i][j] = Urelolohi[idx+i][idx+j];
         Trehihilo[i][j] = Urehihilo[idx+i][idx+j];
         Trelohilo[i][j] = Urelohilo[idx+i][idx+j];
         Trehilolo[i][j] = Urehilolo[idx+i][idx+j];
         Trelololo[i][j] = Urelololo[idx+i][idx+j];
         Timhihihi[i][j] = Uimhihihi[idx+i][idx+j];
         Timlohihi[i][j] = Uimlohihi[idx+i][idx+j];
         Timhilohi[i][j] = Uimhilohi[idx+i][idx+j];
         Timlolohi[i][j] = Uimlolohi[idx+i][idx+j];
         Timhihilo[i][j] = Uimhihilo[idx+i][idx+j];
         Timlohilo[i][j] = Uimlohilo[idx+i][idx+j];
         Timhilolo[i][j] = Uimhilolo[idx+i][idx+j];
         Timlololo[i][j] = Uimlololo[idx+i][idx+j];
      }

   CPU_cmplx8_upper_inverse
      (szt,Trehihihi,Trelohihi,Trehilohi,Trelolohi,
           Trehihilo,Trelohilo,Trehilolo,Trelololo,
           Timhihihi,Timlohihi,Timhilohi,Timlolohi,
           Timhihilo,Timlohilo,Timhilolo,Timlololo,
           invTrehihihi,invTrelohihi,invTrehilohi,invTrelolohi,
           invTrehihilo,invTrelohilo,invTrehilolo,invTrelololo,
           invTimhihihi,invTimlohihi,invTimhilohi,invTimlolohi,
           invTimhihilo,invTimlohilo,invTimhilolo,invTimlololo,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Urehihihi[idx+i][idx+j] = invTrehihihi[i][j];
         Urelohihi[idx+i][idx+j] = invTrelohihi[i][j];
         Urehilohi[idx+i][idx+j] = invTrehilohi[i][j];
         Urelolohi[idx+i][idx+j] = invTrelolohi[i][j];
         Urehihilo[idx+i][idx+j] = invTrehihilo[i][j];
         Urelohilo[idx+i][idx+j] = invTrelohilo[i][j];
         Urehilolo[idx+i][idx+j] = invTrehilolo[i][j];
         Urelololo[idx+i][idx+j] = invTrelololo[i][j];
         Uimhihihi[idx+i][idx+j] = invTimhihihi[i][j];
         Uimlohihi[idx+i][idx+j] = invTimlohihi[i][j];
         Uimhilohi[idx+i][idx+j] = invTimhilohi[i][j];
         Uimlolohi[idx+i][idx+j] = invTimlolohi[i][j];
         Uimhihilo[idx+i][idx+j] = invTimhihilo[i][j];
         Uimlohilo[idx+i][idx+j] = invTimlohilo[i][j];
         Uimhilolo[idx+i][idx+j] = invTimhilolo[i][j];
         Uimlololo[idx+i][idx+j] = invTimlololo[i][j];
      }

   for(int i=0; i<szt; i++)
   {
      xrehihihi[idx+i] = 0.0; xrelohihi[idx+i] = 0.0;
      xrehilohi[idx+i] = 0.0; xrelolohi[idx+i] = 0.0;
      xrehihilo[idx+i] = 0.0; xrelohilo[idx+i] = 0.0;
      xrehilolo[idx+i] = 0.0; xrelololo[idx+i] = 0.0;
      ximhihihi[idx+i] = 0.0; ximlohihi[idx+i] = 0.0;
      ximhilohi[idx+i] = 0.0; ximlolohi[idx+i] = 0.0;
      ximhihilo[idx+i] = 0.0; ximlohilo[idx+i] = 0.0;
      ximhilolo[idx+i] = 0.0; ximlololo[idx+i] = 0.0;

      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         odf_mul(invTrehihihi[i][j],invTrelohihi[i][j],
                 invTrehilohi[i][j],invTrelolohi[i][j],
                 invTrehihilo[i][j],invTrelohilo[i][j],
                 invTrehilolo[i][j],invTrelololo[i][j],
                 brehihihi[idx+j],brelohihi[idx+j],
                 brehilohi[idx+j],brelolohi[idx+j],
                 brehihilo[idx+j],brelohilo[idx+j],
                 brehilolo[idx+j],brelololo[idx+j],
                 &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                 &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
         odf_mul(invTimhihihi[i][j],invTimlohihi[i][j],
                 invTimhilohi[i][j],invTimlolohi[i][j],
                 invTimhihilo[i][j],invTimlohilo[i][j],
                 invTimhilolo[i][j],invTimlololo[i][j],
                 bimhihihi[idx+j],bimlohihi[idx+j],
                 bimhilohi[idx+j],bimlolohi[idx+j],
                 bimhihilo[idx+j],bimlohilo[idx+j],
                 bimhilolo[idx+j],bimlololo[idx+j],
                 &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
                 &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
         odf_mul(invTimhihihi[i][j],invTimlohihi[i][j],
                 invTimhilohi[i][j],invTimlolohi[i][j],
                 invTimhihilo[i][j],invTimlohilo[i][j],
                 invTimhilolo[i][j],invTimlololo[i][j],
                 brehihihi[idx+j],brelohihi[idx+j],
                 brehilohi[idx+j],brelolohi[idx+j],
                 brehihilo[idx+j],brelohilo[idx+j],
                 brehilolo[idx+j],brelololo[idx+j],
                 &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                 &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
         odf_mul(invTrehihihi[i][j],invTrelohihi[i][j],
                 invTrehilohi[i][j],invTrelolohi[i][j],
                 invTrehihilo[i][j],invTrelohilo[i][j],
                 invTrehilolo[i][j],invTrelololo[i][j],
                 bimhihihi[idx+j],bimlohihi[idx+j],
                 bimhilohi[idx+j],bimlolohi[idx+j],
                 bimhihilo[idx+j],bimlohilo[idx+j],
                 bimhilolo[idx+j],bimlololo[idx+j],
                 &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                 &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);
         odf_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                 &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                  acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                  acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
         odf_inc(&xrehihihi[idx+i],&xrelohihi[idx+i],
                 &xrehilohi[idx+i],&xrelolohi[idx+i],
                 &xrehihilo[idx+i],&xrelohilo[idx+i],
                 &xrehilolo[idx+i],&xrelololo[idx+i],
                 acc1hihihi,acc1lohihi,acc1hilohi,acc1lolohi,
                 acc1hihilo,acc1lohilo,acc1hilolo,acc1lololo);
         odf_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                 &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
                  acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
                  acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
         odf_inc(&ximhihihi[idx+i],&ximlohihi[idx+i],
                 &ximhilohi[idx+i],&ximlolohi[idx+i],
                 &ximhihilo[idx+i],&ximlohilo[idx+i],
                 &ximhilolo[idx+i],&ximlololo[idx+i],
                 acc3hihihi,acc3lohihi,acc3hilohi,acc3lolohi,
                 acc3hihilo,acc3lohilo,acc3hilolo,acc3lololo);
      }
   }
   double *wbrehihihi = new double[szt];    // work space for b
   double *wbrelohihi = new double[szt]; 
   double *wbrehilohi = new double[szt];
   double *wbrelolohi = new double[szt];
   double *wbrehihilo = new double[szt];
   double *wbrelohilo = new double[szt]; 
   double *wbrehilolo = new double[szt];
   double *wbrelololo = new double[szt];
   double *wbimhihihi = new double[szt];
   double *wbimlohihi = new double[szt];
   double *wbimhilohi = new double[szt];
   double *wbimlolohi = new double[szt];
   double *wbimhihilo = new double[szt];
   double *wbimlohilo = new double[szt];
   double *wbimhilolo = new double[szt];
   double *wbimlololo = new double[szt];
   double **wTrehihihi = new double*[szt];  // work space for a tile
   double **wTrelohihi = new double*[szt];
   double **wTrehilohi = new double*[szt];
   double **wTrelolohi = new double*[szt];
   double **wTrehihilo = new double*[szt];
   double **wTrelohilo = new double*[szt];
   double **wTrehilolo = new double*[szt];
   double **wTrelololo = new double*[szt];
   double **wTimhihihi = new double*[szt];
   double **wTimlohihi = new double*[szt];
   double **wTimhilohi = new double*[szt]; 
   double **wTimlolohi = new double*[szt];
   double **wTimhihilo = new double*[szt];
   double **wTimlohilo = new double*[szt];
   double **wTimhilolo = new double*[szt]; 
   double **wTimlololo = new double*[szt];
 
   for(int i=0; i<szt; i++)
   {
      wTrehihihi[i] = new double[szt];
      wTrelohihi[i] = new double[szt];
      wTrehilohi[i] = new double[szt];
      wTrelolohi[i] = new double[szt];
      wTrehihilo[i] = new double[szt];
      wTrelohilo[i] = new double[szt];
      wTrehilolo[i] = new double[szt];
      wTrelololo[i] = new double[szt];
      wTimhihihi[i] = new double[szt];
      wTimlohihi[i] = new double[szt];
      wTimhilohi[i] = new double[szt];
      wTimlolohi[i] = new double[szt];
      wTimhihilo[i] = new double[szt];
      wTimlohilo[i] = new double[szt];
      wTimhilolo[i] = new double[szt];
      wTimlololo[i] = new double[szt];
   }
   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt;  // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Trehihihi[i][j] = Urehihihi[idx+i][idx+j];
            Trelohihi[i][j] = Urelohihi[idx+i][idx+j];
            Trehilohi[i][j] = Urehilohi[idx+i][idx+j];
            Trelolohi[i][j] = Urelolohi[idx+i][idx+j];
            Trehihilo[i][j] = Urehihilo[idx+i][idx+j];
            Trelohilo[i][j] = Urelohilo[idx+i][idx+j];
            Trehilolo[i][j] = Urehilolo[idx+i][idx+j];
            Trelololo[i][j] = Urelololo[idx+i][idx+j];
            Timhihihi[i][j] = Uimhihihi[idx+i][idx+j];
            Timlohihi[i][j] = Uimlohihi[idx+i][idx+j];
            Timhilohi[i][j] = Uimhilohi[idx+i][idx+j];
            Timlolohi[i][j] = Uimlolohi[idx+i][idx+j];
            Timhihilo[i][j] = Uimhihilo[idx+i][idx+j];
            Timlohilo[i][j] = Uimlohilo[idx+i][idx+j];
            Timhilolo[i][j] = Uimhilolo[idx+i][idx+j];
            Timlololo[i][j] = Uimlololo[idx+i][idx+j];
         }

      CPU_cmplx8_upper_inverse
         (szt,Trehihihi,Trelohihi,Trehilohi,Trelolohi,
              Trehihilo,Trelohilo,Trehilolo,Trelololo,
              Timhihihi,Timlohihi,Timhilohi,Timlolohi,
              Timhihilo,Timlohilo,Timhilolo,Timlololo,
          invTrehihihi,invTrelohihi,invTrehilohi,invTrelolohi,
          invTrehihilo,invTrelohilo,invTrehilolo,invTrelololo,
          invTimhihihi,invTimlohihi,invTimhilohi,invTimlolohi,
          invTimhihilo,invTimlohilo,invTimhilolo,invTimlololo,&timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Urehihihi[idx+i][idx+j] = invTrehihihi[i][j];
            Urelohihi[idx+i][idx+j] = invTrelohihi[i][j];
            Urehilohi[idx+i][idx+j] = invTrehilohi[i][j];
            Urelolohi[idx+i][idx+j] = invTrelolohi[i][j];
            Urehihilo[idx+i][idx+j] = invTrehihilo[i][j];
            Urelohilo[idx+i][idx+j] = invTrelohilo[i][j];
            Urehilolo[idx+i][idx+j] = invTrehilolo[i][j];
            Urelololo[idx+i][idx+j] = invTrelololo[i][j];
            Uimhihihi[idx+i][idx+j] = invTimhihihi[i][j];
            Uimlohihi[idx+i][idx+j] = invTimlohihi[i][j];
            Uimhilohi[idx+i][idx+j] = invTimhilohi[i][j];
            Uimlolohi[idx+i][idx+j] = invTimlolohi[i][j];
            Uimhihilo[idx+i][idx+j] = invTimhihilo[i][j];
            Uimlohilo[idx+i][idx+j] = invTimlohilo[i][j];
            Uimhilolo[idx+i][idx+j] = invTimhilolo[i][j];
            Uimlololo[idx+i][idx+j] = invTimlololo[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++)  // load the work space
         {
            wbrehihihi[i] = brehihihi[rowidx+i];
            wbrelohihi[i] = brelohihi[rowidx+i];
            wbrehilohi[i] = brehilohi[rowidx+i];
            wbrelolohi[i] = brelolohi[rowidx+i];
            wbrehihilo[i] = brehihilo[rowidx+i];
            wbrelohilo[i] = brelohilo[rowidx+i];
            wbrehilolo[i] = brehilolo[rowidx+i];
            wbrelololo[i] = brelololo[rowidx+i];
            wbimhihihi[i] = bimhihihi[rowidx+i];
            wbimlohihi[i] = bimlohihi[rowidx+i];
            wbimhilohi[i] = bimhilohi[rowidx+i];
            wbimlolohi[i] = bimlolohi[rowidx+i];
            wbimhihilo[i] = bimhihilo[rowidx+i];
            wbimlohilo[i] = bimlohilo[rowidx+i];
            wbimhilolo[i] = bimhilolo[rowidx+i];
            wbimlololo[i] = bimlololo[rowidx+i];

            for(int j=0; j<szt; j++)
            {
               wTrehihihi[i][j] = Urehihihi[rowidx+i][idx+szt+j];
               wTrehilohi[i][j] = Urehilohi[rowidx+i][idx+szt+j];
               wTrelohihi[i][j] = Urelohihi[rowidx+i][idx+szt+j];
               wTrelolohi[i][j] = Urelolohi[rowidx+i][idx+szt+j];
               wTrehihilo[i][j] = Urehihilo[rowidx+i][idx+szt+j];
               wTrehilolo[i][j] = Urehilolo[rowidx+i][idx+szt+j];
               wTrelohilo[i][j] = Urelohilo[rowidx+i][idx+szt+j];
               wTrelololo[i][j] = Urelololo[rowidx+i][idx+szt+j];
               wTimhihihi[i][j] = Uimhihihi[rowidx+i][idx+szt+j];
               wTimhilohi[i][j] = Uimhilohi[rowidx+i][idx+szt+j];
               wTimlohihi[i][j] = Uimlohihi[rowidx+i][idx+szt+j];
               wTimlolohi[i][j] = Uimlolohi[rowidx+i][idx+szt+j];
               wTimhihilo[i][j] = Uimhihilo[rowidx+i][idx+szt+j];
               wTimhilolo[i][j] = Uimhilolo[rowidx+i][idx+szt+j];
               wTimlohilo[i][j] = Uimlohilo[rowidx+i][idx+szt+j];
               wTimlololo[i][j] = Uimlololo[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodrehihihi = 0.0; prodrelohihi = 0.0;
            prodrehilohi = 0.0; prodrelolohi = 0.0;
            prodrehihilo = 0.0; prodrelohilo = 0.0;
            prodrehilolo = 0.0; prodrelololo = 0.0;
            prodimhihihi = 0.0; prodimlohihi = 0.0;
            prodimhilohi = 0.0; prodimlolohi = 0.0;
            prodimhihilo = 0.0; prodimlohilo = 0.0;
            prodimhilolo = 0.0; prodimlololo = 0.0;

            for(int j=0; j<szt; j++) // prod = prod + wT[i][j]*x[idx+szt+j];
            {
               odf_mul(wTrehihihi[i][j],wTrelohihi[i][j],
                       wTrehilohi[i][j],wTrelolohi[i][j],
                       wTrehihilo[i][j],wTrelohilo[i][j],
                       wTrehilolo[i][j],wTrelololo[i][j],
                       xrehihihi[idx+szt+j],xrelohihi[idx+szt+j],
                       xrehilohi[idx+szt+j],xrelolohi[idx+szt+j],
                       xrehihilo[idx+szt+j],xrelohilo[idx+szt+j],
                       xrehilolo[idx+szt+j],xrelololo[idx+szt+j],
                       &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                       &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
               odf_mul(wTimhihihi[i][j],wTimlohihi[i][j],
                       wTimhilohi[i][j],wTimlolohi[i][j],
                       wTimhihilo[i][j],wTimlohilo[i][j],
                       wTimhilolo[i][j],wTimlololo[i][j],
                       ximhihihi[idx+szt+j],ximlohihi[idx+szt+j],
                       ximhilohi[idx+szt+j],ximlolohi[idx+szt+j],
                       ximhihilo[idx+szt+j],ximlohilo[idx+szt+j],
                       ximhilolo[idx+szt+j],ximlololo[idx+szt+j],
                       &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
                       &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
               odf_mul(wTimhihihi[i][j],wTimlohihi[i][j],
                       wTimhilohi[i][j],wTimlolohi[i][j],
                       wTimhihilo[i][j],wTimlohilo[i][j],
                       wTimhilolo[i][j],wTimlololo[i][j],
                       xrehihihi[idx+szt+j],xrelohihi[idx+szt+j],
                       xrehilohi[idx+szt+j],xrelolohi[idx+szt+j],
                       xrehihilo[idx+szt+j],xrelohilo[idx+szt+j],
                       xrehilolo[idx+szt+j],xrelololo[idx+szt+j],
                       &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                       &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
               odf_mul(wTrehihihi[i][j],wTrelohihi[i][j],
                       wTrehilohi[i][j],wTrelolohi[i][j],
                       wTrehihilo[i][j],wTrelohilo[i][j],
                       wTrehilolo[i][j],wTrelololo[i][j],
                       ximhihihi[idx+szt+j],ximlohihi[idx+szt+j],
                       ximhilohi[idx+szt+j],ximlolohi[idx+szt+j],
                       ximhihilo[idx+szt+j],ximlohilo[idx+szt+j],
                       ximhilolo[idx+szt+j],ximlololo[idx+szt+j],
                       &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                       &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);
               odf_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                       &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                        acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                        acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
               odf_inc(&prodrehihihi,&prodrelohihi,
                       &prodrehilohi,&prodrelolohi,
                       &prodrehihilo,&prodrelohilo,
                       &prodrehilolo,&prodrelololo,
                          acc1hihihi,   acc1lohihi, acc1hilohi, acc1lolohi,
                          acc1hihilo,   acc1lohilo, acc1hilolo, acc1lololo);
               odf_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                       &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
                        acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
                        acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
               odf_inc(&prodimhihihi,&prodimlohihi,
                       &prodimhilohi,&prodimlolohi,
                       &prodimhihilo,&prodimlohilo,
                       &prodimhilolo,&prodimlololo,
                          acc3hihihi,   acc3lohihi, acc3hilohi, acc3lolohi,
                          acc3hihilo,   acc3lohilo, acc3hilolo, acc3lololo);
            }
            // wb[i] = wb[i] - prod;
            odf_dec(&wbrehihihi[i],&wbrelohihi[i],
                    &wbrehilohi[i],&wbrelolohi[i],
                    &wbrehihilo[i],&wbrelohilo[i],
                    &wbrehilolo[i],&wbrelololo[i],
                   prodrehihihi,  prodrelohihi,  prodrehilohi,  prodrelolohi,
                   prodrehihilo,  prodrelohilo,  prodrehilolo,  prodrelololo);
            odf_dec(&wbimhihihi[i],&wbimlohihi[i],
                    &wbimhilohi[i],&wbimlolohi[i],
                    &wbimhihilo[i],&wbimlohilo[i],
                    &wbimhilolo[i],&wbimlololo[i],
                   prodimhihihi,  prodimlohihi,  prodimhilohi,  prodimlolohi,
                   prodimhihilo,  prodimlohilo,  prodimhilolo,  prodimlololo);
         }
         for(int i=0; i<szt; i++) // store wb into b for next update
         {
            brehihihi[rowidx+i] = wbrehihihi[i];
            brelohihi[rowidx+i] = wbrelohihi[i];
            brehilohi[rowidx+i] = wbrehilohi[i];
            brelolohi[rowidx+i] = wbrelolohi[i];
            brehihilo[rowidx+i] = wbrehihilo[i];
            brelohilo[rowidx+i] = wbrelohilo[i];
            brehilolo[rowidx+i] = wbrehilolo[i];
            brelololo[rowidx+i] = wbrelololo[i];
            bimhihihi[rowidx+i] = wbimhihihi[i];
            bimlohihi[rowidx+i] = wbimlohihi[i];
            bimhilohi[rowidx+i] = wbimhilohi[i];
            bimlolohi[rowidx+i] = wbimlolohi[i];
            bimhihilo[rowidx+i] = wbimhihilo[i];
            bimlohilo[rowidx+i] = wbimlohilo[i];
            bimhilolo[rowidx+i] = wbimhilolo[i];
            bimlololo[rowidx+i] = wbimlololo[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodrehihihi = 0.0; prodrelohihi = 0.0;
         prodrehilohi = 0.0; prodrelolohi = 0.0;
         prodrehihilo = 0.0; prodrelohilo = 0.0;
         prodrehilolo = 0.0; prodrelololo = 0.0;
         prodimhihihi = 0.0; prodimlohihi = 0.0;
         prodimhilohi = 0.0; prodimlolohi = 0.0;
         prodimhihilo = 0.0; prodimlohilo = 0.0;
         prodimhilolo = 0.0; prodimlololo = 0.0;

         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            odf_mul(invTrehihihi[i][j],invTrelohihi[i][j],
                    invTrehilohi[i][j],invTrelolohi[i][j],
                    invTrehihilo[i][j],invTrelohilo[i][j],
                    invTrehilolo[i][j],invTrelololo[i][j],
                    brehihihi[idx+j],brelohihi[idx+j],
                    brehilohi[idx+j],brelolohi[idx+j],
                    brehihilo[idx+j],brelohilo[idx+j],
                    brehilolo[idx+j],brelololo[idx+j],
                    &acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                    &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo);
            odf_mul(invTimhihihi[i][j],invTimlohihi[i][j],
                    invTimhilohi[i][j],invTimlolohi[i][j],
                    invTimhihilo[i][j],invTimlohilo[i][j],
                    invTimhilolo[i][j],invTimlololo[i][j],
                    bimhihihi[idx+j],bimlohihi[idx+j],
                    bimhilohi[idx+j],bimlolohi[idx+j],
                    bimhihilo[idx+j],bimlohilo[idx+j],
                    bimhilolo[idx+j],bimlololo[idx+j],
                    &acc2hihihi,&acc2lohihi,&acc2hilohi,&acc2lolohi,
                    &acc2hihilo,&acc2lohilo,&acc2hilolo,&acc2lololo);
            odf_mul(invTimhihihi[i][j],invTimlohihi[i][j],
                    invTimhilohi[i][j],invTimlolohi[i][j],
                    invTimhihilo[i][j],invTimlohilo[i][j],
                    invTimhilolo[i][j],invTimlololo[i][j],
                    brehihihi[idx+j],brelohihi[idx+j],
                    brehilohi[idx+j],brelolohi[idx+j],
                    brehihilo[idx+j],brelohilo[idx+j],
                    brehilolo[idx+j],brelololo[idx+j],
                    &acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                    &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo);
            odf_mul(invTrehihihi[i][j],invTrelohihi[i][j],
                    invTrehilohi[i][j],invTrelolohi[i][j],
                    invTrehihilo[i][j],invTrelohilo[i][j],
                    invTrehilolo[i][j],invTrelololo[i][j],
                    bimhihihi[idx+j],bimlohihi[idx+j],
                    bimhilohi[idx+j],bimlolohi[idx+j],
                    bimhihilo[idx+j],bimlohilo[idx+j],
                    bimhilolo[idx+j],bimlololo[idx+j],
                    &acc4hihihi,&acc4lohihi,&acc4hilohi,&acc4lolohi,
                    &acc4hihilo,&acc4lohilo,&acc4hilolo,&acc4lololo);
            odf_dec(&acc1hihihi,&acc1lohihi,&acc1hilohi,&acc1lolohi,
                    &acc1hihilo,&acc1lohilo,&acc1hilolo,&acc1lololo,
                     acc2hihihi, acc2lohihi, acc2hilohi, acc2lolohi,
                     acc2hihilo, acc2lohilo, acc2hilolo, acc2lololo);
            odf_inc(&prodrehihihi,&prodrelohihi,&prodrehilohi,&prodrelolohi,
                    &prodrehihilo,&prodrelohilo,&prodrehilolo,&prodrelololo,
                       acc1hihihi,   acc1lohihi,   acc1hilohi,   acc1lolohi,
                       acc1hihilo,   acc1lohilo,   acc1hilolo,   acc1lololo);
            odf_inc(&acc3hihihi,&acc3lohihi,&acc3hilohi,&acc3lolohi,
                    &acc3hihilo,&acc3lohilo,&acc3hilolo,&acc3lololo,
                     acc4hihihi, acc4lohihi, acc4hilohi, acc4lolohi,
                     acc4hihilo, acc4lohilo, acc4hilolo, acc4lololo);
            odf_inc(&prodimhihihi,&prodimlohihi,&prodimhilohi,&prodimlolohi,
                    &prodimhihilo,&prodimlohilo,&prodimhilolo,&prodimlololo,
                       acc3hihihi,   acc3lohihi,   acc3hilohi,   acc3lolohi,
                       acc3hihilo,   acc3lohilo,   acc3hilolo,   acc3lololo);
         }
         wbrehihihi[i] = prodrehihihi; wbrelohihi[i] = prodrelohihi;
         wbrehilohi[i] = prodrehilohi; wbrelolohi[i] = prodrelolohi;
         wbrehihilo[i] = prodrehihilo; wbrelohilo[i] = prodrelohilo;
         wbrehilolo[i] = prodrehilolo; wbrelololo[i] = prodrelololo;
         wbimhihihi[i] = prodimhihihi; wbimlohihi[i] = prodimlohihi;
         wbimhilohi[i] = prodimhilohi; wbimlolohi[i] = prodimlolohi;
         wbimhihilo[i] = prodimhihilo; wbimlohilo[i] = prodimlohilo;
         wbimhilolo[i] = prodimhilolo; wbimlololo[i] = prodimlololo;
      }
      for(int i=0; i<szt; i++)
      {
         xrehihihi[idx+i] = wbrehihihi[i]; xrelohihi[idx+i] = wbrelohihi[i];
         xrehilohi[idx+i] = wbrehilohi[i]; xrelolohi[idx+i] = wbrelolohi[i];
         xrehihilo[idx+i] = wbrehihilo[i]; xrelohilo[idx+i] = wbrelohilo[i];
         xrehilolo[idx+i] = wbrehilolo[i]; xrelololo[idx+i] = wbrelololo[i];
         ximhihihi[idx+i] = wbimhihihi[i]; ximlohihi[idx+i] = wbimlohihi[i];
         ximhilohi[idx+i] = wbimhilohi[i]; ximlolohi[idx+i] = wbimlolohi[i];
         ximhihilo[idx+i] = wbimhihilo[i]; ximlohilo[idx+i] = wbimlohilo[i];
         ximhilolo[idx+i] = wbimhilolo[i]; ximlololo[idx+i] = wbimlololo[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Trehihihi[i]); free(Trelohihi[i]);
      free(Trehilohi[i]); free(Trelolohi[i]);
      free(Trehihilo[i]); free(Trelohilo[i]);
      free(Trehilolo[i]); free(Trelololo[i]);
      free(Timhihihi[i]); free(Timlohihi[i]);
      free(Timhilohi[i]); free(Timlolohi[i]);
      free(Timhihilo[i]); free(Timlohilo[i]);
      free(Timhilolo[i]); free(Timlololo[i]);
      free(invTrehihihi[i]); free(invTrelohihi[i]);
      free(invTrehilohi[i]); free(invTrelolohi[i]);
      free(invTrehihilo[i]); free(invTrelohilo[i]);
      free(invTrehilolo[i]); free(invTrelololo[i]);
      free(invTimhihihi[i]); free(invTimlohihi[i]);
      free(invTimhilohi[i]); free(invTimlolohi[i]);
      free(invTimhihilo[i]); free(invTimlohilo[i]);
      free(invTimhilolo[i]); free(invTimlololo[i]);
      free(wTrehihihi[i]); free(wTrelohihi[i]);
      free(wTrehilohi[i]); free(wTrelolohi[i]);
      free(wTrehihilo[i]); free(wTrelohilo[i]);
      free(wTrehilolo[i]); free(wTrelololo[i]);
      free(wTimhihihi[i]); free(wTimlohihi[i]);
      free(wTimhilohi[i]); free(wTimlolohi[i]);
      free(wTimhihilo[i]); free(wTimlohilo[i]);
      free(wTimhilolo[i]); free(wTimlololo[i]);
   }
   free(Trehihihi); free(Trelohihi); free(Trehilohi); free(Trelolohi);
   free(Trehihilo); free(Trelohilo); free(Trehilolo); free(Trelololo);
   free(Timhihihi); free(Timlohihi); free(Timhilohi); free(Timlolohi);
   free(Timhihilo); free(Timlohilo); free(Timhilolo); free(Timlololo);
   free(invTrehihihi); free(invTrelohihi); 
   free(invTrehilohi); free(invTrelolohi); 
   free(invTrehihilo); free(invTrelohilo); 
   free(invTrehilolo); free(invTrelololo); 
   free(invTimhihihi); free(invTimlohihi);
   free(invTimhilohi); free(invTimlolohi);
   free(invTimhihilo); free(invTimlohilo);
   free(invTimhilolo); free(invTimlololo);
   free(wbrehihihi); free(wbrelohihi); free(wbrehilohi); free(wbrelolohi);
   free(wbrehihilo); free(wbrelohilo); free(wbrehilolo); free(wbrelololo);
   free(wbimhihihi); free(wbimlohihi); free(wbimhilohi); free(wbimlolohi);
   free(wbimhihilo); free(wbimlohilo); free(wbimhilolo); free(wbimlololo);
   free(wTrehihihi); free(wTrelohihi); free(wTrehilohi); free(wTrelolohi);
   free(wTrehihilo); free(wTrelohilo); free(wTrehilolo); free(wTrelololo);
   free(wTimhihihi); free(wTimlohihi); free(wTimhilohi); free(wTimlolohi);
   free(wTimhihilo); free(wTimlohilo); free(wTimhilolo); free(wTimlololo);
}
