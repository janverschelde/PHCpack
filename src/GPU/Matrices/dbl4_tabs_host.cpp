/* The file dbl4_tabs_host.cpp defines functions specified in
 * the file dbl4_tabs_host.h. */

#include <cstdlib>
#include <ctime>
#include "quad_double_functions.h"
#include "dbl4_factorizations.h"
#include "dbl4_tabs_host.h"

using namespace std;

void CPU_dbl4_backsubs
 ( int dim, double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo )
{
   CPU_dbl4_factors_backward
      (dim,Uhihi,Ulohi,Uhilo,Ulolo,bhihi,blohi,bhilo,blolo,
           xhihi,xlohi,xhilo,xlolo);
}

void CPU_cmplx4_backsubs
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo )
{
   CPU_cmplx4_factors_backward
      (dim,Urehihi,Urelohi,Urehilo,Urelolo,
           Uimhihi,Uimlohi,Uimhilo,Uimlolo,
           brehihi,brelohi,brehilo,brelolo,
           bimhihi,bimlohi,bimhilo,bimlolo,
           xrehihi,xrelohi,xrehilo,xrelolo,
           ximhihi,ximlohi,ximhilo,ximlolo);
}

void CPU_dbl4_upper_inverse
 ( int dim,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double **invUhihi, double **invUlohi, double **invUhilo, double **invUlolo,
   double *lapsec )
{
   double *colhihi = new double[dim];
   double *collohi = new double[dim];
   double *colhilo = new double[dim];
   double *collolo = new double[dim];
   double *rhshihi = new double[dim];
   double *rhslohi = new double[dim];
   double *rhshilo = new double[dim];
   double *rhslolo = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhshihi[i] = 0.0; rhslohi[i] = 0.0;
      rhshilo[i] = 0.0; rhslolo[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhshihi[j] = 1.0;

      CPU_dbl4_backsubs
         (dim,Uhihi,Ulohi,Uhilo,Ulolo,rhshihi,rhslohi,rhshilo,rhslolo,
          colhihi,collohi,colhilo,collolo);

      for(int i=0; i<dim; i++)
      {
         invUhihi[i][j] = colhihi[i]; invUlohi[i][j] = collohi[i];
         invUhilo[i][j] = colhilo[i]; invUlolo[i][j] = collolo[i];
      }
      rhshihi[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhshihi); free(rhslohi); free(rhshilo); free(rhslolo);
   free(colhihi); free(collohi); free(colhilo); free(collolo);
}

void CPU_cmplx4_upper_inverse
 ( int dim,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double **invUrehihi, double **invUrelohi,
   double **invUrehilo, double **invUrelolo,
   double **invUimhihi, double **invUimlohi,
   double **invUimhilo, double **invUimlolo, double *lapsec )
{
   double *colrehihi = new double[dim];
   double *colrelohi = new double[dim];
   double *colrehilo = new double[dim];
   double *colrelolo = new double[dim];
   double *colimhihi = new double[dim];
   double *colimlohi = new double[dim];
   double *colimhilo = new double[dim];
   double *colimlolo = new double[dim];
   double *rhsrehihi = new double[dim];
   double *rhsrelohi = new double[dim];
   double *rhsrehilo = new double[dim];
   double *rhsrelolo = new double[dim];
   double *rhsimhihi = new double[dim];
   double *rhsimlohi = new double[dim];
   double *rhsimhilo = new double[dim];
   double *rhsimlolo = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhsrehihi[i] = 0.0; rhsrelohi[i] = 0.0;
      rhsrehilo[i] = 0.0; rhsrelolo[i] = 0.0;
      rhsimhihi[i] = 0.0; rhsimlohi[i] = 0.0;
      rhsimhilo[i] = 0.0; rhsimlolo[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhsrehihi[j] = 1.0;

      CPU_cmplx4_backsubs
         (dim, Urehihi,  Urelohi,  Urehilo,  Urelolo,
               Uimhihi,  Uimlohi,  Uimhilo,  Uimlolo,
             rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
             rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
             colrehihi,colrelohi,colrehilo,colrelolo,
             colimhihi,colimlohi,colimhilo,colimlolo);

      for(int i=0; i<dim; i++)
      {
         invUrehihi[i][j] = colrehihi[i]; invUrelohi[i][j] = colrelohi[i];
         invUrehilo[i][j] = colrehilo[i]; invUrelolo[i][j] = colrelolo[i];
         invUimhihi[i][j] = colimhihi[i]; invUimlohi[i][j] = colimlohi[i];
         invUimhilo[i][j] = colimhilo[i]; invUimlolo[i][j] = colimlolo[i];
      }
      rhsrehihi[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhsrehihi); free(rhsrelohi); free(rhsrehilo); free(rhsrelolo);
   free(colrehihi); free(colrelohi); free(colrehilo); free(colrelolo);
   free(rhsimhihi); free(rhsimlohi); free(rhsimhilo); free(rhsimlolo);
   free(colimhihi); free(colimlohi); free(colimhilo); free(colimlolo);
}

void CPU_dbl4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Uhihi, double **Ulohi, double **Uhilo, double **Ulolo,
   double *bhihi, double *blohi, double *bhilo, double *blolo,
   double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *lapsec )
{
   double acchihi,acclohi,acchilo,acclolo;
   double prodhihi,prodlohi,prodhilo,prodlolo;
   double **Thihi = new double*[szt];
   double **Tlohi = new double*[szt];
   double **Thilo = new double*[szt];
   double **Tlolo = new double*[szt];
   double **invThihi = new double*[szt];
   double **invTlohi = new double*[szt];
   double **invThilo = new double*[szt];
   double **invTlolo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Thihi[i] = new double[szt];
      Tlohi[i] = new double[szt];
      Thilo[i] = new double[szt];
      Tlolo[i] = new double[szt];
      invThihi[i] = new double[szt];
      invTlohi[i] = new double[szt];
      invThilo[i] = new double[szt];
      invTlolo[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Thihi[i][j] = Uhihi[idx+i][idx+j];
         Tlohi[i][j] = Ulohi[idx+i][idx+j];
         Thilo[i][j] = Uhilo[idx+i][idx+j];
         Tlolo[i][j] = Ulolo[idx+i][idx+j];
      }

   CPU_dbl4_upper_inverse
      (szt,Thihi,Tlohi,Thilo,Tlolo,
       invThihi,invTlohi,invThilo,invTlolo,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Uhihi[idx+i][idx+j] = invThihi[i][j];
         Ulohi[idx+i][idx+j] = invTlohi[i][j];
         Uhilo[idx+i][idx+j] = invThilo[i][j];
         Ulolo[idx+i][idx+j] = invTlolo[i][j];
      }

   for(int i=0; i<szt; i++)
   {
      xhihi[idx+i] = 0.0; xlohi[idx+i] = 0.0;
      xhilo[idx+i] = 0.0; xlolo[idx+i] = 0.0;

      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         qdf_mul(invThihi[i][j],invTlohi[i][j],invThilo[i][j],invTlolo[i][j],
                    bhihi[idx+j],  blohi[idx+j],  bhilo[idx+j],  blolo[idx+j],
                 &acchihi,      &acclohi,      &acchilo,      &acclolo);
         qdf_inc(&xhihi[idx+i],&xlohi[idx+i],&xhilo[idx+i],&xlolo[idx+i],
                 acchihi,      acclohi,      acchilo,      acclolo);
      }
   }
   double *wbhihi = new double[szt];    // work space for bhi
   double *wblohi = new double[szt];    // work space for bhi
   double *wbhilo = new double[szt];    // work space for blo
   double *wblolo = new double[szt];    // work space for blo
   double **wThihi = new double*[szt];  // work space for a tile
   double **wTlohi = new double*[szt];  // work space for a tile
   double **wThilo = new double*[szt];  // work space for a tile
   double **wTlolo = new double*[szt];  // work space for a tile

   for(int i=0; i<szt; i++)
   {
      wThihi[i] = new double[szt];
      wTlohi[i] = new double[szt];
      wThilo[i] = new double[szt];
      wTlolo[i] = new double[szt];
   }
   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt;  // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Thihi[i][j] = Uhihi[idx+i][idx+j];
            Tlohi[i][j] = Ulohi[idx+i][idx+j];
            Thilo[i][j] = Uhilo[idx+i][idx+j];
            Tlolo[i][j] = Ulolo[idx+i][idx+j];
         }

      CPU_dbl4_upper_inverse
         (szt,Thihi,Tlohi,Thilo,Tlolo,invThihi,invTlohi,invThilo,invTlolo,
          &timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhihi[idx+i][idx+j] = invThihi[i][j];
            Ulohi[idx+i][idx+j] = invTlohi[i][j];
            Uhilo[idx+i][idx+j] = invThilo[i][j];
            Ulolo[idx+i][idx+j] = invTlolo[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++)  // load the work space
         {
            wbhihi[i] = bhihi[rowidx+i];
            wblohi[i] = blohi[rowidx+i];
            wbhilo[i] = bhilo[rowidx+i];
            wblolo[i] = blolo[rowidx+i];

            for(int j=0; j<szt; j++)
            {
               wThihi[i][j] = Uhihi[rowidx+i][idx+szt+j];
               wTlohi[i][j] = Ulohi[rowidx+i][idx+szt+j];
               wThilo[i][j] = Uhilo[rowidx+i][idx+szt+j];
               wTlolo[i][j] = Ulolo[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodhihi = 0.0; prodlohi = 0.0;
            prodhilo = 0.0; prodlolo = 0.0;

            for(int j=0; j<szt; j++)
            {  // prod = prod + wT[i][j]*x[idx+szt+j];
               qdf_mul(wThihi[i][j],wTlohi[i][j],wThilo[i][j],wTlolo[i][j],
                        xhihi[idx+szt+j],xlohi[idx+szt+j],
                        xhilo[idx+szt+j],xlolo[idx+szt+j],
                     &acchihi,&acclohi,&acchilo,&acclolo);
               qdf_inc(&prodhihi,&prodlohi,&prodhilo,&prodlolo,
                         acchihi,  acclohi,  acchilo,  acclolo);
            }
            // wb[i] = wb[i] - prod;
            qdf_dec(&wbhihi[i],&wblohi[i],&wbhilo[i],&wblolo[i],
                   prodhihi,  prodlohi,  prodhilo,  prodlolo);
         }
         for(int i=0; i<szt; i++) // store wb into b for next update
         {
            bhihi[rowidx+i] = wbhihi[i];
            blohi[rowidx+i] = wblohi[i];
            bhilo[rowidx+i] = wbhilo[i];
            blolo[rowidx+i] = wblolo[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodhihi = 0.0; prodlohi = 0.0;
         prodhilo = 0.0; prodlolo = 0.0;

         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            qdf_mul(invThihi[i][j],invTlohi[i][j],
                    invThilo[i][j],invTlolo[i][j],
                    bhihi[idx+j],blohi[idx+j],bhilo[idx+j],blolo[idx+j],
                 &acchihi,    &acclohi,    &acchilo,    &acclolo);
            qdf_inc(&prodhihi,&prodlohi,&prodhilo,&prodlolo,
                      acchihi,  acclohi,  acchilo,  acclolo);
         }
         wbhihi[i] = prodhihi; wblohi[i] = prodlohi;
         wbhilo[i] = prodhilo; wblolo[i] = prodlolo;
      }
      for(int i=0; i<szt; i++)
      {
         xhihi[idx+i] = wbhihi[i]; xlohi[idx+i] = wblohi[i];
         xhilo[idx+i] = wbhilo[i]; xlolo[idx+i] = wblolo[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Thihi[i]); free(Tlohi[i]); free(Thilo[i]); free(Tlolo[i]);
      free(invThihi[i]); free(invTlohi[i]);
      free(invThilo[i]); free(invTlolo[i]);
      free(wThihi[i]); free(wTlohi[i]); free(wThilo[i]); free(wTlolo[i]);
   }
   free(Thihi); free(Tlohi); free(Thilo); free(Tlolo);
   free(invThihi); free(invTlohi); free(invThilo); free(invTlolo);
   free(wbhihi); free(wblohi); free(wbhilo); free(wblolo);
   free(wThihi); free(wTlohi); free(wThilo); free(wTlolo);
}

void CPU_cmplx4_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehihi, double **Urelohi, double **Urehilo, double **Urelolo,
   double **Uimhihi, double **Uimlohi, double **Uimhilo, double **Uimlolo,
   double *brehihi, double *brelohi, double *brehilo, double *brelolo,
   double *bimhihi, double *bimlohi, double *bimhilo, double *bimlolo,
   double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *lapsec )
{
   double acc1hihi,acc1lohi,acc1hilo,acc1lolo;
   double acc2hihi,acc2lohi,acc2hilo,acc2lolo;
   double acc3hihi,acc3lohi,acc3hilo,acc3lolo;
   double acc4hihi,acc4lohi,acc4hilo,acc4lolo;
   double prodrehihi,prodrelohi,prodrehilo,prodrelolo;
   double prodimhihi,prodimlohi,prodimhilo,prodimlolo;
   double **Trehihi = new double*[szt];
   double **Trelohi = new double*[szt];
   double **Trehilo = new double*[szt];
   double **Trelolo = new double*[szt];
   double **Timhihi = new double*[szt];
   double **Timlohi = new double*[szt];
   double **Timhilo = new double*[szt];
   double **Timlolo = new double*[szt];
   double **invTrehihi = new double*[szt];
   double **invTrelohi = new double*[szt];
   double **invTrehilo = new double*[szt];
   double **invTrelolo = new double*[szt];
   double **invTimhihi = new double*[szt];
   double **invTimlohi = new double*[szt];
   double **invTimhilo = new double*[szt];
   double **invTimlolo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Trehihi[i] = new double[szt];
      Trelohi[i] = new double[szt];
      Trehilo[i] = new double[szt];
      Trelolo[i] = new double[szt];
      Timhihi[i] = new double[szt];
      Timlohi[i] = new double[szt];
      Timhilo[i] = new double[szt];
      Timlolo[i] = new double[szt];
      invTrehihi[i] = new double[szt];
      invTrelohi[i] = new double[szt];
      invTrehilo[i] = new double[szt];
      invTrelolo[i] = new double[szt];
      invTimhihi[i] = new double[szt];
      invTimlohi[i] = new double[szt];
      invTimhilo[i] = new double[szt];
      invTimlolo[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Trehihi[i][j] = Urehihi[idx+i][idx+j];
         Trelohi[i][j] = Urelohi[idx+i][idx+j];
         Trehilo[i][j] = Urehilo[idx+i][idx+j];
         Trelolo[i][j] = Urelolo[idx+i][idx+j];
         Timhihi[i][j] = Uimhihi[idx+i][idx+j];
         Timlohi[i][j] = Uimlohi[idx+i][idx+j];
         Timhilo[i][j] = Uimhilo[idx+i][idx+j];
         Timlolo[i][j] = Uimlolo[idx+i][idx+j];
      }

   CPU_cmplx4_upper_inverse
      (szt,Trehihi,Trelohi,Trehilo,Trelolo,Timhihi,Timlohi,Timhilo,Timlolo,
           invTrehihi,invTrelohi,invTrehilo,invTrelolo,
           invTimhihi,invTimlohi,invTimhilo,invTimlolo,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Urehihi[idx+i][idx+j] = invTrehihi[i][j];
         Urelohi[idx+i][idx+j] = invTrelohi[i][j];
         Urehilo[idx+i][idx+j] = invTrehilo[i][j];
         Urelolo[idx+i][idx+j] = invTrelolo[i][j];
         Uimhihi[idx+i][idx+j] = invTimhihi[i][j];
         Uimlohi[idx+i][idx+j] = invTimlohi[i][j];
         Uimhilo[idx+i][idx+j] = invTimhilo[i][j];
         Uimlolo[idx+i][idx+j] = invTimlolo[i][j];
      }

   for(int i=0; i<szt; i++)
   {
      xrehihi[idx+i] = 0.0; xrelohi[idx+i] = 0.0;
      xrehilo[idx+i] = 0.0; xrelolo[idx+i] = 0.0;
      ximhihi[idx+i] = 0.0; ximlohi[idx+i] = 0.0;
      ximhilo[idx+i] = 0.0; ximlolo[idx+i] = 0.0;

      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         qdf_mul(invTrehihi[i][j],invTrelohi[i][j],
                 invTrehilo[i][j],invTrelolo[i][j],
                 brehihi[idx+j],brelohi[idx+j],brehilo[idx+j],brelolo[idx+j],
                 &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
         qdf_mul(invTimhihi[i][j],invTimlohi[i][j],
                 invTimhilo[i][j],invTimlolo[i][j],
                 bimhihi[idx+j],bimlohi[idx+j],bimhilo[idx+j],bimlolo[idx+j],
                 &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
         qdf_mul(invTimhihi[i][j],invTimlohi[i][j],
                 invTimhilo[i][j],invTimlolo[i][j],
                 brehihi[idx+j],brelohi[idx+j],brehilo[idx+j],brelolo[idx+j],
                 &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
         qdf_mul(invTrehihi[i][j],invTrelohi[i][j],
                 invTrehilo[i][j],invTrelolo[i][j],
                 bimhihi[idx+j],bimlohi[idx+j],bimhilo[idx+j],bimlolo[idx+j],
                 &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
         qdf_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                  acc2hihi, acc2lohi, acc2hilo, acc2lolo);
         qdf_inc(&xrehihi[idx+i],&xrelohi[idx+i],
                 &xrehilo[idx+i],&xrelolo[idx+i],
                 acc1hihi,acc1lohi,acc1hilo,acc1lolo);
         qdf_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
                  acc4hihi, acc4lohi, acc4hilo, acc4lolo);
         qdf_inc(&ximhihi[idx+i],&ximlohi[idx+i],
                 &ximhilo[idx+i],&ximlolo[idx+i],
                 acc3hihi,acc3lohi,acc3hilo,acc3lolo);
      }
   }
   double *wbrehihi = new double[szt];    // work space for brehi
   double *wbrelohi = new double[szt];    // work space for brehi
   double *wbrehilo = new double[szt];    // work space for brelo
   double *wbrelolo = new double[szt];    // work space for brelo
   double *wbimhihi = new double[szt];    // work space for bimhi
   double *wbimlohi = new double[szt];    // work space for bimhi
   double *wbimhilo = new double[szt];    // work space for bimlo
   double *wbimlolo = new double[szt];    // work space for bimlo
   double **wTrehihi = new double*[szt];  // work space for a tile
   double **wTrelohi = new double*[szt];  // work space for a tile
   double **wTrehilo = new double*[szt];  // work space for a tile
   double **wTrelolo = new double*[szt];  // work space for a tile
   double **wTimhihi = new double*[szt];  // work space for a tile
   double **wTimlohi = new double*[szt];  // work space for a tile
   double **wTimhilo = new double*[szt];  // work space for a tile
   double **wTimlolo = new double*[szt];  // work space for a tile
 
   for(int i=0; i<szt; i++)
   {
      wTrehihi[i] = new double[szt];
      wTrelohi[i] = new double[szt];
      wTrehilo[i] = new double[szt];
      wTrelolo[i] = new double[szt];
      wTimhihi[i] = new double[szt];
      wTimlohi[i] = new double[szt];
      wTimhilo[i] = new double[szt];
      wTimlolo[i] = new double[szt];
   }
   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt;  // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Trehihi[i][j] = Urehihi[idx+i][idx+j];
            Trelohi[i][j] = Urelohi[idx+i][idx+j];
            Trehilo[i][j] = Urehilo[idx+i][idx+j];
            Trelolo[i][j] = Urelolo[idx+i][idx+j];
            Timhihi[i][j] = Uimhihi[idx+i][idx+j];
            Timlohi[i][j] = Uimlohi[idx+i][idx+j];
            Timhilo[i][j] = Uimhilo[idx+i][idx+j];
            Timlolo[i][j] = Uimlolo[idx+i][idx+j];
         }

      CPU_cmplx4_upper_inverse
         (szt,Trehihi,Trelohi,Trehilo,Trelolo,
              Timhihi,Timlohi,Timhilo,Timlolo,
          invTrehihi,invTrelohi,invTrehilo,invTrelolo,
          invTimhihi,invTimlohi,invTimhilo,invTimlolo,&timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Urehihi[idx+i][idx+j] = invTrehihi[i][j];
            Urelohi[idx+i][idx+j] = invTrelohi[i][j];
            Urehilo[idx+i][idx+j] = invTrehilo[i][j];
            Urelolo[idx+i][idx+j] = invTrelolo[i][j];
            Uimhihi[idx+i][idx+j] = invTimhihi[i][j];
            Uimlohi[idx+i][idx+j] = invTimlohi[i][j];
            Uimhilo[idx+i][idx+j] = invTimhilo[i][j];
            Uimlolo[idx+i][idx+j] = invTimlolo[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++)  // load the work space
         {
            wbrehihi[i] = brehihi[rowidx+i]; wbrelohi[i] = brelohi[rowidx+i];
            wbrehilo[i] = brehilo[rowidx+i]; wbrelolo[i] = brelolo[rowidx+i];
            wbimhihi[i] = bimhihi[rowidx+i]; wbimlohi[i] = bimlohi[rowidx+i];
            wbimhilo[i] = bimhilo[rowidx+i]; wbimlolo[i] = bimlolo[rowidx+i];

            for(int j=0; j<szt; j++)
            {
               wTrehihi[i][j] = Urehihi[rowidx+i][idx+szt+j];
               wTrehilo[i][j] = Urehilo[rowidx+i][idx+szt+j];
               wTrelohi[i][j] = Urelohi[rowidx+i][idx+szt+j];
               wTrelolo[i][j] = Urelolo[rowidx+i][idx+szt+j];
               wTimhihi[i][j] = Uimhihi[rowidx+i][idx+szt+j];
               wTimhilo[i][j] = Uimhilo[rowidx+i][idx+szt+j];
               wTimlohi[i][j] = Uimlohi[rowidx+i][idx+szt+j];
               wTimlolo[i][j] = Uimlolo[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodrehihi = 0.0; prodrelohi = 0.0;
            prodrehilo = 0.0; prodrelolo = 0.0;
            prodimhihi = 0.0; prodimlohi = 0.0;
            prodimhilo = 0.0; prodimlolo = 0.0;

            for(int j=0; j<szt; j++) // prod = prod + wT[i][j]*x[idx+szt+j];
            {
               qdf_mul(wTrehihi[i][j],wTrelohi[i][j],
                       wTrehilo[i][j],wTrelolo[i][j],
                       xrehihi[idx+szt+j],xrelohi[idx+szt+j],
                       xrehilo[idx+szt+j],xrelolo[idx+szt+j],
                       &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
               qdf_mul(wTimhihi[i][j],wTimlohi[i][j],
                       wTimhilo[i][j],wTimlolo[i][j],
                       ximhihi[idx+szt+j],ximlohi[idx+szt+j],
                       ximhilo[idx+szt+j],ximlolo[idx+szt+j],
                       &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
               qdf_mul(wTimhihi[i][j],wTimlohi[i][j],
                       wTimhilo[i][j],wTimlolo[i][j],
                       xrehihi[idx+szt+j],xrelohi[idx+szt+j],
                       xrehilo[idx+szt+j],xrelolo[idx+szt+j],
                       &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
               qdf_mul(wTrehihi[i][j],wTrelohi[i][j],
                       wTrehilo[i][j],wTrelolo[i][j],
                       ximhihi[idx+szt+j],ximlohi[idx+szt+j],
                       ximhilo[idx+szt+j],ximlolo[idx+szt+j],
                       &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
               qdf_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                        acc2hihi, acc2lohi, acc2hilo, acc2lolo);
               qdf_inc(&prodrehihi,&prodrelohi,&prodrehilo,&prodrelolo,
                          acc1hihi,   acc1lohi,   acc1hilo,   acc1lolo);
               qdf_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
                        acc4hihi, acc4lohi, acc4hilo, acc4lolo);
               qdf_inc(&prodimhihi,&prodimlohi,&prodimhilo,&prodimlolo,
                          acc3hihi,   acc3lohi,   acc3hilo,   acc3lolo);
            }
            // wb[i] = wb[i] - prod;
            qdf_dec(&wbrehihi[i],&wbrelohi[i],&wbrehilo[i],&wbrelolo[i],
                   prodrehihi,  prodrelohi,  prodrehilo,  prodrelolo);
            qdf_dec(&wbimhihi[i],&wbimlohi[i],&wbimhilo[i],&wbimlolo[i],
                   prodimhihi,  prodimlohi,  prodimhilo,  prodimlolo);
         }
         for(int i=0; i<szt; i++) // store wb into b for next update
         {
            brehihi[rowidx+i] = wbrehihi[i]; brelohi[rowidx+i] = wbrelohi[i];
            brehilo[rowidx+i] = wbrehilo[i]; brelolo[rowidx+i] = wbrelolo[i];
            bimhihi[rowidx+i] = wbimhihi[i]; bimlohi[rowidx+i] = wbimlohi[i];
            bimhilo[rowidx+i] = wbimhilo[i]; bimlolo[rowidx+i] = wbimlolo[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodrehihi = 0.0; prodrelohi = 0.0;
         prodrehilo = 0.0; prodrelolo = 0.0;
         prodimhihi = 0.0; prodimlohi = 0.0;
         prodimhilo = 0.0; prodimlolo = 0.0;

         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            qdf_mul(invTrehihi[i][j],invTrelohi[i][j],
                    invTrehilo[i][j],invTrelolo[i][j],
                    brehihi[idx+j],brelohi[idx+j],
                    brehilo[idx+j],brelolo[idx+j],
                    &acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo);
            qdf_mul(invTimhihi[i][j],invTimlohi[i][j],
                    invTimhilo[i][j],invTimlolo[i][j],
                    bimhihi[idx+j],bimlohi[idx+j],
                    bimhilo[idx+j],bimlolo[idx+j],
                    &acc2hihi,&acc2lohi,&acc2hilo,&acc2lolo);
            qdf_mul(invTimhihi[i][j],invTimlohi[i][j],
                    invTimhilo[i][j],invTimlolo[i][j],
                    brehihi[idx+j],brelohi[idx+j],
                    brehilo[idx+j],brelolo[idx+j],
                    &acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo);
            qdf_mul(invTrehihi[i][j],invTrelohi[i][j],
                    invTrehilo[i][j],invTrelolo[i][j],
                    bimhihi[idx+j],bimlohi[idx+j],
                    bimhilo[idx+j],bimlolo[idx+j],
                    &acc4hihi,&acc4lohi,&acc4hilo,&acc4lolo);
            qdf_dec(&acc1hihi,&acc1lohi,&acc1hilo,&acc1lolo,
                     acc2hihi, acc2lohi, acc2hilo, acc2lolo);
            qdf_inc(&prodrehihi,&prodrelohi,&prodrehilo,&prodrelolo,
                       acc1hihi,   acc1lohi,   acc1hilo,   acc1lolo);
            qdf_inc(&acc3hihi,&acc3lohi,&acc3hilo,&acc3lolo,
                     acc4hihi, acc4lohi, acc4hilo, acc4lolo);
            qdf_inc(&prodimhihi,&prodimlohi,&prodimhilo,&prodimlolo,
                       acc3hihi,   acc3lohi,   acc3hilo,   acc3lolo);
         }
         wbrehihi[i] = prodrehihi; wbrelohi[i] = prodrelohi;
         wbrehilo[i] = prodrehilo; wbrelolo[i] = prodrelolo;
         wbimhihi[i] = prodimhihi; wbimlohi[i] = prodimlohi;
         wbimhilo[i] = prodimhilo; wbimlolo[i] = prodimlolo;
      }
      for(int i=0; i<szt; i++)
      {
         xrehihi[idx+i] = wbrehihi[i]; xrelohi[idx+i] = wbrelohi[i];
         xrehilo[idx+i] = wbrehilo[i]; xrelolo[idx+i] = wbrelolo[i];
         ximhihi[idx+i] = wbimhihi[i]; ximlohi[idx+i] = wbimlohi[i];
         ximhilo[idx+i] = wbimhilo[i]; ximlolo[idx+i] = wbimlolo[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Trehihi[i]); free(Trelohi[i]); free(Trehilo[i]); free(Trelolo[i]);
      free(Timhihi[i]); free(Timlohi[i]); free(Timhilo[i]); free(Timlolo[i]);
      free(invTrehihi[i]); free(invTrelohi[i]);
      free(invTrehilo[i]); free(invTrelolo[i]);
      free(invTimhihi[i]); free(invTimlohi[i]);
      free(invTimhilo[i]); free(invTimlolo[i]);
      free(wTrehihi[i]); free(wTrelohi[i]);
      free(wTrehilo[i]); free(wTrelolo[i]);
      free(wTimhihi[i]); free(wTimlohi[i]);
      free(wTimhilo[i]); free(wTimlolo[i]);
   }
   free(Trehihi); free(Trelohi); free(Trehilo); free(Trelolo);
   free(Timhihi); free(Timlohi); free(Timhilo); free(Timlolo);
   free(invTrehihi); free(invTrelohi); 
   free(invTrehilo); free(invTrelolo); 
   free(invTimhihi); free(invTimlohi);
   free(invTimhilo); free(invTimlolo);
   free(wbrehihi); free(wbrelohi); free(wbrehilo); free(wbrelolo);
   free(wbimhihi); free(wbimlohi); free(wbimhilo); free(wbimlolo);
   free(wTrehihi); free(wTrelohi); free(wTrehilo); free(wTrelolo);
   free(wTimhihi); free(wTimlohi); free(wTimhilo); free(wTimlolo);
}
