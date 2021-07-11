/* The file dbl2_tabs_host.cpp defines functions specified in
 * the file dbl2_tabs_host.h. */

#include <cstdlib>
#include <ctime>
#include "double_double_functions.h"
#include "dbl2_factorizations.h"
#include "dbl2_tabs_host.h"

using namespace std;

void CPU_dbl2_backsubs
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo )
{
   CPU_dbl2_factors_backward(dim,Uhi,Ulo,bhi,blo,xhi,xlo);
}

void CPU_cmplx2_backsubs
 ( int dim, double **Urehi, double **Urelo,
            double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo )
{
   CPU_cmplx2_factors_backward
      (dim,Urehi,Urelo,Uimhi,Uimlo,
           brehi,brelo,bimhi,bimlo,
           xrehi,xrelo,ximhi,ximlo);
}

void CPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo,
   double *lapsec )
{
   double *colhi = new double[dim];
   double *collo = new double[dim];
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhshi[i] = 0.0;
      rhslo[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhshi[j] = 1.0;
      CPU_dbl2_backsubs(dim,Uhi,Ulo,rhshi,rhslo,colhi,collo);
      for(int i=0; i<dim; i++)
      {
         invUhi[i][j] = colhi[i];
         invUlo[i][j] = collo[i];
      }
      rhshi[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhshi); free(colhi);
   free(rhslo); free(collo);
}

void CPU_cmplx2_upper_inverse
 ( int dim, double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double **invUrehi, double **invUrelo,
   double **invUimhi, double **invUimlo, double *lapsec )
{
   double *colrehi = new double[dim];
   double *colrelo = new double[dim];
   double *colimhi = new double[dim];
   double *colimlo = new double[dim];
   double *rhsrehi = new double[dim];
   double *rhsrelo = new double[dim];
   double *rhsimhi = new double[dim];
   double *rhsimlo = new double[dim];

   clock_t start = clock();

   for(int i=0; i<dim; i++)
   {
      rhsrehi[i] = 0.0; rhsrelo[i] = 0.0;
      rhsimhi[i] = 0.0; rhsimlo[i] = 0.0;
   }
   for(int j=0; j<dim; j++) // compute j-th column of the inverse
   {
      rhsrehi[j] = 1.0;
      CPU_cmplx2_backsubs(dim,  Urehi,  Urelo,  Uimhi,  Uimlo,
                              rhsrehi,rhsrelo,rhsimhi,rhsimlo,
                              colrehi,colrelo,colimhi,colimlo);
      for(int i=0; i<dim; i++)
      {
         invUrehi[i][j] = colrehi[i]; invUrelo[i][j] = colrelo[i];
         invUimhi[i][j] = colimhi[i]; invUimlo[i][j] = colimlo[i];
      }
      rhsrehi[j] = 0.0;
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   free(rhsrehi); free(colrehi);
   free(rhsrelo); free(colrelo);
   free(rhsimhi); free(colimhi);
   free(rhsimlo); free(colimlo);
}

void CPU_dbl2_matmatmul
 ( int dim, double **Ahi, double **Alo, double **Fhi, double **Flo )
{
   double **resulthi = new double*[dim];
   double **resultlo = new double*[dim];
   double acchi,acclo;

   for(int i=0; i<dim; i++)
   {
      resulthi[i] = new double[dim];
      resultlo[i] = new double[dim];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         resulthi[i][j] = 0.0;
         resultlo[i][j] = 0.0;
         for(int k=0; k<dim; k++)
         {  // result[i][j] = result[i][j] + F[i][k]*A[k][j];
            ddf_mul(Fhi[i][k],Flo[i][k],Ahi[k][j],Alo[k][j],&acchi,&acclo);
            ddf_inc(&resulthi[i][j],&resultlo[i][j],acchi,acclo);
         }
      }

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
      {
         Ahi[i][j] = resulthi[i][j];
         Alo[i][j] = resultlo[i][j];
      }

   for(int i=0; i<dim; i++)
   {
      free(resulthi[i]);
      free(resultlo[i]);
   }
   free(resulthi); free(resultlo);
}

void CPU_dbl2_upper_tiled_solver
 ( int dim, int szt, int nbt, double **Uhi, double **Ulo,
   double *bhi, double *blo, double *xhi, double *xlo, double *lapsec )
{
   double acchi,acclo,prodhi,prodlo;
   double **Thi = new double*[szt];
   double **Tlo = new double*[szt];
   double **invThi = new double*[szt];
   double **invTlo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Thi[i] = new double[szt];
      Tlo[i] = new double[szt];
      invThi[i] = new double[szt];
      invTlo[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Thi[i][j] = Uhi[idx+i][idx+j];
         Tlo[i][j] = Ulo[idx+i][idx+j];
      }

   CPU_dbl2_upper_inverse(szt,Thi,Tlo,invThi,invTlo,&timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Uhi[idx+i][idx+j] = invThi[i][j];
         Ulo[idx+i][idx+j] = invTlo[i][j];
      }

   for(int i=0; i<szt; i++)
   {
      xhi[idx+i] = 0.0;
      xlo[idx+i] = 0.0;
      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         ddf_mul(invThi[i][j],invTlo[i][j],bhi[idx+j],blo[idx+j],
                 &acchi,&acclo);
         ddf_inc(&xhi[idx+i],&xlo[idx+i],acchi,acclo);
      }
   }
   double *wbhi = new double[szt];    // work space for bhi
   double *wblo = new double[szt];    // work space for blo
   double **wThi = new double*[szt];  // work space for a tile
   double **wTlo = new double*[szt];  // work space for a tile
   for(int i=0; i<szt; i++)
   {
      wThi[i] = new double[szt];
      wTlo[i] = new double[szt];
   }
   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt;  // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Thi[i][j] = Uhi[idx+i][idx+j];
            Tlo[i][j] = Ulo[idx+i][idx+j];
         }

      CPU_dbl2_upper_inverse(szt,Thi,Tlo,invThi,invTlo,&timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Uhi[idx+i][idx+j] = invThi[i][j];
            Ulo[idx+i][idx+j] = invTlo[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++)  // load the work space
         {
            wbhi[i] = bhi[rowidx+i];
            wblo[i] = blo[rowidx+i];
            for(int j=0; j<szt; j++)
            {
               wThi[i][j] = Uhi[rowidx+i][idx+szt+j];
               wTlo[i][j] = Ulo[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodhi = 0.0;
            prodlo = 0.0;
            for(int j=0; j<szt; j++)
            {  // prod = prod + wT[i][j]*x[idx+szt+j];
               ddf_mul(wThi[i][j],wTlo[i][j],xhi[idx+szt+j],
                       xlo[idx+szt+j],&acchi,&acclo);
               ddf_inc(&prodhi,&prodlo,acchi,acclo);
            }
            // wb[i] = wb[i] - prod;
            ddf_dec(&wbhi[i],&wblo[i],prodhi,prodlo);
         }
         for(int i=0; i<szt; i++) // store wb into b for next update
         {
            bhi[rowidx+i] = wbhi[i];
            blo[rowidx+i] = wblo[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodhi = 0.0;
         prodlo = 0.0;
         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            ddf_mul(invThi[i][j],invTlo[i][j],bhi[idx+j],blo[idx+j],
                    &acchi,&acclo);
            ddf_inc(&prodhi,&prodlo,acchi,acclo);
         }
         wbhi[i] = prodhi;
         wblo[i] = prodlo;
      }
      for(int i=0; i<szt; i++)
      {
         xhi[idx+i] = wbhi[i];
         xlo[idx+i] = wblo[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Thi[i]); free(invThi[i]); free(wThi[i]);
      free(Tlo[i]); free(invTlo[i]); free(wTlo[i]);
   }
   free(Thi); free(invThi); free(wbhi); free(wThi);
   free(Tlo); free(invTlo); free(wblo); free(wTlo);
}

void CPU_cmplx2_upper_tiled_solver
 ( int dim, int szt, int nbt,
   double **Urehi, double **Urelo, double **Uimhi, double **Uimlo,
   double *brehi, double *brelo, double *bimhi, double *bimlo,
   double *xrehi, double *xrelo, double *ximhi, double *ximlo,
   double *lapsec )
{
   double acc1hi,acc1lo,acc2hi,acc2lo;
   double acc3hi,acc3lo,acc4hi,acc4lo;
   double prodrehi,prodrelo,prodimhi,prodimlo;
   double **Trehi = new double*[szt];
   double **Trelo = new double*[szt];
   double **Timhi = new double*[szt];
   double **Timlo = new double*[szt];
   double **invTrehi = new double*[szt];
   double **invTrelo = new double*[szt];
   double **invTimhi = new double*[szt];
   double **invTimlo = new double*[szt];

   for(int i=0; i<szt; i++)
   {
      Trehi[i] = new double[szt];
      Trelo[i] = new double[szt];
      Timhi[i] = new double[szt];
      Timlo[i] = new double[szt];
      invTrehi[i] = new double[szt];
      invTrelo[i] = new double[szt];
      invTimhi[i] = new double[szt];
      invTimlo[i] = new double[szt];
   }
   double timelapsed;

   clock_t start = clock();

   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Trehi[i][j] = Urehi[idx+i][idx+j];
         Trelo[i][j] = Urelo[idx+i][idx+j];
         Timhi[i][j] = Uimhi[idx+i][idx+j];
         Timlo[i][j] = Uimlo[idx+i][idx+j];
      }

   CPU_cmplx2_upper_inverse
      (szt,Trehi,Trelo,Timhi,Timlo,invTrehi,invTrelo,invTimhi,invTimlo,
       &timelapsed);

   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Urehi[idx+i][idx+j] = invTrehi[i][j];
         Urelo[idx+i][idx+j] = invTrelo[i][j];
         Uimhi[idx+i][idx+j] = invTimhi[i][j];
         Uimlo[idx+i][idx+j] = invTimlo[i][j];
      }

   for(int i=0; i<szt; i++)
   {
      xrehi[idx+i] = 0.0; xrelo[idx+i] = 0.0;
      ximhi[idx+i] = 0.0; ximlo[idx+i] = 0.0;

      for(int j=0; j<szt; j++) // x[idx+i] = x[idx+i] + invT[i][j]*b[idx+j];
      {
         ddf_mul(invTrehi[i][j],invTrelo[i][j],
                    brehi[idx+j],  brelo[idx+j],&acc1hi,&acc1lo);
         ddf_mul(invTimhi[i][j],invTimlo[i][j],
                    bimhi[idx+j],  bimlo[idx+j],&acc2hi,&acc2lo);
         ddf_mul(invTimhi[i][j],invTimlo[i][j],
                    brehi[idx+j],  brelo[idx+j],&acc3hi,&acc3lo);
         ddf_mul(invTrehi[i][j],invTrelo[i][j],
                    bimhi[idx+j],  bimlo[idx+j],&acc4hi,&acc4lo);
         ddf_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
         ddf_inc(&xrehi[idx+i],&xrelo[idx+i],acc1hi,acc1lo);
         ddf_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
         ddf_inc(&ximhi[idx+i],&ximlo[idx+i],acc3hi,acc3lo);
      }
   }
   double *wbrehi = new double[szt];    // work space for brehi
   double *wbrelo = new double[szt];    // work space for brelo
   double *wbimhi = new double[szt];    // work space for bimhi
   double *wbimlo = new double[szt];    // work space for bimlo
   double **wTrehi = new double*[szt];  // work space for a tile
   double **wTrelo = new double*[szt];  // work space for a tile
   double **wTimhi = new double*[szt];  // work space for a tile
   double **wTimlo = new double*[szt];  // work space for a tile

   for(int i=0; i<szt; i++)
   {
      wTrehi[i] = new double[szt];
      wTrelo[i] = new double[szt];
      wTimhi[i] = new double[szt];
      wTimlo[i] = new double[szt];
   }
   for(int k=nbt-1; k>0; k--)  // update with solution tile k
   {
      idx = idx - szt;  // idx is start index of diagonal tile

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Trehi[i][j] = Urehi[idx+i][idx+j];
            Trelo[i][j] = Urelo[idx+i][idx+j];
            Timhi[i][j] = Uimhi[idx+i][idx+j];
            Timlo[i][j] = Uimlo[idx+i][idx+j];
         }

      CPU_cmplx2_upper_inverse
         (szt,Trehi,Trelo,Timhi,Timlo,invTrehi,invTrelo,invTimhi,invTimlo,
          &timelapsed);

      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
         {
            Urehi[idx+i][idx+j] = invTrehi[i][j];
            Urelo[idx+i][idx+j] = invTrelo[i][j];
            Uimhi[idx+i][idx+j] = invTimhi[i][j];
            Uimlo[idx+i][idx+j] = invTimlo[i][j];
         }

      for(int L=0; L<k; L++)   // update wb as many times as k
      {
         int rowidx = L*szt;

         for(int i=0; i<szt; i++)  // load the work space
         {
            wbrehi[i] = brehi[rowidx+i]; wbrelo[i] = brelo[rowidx+i];
            wbimhi[i] = bimhi[rowidx+i]; wbimlo[i] = bimlo[rowidx+i];

            for(int j=0; j<szt; j++)
            {
               wTrehi[i][j] = Urehi[rowidx+i][idx+szt+j];
               wTrelo[i][j] = Urelo[rowidx+i][idx+szt+j];
               wTimhi[i][j] = Uimhi[rowidx+i][idx+szt+j];
               wTimlo[i][j] = Uimlo[rowidx+i][idx+szt+j];
            }
         }
         for(int i=0; i<szt; i++) // update wb
         {
            prodrehi = 0.0; prodrelo = 0.0;
            prodimhi = 0.0; prodimlo = 0.0;

            for(int j=0; j<szt; j++) // prod = prod + wT[i][j]*x[idx+szt+j];
            {
               ddf_mul(wTrehi[i][j],wTrelo[i][j],
                       xrehi[idx+szt+j],xrelo[idx+szt+j],&acc1hi,&acc1lo);
               ddf_mul(wTimhi[i][j],wTimlo[i][j],
                       ximhi[idx+szt+j],ximlo[idx+szt+j],&acc2hi,&acc2lo);
               ddf_mul(wTimhi[i][j],wTimlo[i][j],
                       xrehi[idx+szt+j],xrelo[idx+szt+j],&acc3hi,&acc3lo);
               ddf_mul(wTrehi[i][j],wTrelo[i][j],
                       ximhi[idx+szt+j],ximlo[idx+szt+j],&acc4hi,&acc4lo);
               ddf_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
               ddf_inc(&prodrehi,&prodrelo,acc1hi,acc1lo);
               ddf_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
               ddf_inc(&prodimhi,&prodimlo,acc3hi,acc3lo);
            }
            // wb[i] = wb[i] - prod;
            ddf_dec(&wbrehi[i],&wbrelo[i],prodrehi,prodrelo);
            ddf_dec(&wbimhi[i],&wbimlo[i],prodimhi,prodimlo);
         }
         for(int i=0; i<szt; i++) // store wb into b for next update
         {
            brehi[rowidx+i] = wbrehi[i]; brelo[rowidx+i] = wbrelo[i];
            bimhi[rowidx+i] = wbimhi[i]; bimlo[rowidx+i] = wbimlo[i];
         }
      }
      for(int i=0; i<szt; i++)   // wb = invT*b
      {
         prodrehi = 0.0; prodrelo = 0.0;
         prodimhi = 0.0; prodimlo = 0.0;

         for(int j=0; j<szt; j++) // prod = prod + invT[i][j]*b[idx+j];
         {
            ddf_mul(invTrehi[i][j],invTrelo[i][j],
                    brehi[idx+j],brelo[idx+j],&acc1hi,&acc1lo);
            ddf_mul(invTimhi[i][j],invTimlo[i][j],
                    bimhi[idx+j],bimlo[idx+j],&acc2hi,&acc2lo);
            ddf_mul(invTimhi[i][j],invTimlo[i][j],
                    brehi[idx+j],brelo[idx+j],&acc3hi,&acc3lo);
            ddf_mul(invTrehi[i][j],invTrelo[i][j],
                    bimhi[idx+j],bimlo[idx+j],&acc4hi,&acc4lo);
            ddf_dec(&acc1hi,&acc1lo,acc2hi,acc2lo);
            ddf_inc(&prodrehi,&prodrelo,acc1hi,acc1lo);
            ddf_inc(&acc3hi,&acc3lo,acc4hi,acc4lo);
            ddf_inc(&prodimhi,&prodimlo,acc3hi,acc3lo);
         }
         wbrehi[i] = prodrehi; wbrelo[i] = prodrelo;
         wbimhi[i] = prodimhi; wbimlo[i] = prodimlo;
      }
      for(int i=0; i<szt; i++)
      {
         xrehi[idx+i] = wbrehi[i]; xrelo[idx+i] = wbrelo[i];
         ximhi[idx+i] = wbimhi[i]; ximlo[idx+i] = wbimlo[i];
      }
   }
   clock_t end = clock();
   *lapsec = double(end - start)/CLOCKS_PER_SEC;

   for(int i=0; i<szt; i++)
   {
      free(Trehi[i]); free(invTrehi[i]); free(wTrehi[i]);
      free(Trelo[i]); free(invTrelo[i]); free(wTrelo[i]);
      free(Timhi[i]); free(invTimhi[i]); free(wTimhi[i]);
      free(Timlo[i]); free(invTimlo[i]); free(wTimlo[i]);
   }
   free(Trehi); free(invTrehi); free(wbrehi); free(wTrehi);
   free(Trelo); free(invTrelo); free(wbrelo); free(wTrelo);
   free(Timhi); free(invTimhi); free(wbimhi); free(wTimhi);
   free(Timlo); free(invTimlo); free(wbimlo); free(wTimlo);
}
