/* The file dbl2_tabs_host.cpp defines functions specified in
 * the file dbl2_tabs_host.h. */

#include <cstdlib>
#include "double_double_functions.h"
#include "dbl2_tabs_host.h"

using namespace std;

void CPU_dbl2_backsubs
 ( int dim, double **Uhi, double **Ulo, double *bhi, double *blo,
   double *xhi, double *xlo )
{
   double acchi,acclo;
   int idx = dim-1;
   
   // x[idx] = b[idx]/U[idx][idx];
   ddf_div(bhi[idx],blo[idx],Uhi[idx][idx],Ulo[idx][idx],&xhi[idx],&xlo[idx]);

   for(int i=idx-1; i>=0; i--)
   {
      xhi[i] = bhi[i]; xlo[i] = blo[i];

      for(int j=i+1; j<dim; j++)     // x[i] = x[i] - U[i][j]*x[j];
      {
         ddf_mul(Uhi[i][j],Ulo[i][j],xhi[j],xlo[j],&acchi,&acclo);
         ddf_dec(&xhi[i],&xlo[i],acchi,acclo);
      }
      // x[i] = x[i]/U[i][i];
      ddf_div(xhi[i],xlo[i],Uhi[i][i],Ulo[i][i],&acchi,&acclo);
      xhi[i] = acchi; xlo[i] = acclo;
   }
}

void CPU_dbl2_upper_inverse
 ( int dim, double **Uhi, double **Ulo, double **invUhi, double **invUlo )
{
   double *colhi = new double[dim];
   double *collo = new double[dim];
   double *rhshi = new double[dim];
   double *rhslo = new double[dim];

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
   free(rhshi); free(colhi);
   free(rhslo); free(collo);
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
   double *bhi, double *blo, double *xhi, double *xlo )
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
   int idx = (nbt-1)*szt;
   for(int i=0; i<szt; i++)
      for(int j=0; j<szt; j++)
      {
         Thi[i][j] = Uhi[idx+i][idx+j];
         Tlo[i][j] = Ulo[idx+i][idx+j];
      }

   CPU_dbl2_upper_inverse(szt,Thi,Tlo,invThi,invTlo);

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

      CPU_dbl2_upper_inverse(szt,Thi,Tlo,invThi,invTlo);

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
   for(int i=0; i<szt; i++)
   {
      free(Thi[i]); free(invThi[i]); free(wThi[i]);
      free(Tlo[i]); free(invTlo[i]); free(wTlo[i]);
   }
   free(Thi); free(invThi); free(wbhi); free(wThi);
   free(Tlo); free(invTlo); free(wblo); free(wTlo);
}
