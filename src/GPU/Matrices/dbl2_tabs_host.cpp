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
