// The file dbl4_test_utilities.cpp defines the functions specified in
// the file dbl4_test_utilities.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "quad_double_functions.h"
#include "random4_matrices.h"
#include "dbl4_factorizations.h"
#include "dbl4_test_utilities.h"

using namespace std;

double dbl4_Difference_Sum
 ( int n, double *xhihi, double *xlohi, double *xhilo, double *xlolo,
   double *yhihi, double *ylohi, double *yhilo, double *ylolo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xhihi[i] - yhihi[i]) + abs(xlohi[i] - ylohi[i])
                      + abs(xhilo[i] - yhilo[i]) + abs(xlolo[i] - ylolo[i]);

   return result;
}

double cmplx4_Difference_Sum
 ( int n, double *xrehihi, double *xrelohi, double *xrehilo, double *xrelolo,
   double *ximhihi, double *ximlohi, double *ximhilo, double *ximlolo,
   double *yrehihi, double *yrelohi, double *yrehilo, double *yrelolo,
   double *yimhihi, double *yimlohi, double *yimhilo, double *yimlolo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xrehihi[i] - yrehihi[i])
                      + abs(xrelohi[i] - yrelohi[i])
                      + abs(xrehilo[i] - yrehilo[i])
                      + abs(xrelolo[i] - yrelolo[i])
                      + abs(ximhihi[i] - yimhihi[i])
                      + abs(ximlohi[i] - yimlohi[i])
                      + abs(ximhilo[i] - yimhilo[i])
                      + abs(ximlolo[i] - yimlolo[i]);

   return result;
}

double dbl4_Column_Sum
 ( int dim, int col,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo )
{
   double resulthihi = 0.0;
   double resultlohi = 0.0;
   double resulthilo = 0.0;
   double resultlolo = 0.0;

   for(int i=0; i<dim; i++)
      qdf_inc(&resulthihi,&resultlohi,&resulthilo,&resultlolo,
              abs(Ahihi[i][col]),abs(Alohi[i][col]),
              abs(Ahilo[i][col]),abs(Alolo[i][col]));

   return resulthihi;
}

double cmplx4_Column_Sum
 ( int dim, int col,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo )
{
   double resultrehihi = 0.0;
   double resultrelohi = 0.0;
   double resultrehilo = 0.0;
   double resultrelolo = 0.0;
   double resultimhihi = 0.0;
   double resultimlohi = 0.0;
   double resultimhilo = 0.0;
   double resultimlolo = 0.0;

   for(int i=0; i<dim; i++)
   {
      qdf_inc(&resultrehihi,&resultrelohi,&resultrehilo,&resultrelolo,
              abs(Arehihi[i][col]),abs(Arelohi[i][col]),
              abs(Arehilo[i][col]),abs(Arelolo[i][col]));
      qdf_inc(&resultimhihi,&resultimlohi,&resultimhilo,&resultimlolo,
              abs(Aimhihi[i][col]),abs(Aimlohi[i][col]),
              abs(Aimhilo[i][col]),abs(Aimlolo[i][col]));
   }
   return (resultrehihi + resultimhihi);
}

double dbl4_Max_Column_Sum
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo )
{
   double result = dbl4_Column_Sum(dim,0,Ahihi,Alohi,Ahilo,Alolo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = dbl4_Column_Sum(dim,j,Ahihi,Alohi,Ahilo,Alolo);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double cmplx4_Max_Column_Sum
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo )
{
   double result = cmplx4_Column_Sum(dim,0,Arehihi,Arelohi,Arehilo,Arelolo,
                                           Aimhihi,Aimlohi,Aimhilo,Aimlolo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = cmplx4_Column_Sum(dim,j,Arehihi,Arelohi,Arehilo,Arelolo,
                                       Aimhihi,Aimlohi,Aimhilo,Aimlolo);

      if(colsum > result) result = colsum;
   }
   return result;  
}

double dbl4_condition
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **invAhihi, double **invAlohi,
   double **invAhilo, double **invAlolo )
{
   double Amaxcolsum = dbl4_Max_Column_Sum(dim,Ahihi,Alohi,Ahilo,Alolo);
   double invAmaxcolsum
      = dbl4_Max_Column_Sum(dim,invAhihi,invAlohi,invAhilo,invAlolo);

   return Amaxcolsum*invAmaxcolsum;
}

double cmplx4_condition
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **invArehihi, double **invArelohi,
   double **invArehilo, double **invArelolo,
   double **invAimhihi, double **invAimlohi,
   double **invAimhilo, double **invAimlolo )
{
   double Amaxcolsum
      = cmplx4_Max_Column_Sum(dim,Arehihi,Arelohi,Arehilo,Arelolo,
                                  Aimhihi,Aimlohi,Aimhilo,Aimlolo);
   double invAmaxcolsum
      = cmplx4_Max_Column_Sum(dim,
                              invArehihi,invArelohi,invArehilo,invArelolo,
                              invAimhihi,invAimlohi,invAimhilo,invAimlolo);

   return Amaxcolsum*invAmaxcolsum;
}

double dbl4_Matrix_Difference_Sum
 ( int n, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Ahihi[i][j] - Bhihi[i][j])
                         + abs(Alohi[i][j] - Blohi[i][j])
                         + abs(Ahilo[i][j] - Bhilo[i][j])
                         + abs(Alolo[i][j] - Blolo[i][j]);

   return result;
}

double cmplx4_Matrix_Difference_Sum
 ( int n,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Brehihi, double **Brelohi, double **Brehilo, double **Brelolo,
   double **Bimhihi, double **Bimlohi, double **Bimhilo, double **Bimlolo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Arehihi[i][j] - Brehihi[i][j])
                         + abs(Arelohi[i][j] - Brelohi[i][j])
                         + abs(Arehilo[i][j] - Brehilo[i][j])
                         + abs(Arelolo[i][j] - Brelolo[i][j])
                         + abs(Aimhihi[i][j] - Bimhihi[i][j])
                         + abs(Aimlohi[i][j] - Bimlohi[i][j])
                         + abs(Aimhilo[i][j] - Bimhilo[i][j])
                         + abs(Aimlolo[i][j] - Bimlolo[i][j]);

   return result;
}

double dbl4_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Ahihi, double **Alohi, double **Ahilo, double **Alolo,
   double **Bhihi, double **Blohi, double **Bhilo, double **Blolo )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Ahihi[offset+i][offset+j]
                                - Bhihi[offset+i][offset+j])
                            + abs(Alohi[offset+i][offset+j]
                                - Blohi[offset+i][offset+j])
                            + abs(Ahilo[offset+i][offset+j]
                                - Bhilo[offset+i][offset+j])
                            + abs(Alolo[offset+i][offset+j]
                                - Blolo[offset+i][offset+j]);
   }
   return result;
}

double cmplx4_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo,
   double **Brehihi, double **Brelohi, double **Brehilo, double **Brelolo,
   double **Bimhihi, double **Bimlohi, double **Bimhilo, double **Bimlolo )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Arehihi[offset+i][offset+j]
                                - Brehihi[offset+i][offset+j])
                            + abs(Arelohi[offset+i][offset+j]
                                - Brelohi[offset+i][offset+j])
                            + abs(Arehilo[offset+i][offset+j]
                                - Brehilo[offset+i][offset+j])
                            + abs(Arelolo[offset+i][offset+j]
                                - Brelolo[offset+i][offset+j])
                            + abs(Aimhihi[offset+i][offset+j]
                                - Bimhihi[offset+i][offset+j])
                            + abs(Aimlohi[offset+i][offset+j]
                                - Bimlohi[offset+i][offset+j])
                            + abs(Aimhilo[offset+i][offset+j]
                                - Bimhilo[offset+i][offset+j])
                            + abs(Aimlolo[offset+i][offset+j]
                                - Bimlolo[offset+i][offset+j]);
   }
   return result;
}

void dbl4_random_upper_factor
 ( int dim, double **Ahihi, double **Alohi, double **Ahilo, double **Alolo )
{
   random_dbl4_matrix(dim,dim,Ahihi,Alohi,Ahilo,Alolo);

   int *pivots = new int[dim];

   CPU_dbl4_factors_lufac(dim,Ahihi,Alohi,Ahilo,Alolo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Ahihi[i][j] = 0.0; Alohi[i][j] = 0.0;
         Ahilo[i][j] = 0.0; Alolo[i][j] = 0.0;
      }

   free(pivots);
}

void cmplx4_random_upper_factor
 ( int dim,
   double **Arehihi, double **Arelohi, double **Arehilo, double **Arelolo,
   double **Aimhihi, double **Aimlohi, double **Aimhilo, double **Aimlolo )
{
   random_cmplx4_matrix
      (dim,dim,Arehihi,Arelohi,Arehilo,Arelolo,
               Aimhihi,Aimlohi,Aimhilo,Aimlolo);

   int *pivots = new int[dim];

   CPU_cmplx4_factors_lufac(dim,Arehihi,Arelohi,Arehilo,Arelolo,
                                Aimhihi,Aimlohi,Aimhilo,Aimlolo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Arehihi[i][j] = 0.0; Arelohi[i][j] = 0.0;
         Arehilo[i][j] = 0.0; Arelolo[i][j] = 0.0;
         Aimhihi[i][j] = 0.0; Aimlohi[i][j] = 0.0;
         Aimhilo[i][j] = 0.0; Aimlolo[i][j] = 0.0;
      }

   free(pivots);
}
