// The file dbl8_test_utilities.cpp defines the functions specified in
// the file dbl8_test_utilities.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include "octo_double_functions.h"
#include "random8_matrices.h"
#include "dbl8_factorizations.h"
#include "dbl8_test_utilities.h"

using namespace std;

double dbl8_Difference_Sum
 ( int n,
   double *xhihihi, double *xlohihi, double *xhilohi, double *xlolohi,
   double *xhihilo, double *xlohilo, double *xhilolo, double *xlololo,
   double *yhihihi, double *ylohihi, double *yhilohi, double *ylolohi,
   double *yhihilo, double *ylohilo, double *yhilolo, double *ylololo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xhihihi[i] - yhihihi[i])
                      + abs(xlohihi[i] - ylohihi[i])
                      + abs(xhilohi[i] - yhilohi[i])
                      + abs(xlolohi[i] - ylolohi[i])
                      + abs(xhihilo[i] - yhihilo[i])
                      + abs(xlohilo[i] - ylohilo[i])
                      + abs(xhilolo[i] - yhilolo[i])
                      + abs(xlololo[i] - ylololo[i]);

   return result;
}

double cmplx8_Difference_Sum
 ( int n,
   double *xrehihihi, double *xrelohihi,
   double *xrehilohi, double *xrelolohi,
   double *xrehihilo, double *xrelohilo,
   double *xrehilolo, double *xrelololo,
   double *ximhihihi, double *ximlohihi,
   double *ximhilohi, double *ximlolohi,
   double *ximhihilo, double *ximlohilo,
   double *ximhilolo, double *ximlololo,
   double *yrehihihi, double *yrelohihi,
   double *yrehilohi, double *yrelolohi,
   double *yrehihilo, double *yrelohilo,
   double *yrehilolo, double *yrelololo,
   double *yimhihihi, double *yimlohihi,
   double *yimhilohi, double *yimlolohi,
   double *yimhihilo, double *yimlohilo,
   double *yimhilolo, double *yimlololo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      result = result + abs(xrehihihi[i] - yrehihihi[i])
                      + abs(xrelohihi[i] - yrelohihi[i])
                      + abs(xrehilohi[i] - yrehilohi[i])
                      + abs(xrelolohi[i] - yrelolohi[i])
                      + abs(xrehihilo[i] - yrehihilo[i])
                      + abs(xrelohilo[i] - yrelohilo[i])
                      + abs(xrehilolo[i] - yrehilolo[i])
                      + abs(xrelololo[i] - yrelololo[i])
                      + abs(ximhihihi[i] - yimhihihi[i])
                      + abs(ximlohihi[i] - yimlohihi[i])
                      + abs(ximhilohi[i] - yimhilohi[i])
                      + abs(ximlolohi[i] - yimlolohi[i])
                      + abs(ximhihilo[i] - yimhihilo[i])
                      + abs(ximlohilo[i] - yimlohilo[i])
                      + abs(ximhilolo[i] - yimhilolo[i])
                      + abs(ximlololo[i] - yimlololo[i]);

   return result;
}

double dbl8_Column_Sum
 ( int dim, int col,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo )
{
   double resulthihihi = 0.0;
   double resultlohihi = 0.0;
   double resulthilohi = 0.0;
   double resultlolohi = 0.0;
   double resulthihilo = 0.0;
   double resultlohilo = 0.0;
   double resulthilolo = 0.0;
   double resultlololo = 0.0;

   for(int i=0; i<dim; i++)
      odf_inc(&resulthihihi,&resultlohihi,&resulthilohi,&resultlolohi,
              &resulthihilo,&resultlohilo,&resulthilolo,&resultlololo,
              abs(Ahihihi[i][col]),abs(Alohihi[i][col]),
              abs(Ahilohi[i][col]),abs(Alolohi[i][col]),
              abs(Ahihilo[i][col]),abs(Alohilo[i][col]),
              abs(Ahilolo[i][col]),abs(Alololo[i][col]));

   return resulthihihi;
}

double cmplx8_Column_Sum
 ( int dim, int col,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo )
{
   double resultrehihihi = 0.0;
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

   for(int i=0; i<dim; i++)
   {
      odf_inc(&resultrehihihi,&resultrelohihi,&resultrehilohi,&resultrelolohi,
              &resultrehihilo,&resultrelohilo,&resultrehilolo,&resultrelololo,
              abs(Arehihihi[i][col]),abs(Arelohihi[i][col]),
              abs(Arehilohi[i][col]),abs(Arelolohi[i][col]),
              abs(Arehihilo[i][col]),abs(Arelohilo[i][col]),
              abs(Arehilolo[i][col]),abs(Arelololo[i][col]));
      odf_inc(&resultimhihihi,&resultimlohihi,&resultimhilohi,&resultimlolohi,
              &resultimhihilo,&resultimlohilo,&resultimhilolo,&resultimlololo,
              abs(Aimhihihi[i][col]),abs(Aimlohihi[i][col]),
              abs(Aimhilohi[i][col]),abs(Aimlolohi[i][col]),
              abs(Aimhihilo[i][col]),abs(Aimlohilo[i][col]),
              abs(Aimhilolo[i][col]),abs(Aimlololo[i][col]));
   }
   return (resultrehihihi + resultimhihihi);
}

double dbl8_Max_Column_Sum
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo )
{
   double result = dbl8_Column_Sum(dim,0,Ahihihi,Alohihi,Ahilohi,Alolohi,
                                         Ahihilo,Alohilo,Ahilolo,Alololo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = dbl8_Column_Sum(dim,j,Ahihihi,Alohihi,Ahilohi,Alolohi,
                                     Ahihilo,Alohilo,Ahilolo,Alololo);
      if(colsum > result) result = colsum;
   }
   return result;  
}

double cmplx8_Max_Column_Sum
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo )
{
   double result = cmplx8_Column_Sum
                      (dim,0,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                             Arehihilo,Arelohilo,Arehilolo,Arelololo,
                             Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                             Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);
   double colsum;
   
   for(int j=1; j<dim; j++)
   {
      colsum = cmplx8_Column_Sum
                  (dim,j,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                         Arehihilo,Arelohilo,Arehilolo,Arelololo,
                         Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                         Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);

      if(colsum > result) result = colsum;
   }
   return result;  
}

double dbl8_condition
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **invAhihihi, double **invAlohihi,
   double **invAhilohi, double **invAlolohi,
   double **invAhihilo, double **invAlohilo,
   double **invAhilolo, double **invAlololo )
{
   double Amaxcolsum = dbl8_Max_Column_Sum
                           (dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
                                Ahihilo,Alohilo,Ahilolo,Alololo);
   double invAmaxcolsum
      = dbl8_Max_Column_Sum(dim,invAhihihi,invAlohihi,invAhilohi,invAlolohi,
                                invAhihilo,invAlohilo,invAhilolo,invAlololo);

   return Amaxcolsum*invAmaxcolsum;
}

double cmplx8_condition
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **invArehihihi, double **invArelohihi,
   double **invArehilohi, double **invArelolohi,
   double **invArehihilo, double **invArelohilo,
   double **invArehilolo, double **invArelololo,
   double **invAimhihihi, double **invAimlohihi,
   double **invAimhilohi, double **invAimlolohi,
   double **invAimhihilo, double **invAimlohilo,
   double **invAimhilolo, double **invAimlololo )
{
   double Amaxcolsum
      = cmplx8_Max_Column_Sum(dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
                                  Arehihilo,Arelohilo,Arehilolo,Arelololo,
                                  Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
                                  Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);
   double invAmaxcolsum
      = cmplx8_Max_Column_Sum
           (dim,invArehihihi,invArelohihi,invArehilohi,invArelolohi,
                invArehihilo,invArelohilo,invArehilolo,invArelololo,
                invAimhihihi,invAimlohihi,invAimhilohi,invAimlolohi,
                invAimhihilo,invAimlohilo,invAimhilolo,invAimlololo);

   return Amaxcolsum*invAmaxcolsum;
}

double dbl8_Matrix_Difference_Sum
 ( int n,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Bhihihi, double **Blohihi, double **Bhilohi, double **Blolohi,
   double **Bhihilo, double **Blohilo, double **Bhilolo, double **Blololo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Ahihihi[i][j] - Bhihihi[i][j])
                         + abs(Alohihi[i][j] - Blohihi[i][j])
                         + abs(Ahilohi[i][j] - Bhilohi[i][j])
                         + abs(Alolohi[i][j] - Blolohi[i][j])
                         + abs(Ahihilo[i][j] - Bhihilo[i][j])
                         + abs(Alohilo[i][j] - Blohilo[i][j])
                         + abs(Ahilolo[i][j] - Bhilolo[i][j])
                         + abs(Alololo[i][j] - Blololo[i][j]);

   return result;
}

double cmplx8_Matrix_Difference_Sum
 ( int n,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Brehihihi, double **Brelohihi,
   double **Brehilohi, double **Brelolohi,
   double **Brehihilo, double **Brelohilo,
   double **Brehilolo, double **Brelololo,
   double **Bimhihihi, double **Bimlohihi,
   double **Bimhilohi, double **Bimlolohi,
   double **Bimhihilo, double **Bimlohilo,
   double **Bimhilolo, double **Bimlololo )
{
   double result = 0.0;

   for(int i=0; i<n; i++)
      for(int j=0; j<n; j++)
         result = result + abs(Arehihihi[i][j] - Brehihihi[i][j])
                         + abs(Arelohihi[i][j] - Brelohihi[i][j])
                         + abs(Arehilohi[i][j] - Brehilohi[i][j])
                         + abs(Arelolohi[i][j] - Brelolohi[i][j])
                         + abs(Arehihilo[i][j] - Brehihilo[i][j])
                         + abs(Arelohilo[i][j] - Brelohilo[i][j])
                         + abs(Arehilolo[i][j] - Brehilolo[i][j])
                         + abs(Arelololo[i][j] - Brelololo[i][j])
                         + abs(Aimhihihi[i][j] - Bimhihihi[i][j])
                         + abs(Aimlohihi[i][j] - Bimlohihi[i][j])
                         + abs(Aimhilohi[i][j] - Bimhilohi[i][j])
                         + abs(Aimlolohi[i][j] - Bimlolohi[i][j])
                         + abs(Aimhihilo[i][j] - Bimhihilo[i][j])
                         + abs(Aimlohilo[i][j] - Bimlohilo[i][j])
                         + abs(Aimhilolo[i][j] - Bimhilolo[i][j])
                         + abs(Aimlololo[i][j] - Bimlololo[i][j]);

   return result;
}

double dbl8_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo,
   double **Bhihihi, double **Blohihi, double **Bhilohi, double **Blolohi,
   double **Bhihilo, double **Blohilo, double **Bhilolo, double **Blololo )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Ahihihi[offset+i][offset+j]
                                - Bhihihi[offset+i][offset+j])
                            + abs(Alohihi[offset+i][offset+j]
                                - Blohihi[offset+i][offset+j])
                            + abs(Ahilohi[offset+i][offset+j]
                                - Bhilohi[offset+i][offset+j])
                            + abs(Alolohi[offset+i][offset+j]
                                - Blolohi[offset+i][offset+j])
                            + abs(Ahihilo[offset+i][offset+j]
                                - Bhihilo[offset+i][offset+j])
                            + abs(Alohilo[offset+i][offset+j]
                                - Blohilo[offset+i][offset+j])
                            + abs(Ahilolo[offset+i][offset+j]
                                - Bhilolo[offset+i][offset+j])
                            + abs(Alololo[offset+i][offset+j]
                                - Blololo[offset+i][offset+j]);
   }
   return result;
}

double cmplx8_Diagonal_Difference_Sum
 ( int nbt, int szt,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo,
   double **Brehihihi, double **Brelohihi,
   double **Brehilohi, double **Brelolohi,
   double **Brehihilo, double **Brelohilo,
   double **Brehilolo, double **Brelololo,
   double **Bimhihihi, double **Bimlohihi,
   double **Bimhilohi, double **Bimlolohi,
   double **Bimhihilo, double **Bimlohilo,
   double **Bimhilolo, double **Bimlololo )
{
   double result = 0.0;
   int offset;

   for(int k=0; k<nbt; k++) // difference between k-th tiles
   {
      offset = k*szt;
      for(int i=0; i<szt; i++)
         for(int j=0; j<szt; j++)
            result = result + abs(Arehihihi[offset+i][offset+j]
                                - Brehihihi[offset+i][offset+j])
                            + abs(Arelohihi[offset+i][offset+j]
                                - Brelohihi[offset+i][offset+j])
                            + abs(Arehilohi[offset+i][offset+j]
                                - Brehilohi[offset+i][offset+j])
                            + abs(Arelolohi[offset+i][offset+j]
                                - Brelolohi[offset+i][offset+j])
                            + abs(Arehihilo[offset+i][offset+j]
                                - Brehihilo[offset+i][offset+j])
                            + abs(Arelohilo[offset+i][offset+j]
                                - Brelohilo[offset+i][offset+j])
                            + abs(Arehilolo[offset+i][offset+j]
                                - Brehilolo[offset+i][offset+j])
                            + abs(Arelololo[offset+i][offset+j]
                                - Brelololo[offset+i][offset+j])
                            + abs(Aimhihihi[offset+i][offset+j]
                                - Bimhihihi[offset+i][offset+j])
                            + abs(Aimlohihi[offset+i][offset+j]
                                - Bimlohihi[offset+i][offset+j])
                            + abs(Aimhilohi[offset+i][offset+j]
                                - Bimhilohi[offset+i][offset+j])
                            + abs(Aimlolohi[offset+i][offset+j]
                                - Bimlolohi[offset+i][offset+j])
                            + abs(Aimhihilo[offset+i][offset+j]
                                - Bimhihilo[offset+i][offset+j])
                            + abs(Aimlohilo[offset+i][offset+j]
                                - Bimlohilo[offset+i][offset+j])
                            + abs(Aimhilolo[offset+i][offset+j]
                                - Bimhilolo[offset+i][offset+j])
                            + abs(Aimlololo[offset+i][offset+j]
                                - Bimlololo[offset+i][offset+j]);
   }
   return result;
}

void dbl8_random_upper_factor
 ( int dim,
   double **Ahihihi, double **Alohihi, double **Ahilohi, double **Alolohi,
   double **Ahihilo, double **Alohilo, double **Ahilolo, double **Alololo )
{
   random_dbl8_matrix(dim,dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
                              Ahihilo,Alohilo,Ahilolo,Alololo);

   int *pivots = new int[dim];

   CPU_dbl8_factors_lufac(dim,Ahihihi,Alohihi,Ahilohi,Alolohi,
                              Ahihilo,Alohilo,Ahilolo,Alololo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Ahihihi[i][j] = 0.0; Alohihi[i][j] = 0.0;
         Ahilohi[i][j] = 0.0; Alolohi[i][j] = 0.0;
         Ahihilo[i][j] = 0.0; Alohilo[i][j] = 0.0;
         Ahilolo[i][j] = 0.0; Alololo[i][j] = 0.0;
      }

   free(pivots);
}

void cmplx8_random_upper_factor
 ( int dim,
   double **Arehihihi, double **Arelohihi,
   double **Arehilohi, double **Arelolohi,
   double **Arehihilo, double **Arelohilo,
   double **Arehilolo, double **Arelololo,
   double **Aimhihihi, double **Aimlohihi,
   double **Aimhilohi, double **Aimlolohi,
   double **Aimhihilo, double **Aimlohilo,
   double **Aimhilolo, double **Aimlololo )
{
   random_cmplx8_matrix
      (dim,dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
               Arehihilo,Arelohilo,Arehilolo,Arelololo,
               Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
               Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo);

   int *pivots = new int[dim];

   CPU_cmplx8_factors_lufac
      (dim,Arehihihi,Arelohihi,Arehilohi,Arelolohi,
           Arehihilo,Arelohilo,Arehilolo,Arelololo,
           Aimhihihi,Aimlohihi,Aimhilohi,Aimlolohi,
           Aimhihilo,Aimlohilo,Aimhilolo,Aimlololo,pivots);

   for(int i=0; i<dim; i++)
      for(int j=0; j<i; j++)
      {
         Arehihihi[i][j] = 0.0; Arelohihi[i][j] = 0.0;
         Arehilohi[i][j] = 0.0; Arelolohi[i][j] = 0.0;
         Arehihilo[i][j] = 0.0; Arelohilo[i][j] = 0.0;
         Arehilolo[i][j] = 0.0; Arelololo[i][j] = 0.0;
         Aimhihihi[i][j] = 0.0; Aimlohihi[i][j] = 0.0;
         Aimhilohi[i][j] = 0.0; Aimlolohi[i][j] = 0.0;
         Aimhihilo[i][j] = 0.0; Aimlohilo[i][j] = 0.0;
         Aimhilolo[i][j] = 0.0; Aimlololo[i][j] = 0.0;
      }

   free(pivots);
}
