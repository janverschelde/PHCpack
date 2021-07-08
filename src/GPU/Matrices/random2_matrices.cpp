// The file random2_matrices.cpp defines the functions specified in
// the file random2_matrices.h.

#include "random2_vectors.h"
#include "random2_matrices.h"
#include "double_double_functions.h"

void random_dbl2_upper_matrix
 ( int rows, int cols, double **Ahi, double **Alo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Ahi[i][j] = 0.0;
         Alo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
         random_double_double(&Ahi[i][j],&Alo[i][j]);
   }
}

void random_cmplx2_upper_matrix
 ( int rows, int cols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo )
{
   double rnd_hi,rnd_lo;
   double cosrnd_hi,cosrnd_lo,sinrnd_hi,sinrnd_lo;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         Arehi[i][j] = 0.0; Arelo[i][j] = 0.0;
         Aimhi[i][j] = 0.0; Aimlo[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
      {
         random_double_double(&rnd_hi,&rnd_lo);
         sinrnd_hi = rnd_hi; sinrnd_lo = rnd_lo;

         double y_hi,y_lo;          // work around to compute cos

         ddf_sqr(sinrnd_hi,sinrnd_lo,&y_hi,&y_lo);
         ddf_minus(&y_hi,&y_lo);                    // y = -sin^2
         ddf_inc_d(&y_hi,&y_lo,1.0);                // y = 1 - sin^2
         ddf_sqrt(y_hi,y_lo,&cosrnd_hi,&cosrnd_lo); // cos is sqrt(1-sin^2)

         Arehi[i][j] = cosrnd_hi; Arelo[i][j] = cosrnd_lo;
         Aimhi[i][j] = sinrnd_hi; Aimlo[i][j] = sinrnd_lo;
      }
   }
}

void random_dbl2_matrix
 ( int rows, int cols, double **Ahi, double **Alo )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
         random_double_double(&Ahi[i][j],&Alo[i][j]);
   }
}

void random_cmplx2_matrix
 ( int rows, int cols,
   double **Arehi, double **Arelo, double **Aimhi, double **Aimlo )
{
   double rnd_hi,rnd_lo;
   double cosrnd_hi,cosrnd_lo,sinrnd_hi,sinrnd_lo;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<cols; j++)
      {
         random_double_double(&rnd_hi,&rnd_lo);
         sinrnd_hi = rnd_hi; sinrnd_lo = rnd_lo;

         double y_hi,y_lo;          // work around to compute cos

         ddf_sqr(sinrnd_hi,sinrnd_lo,&y_hi,&y_lo);
         ddf_minus(&y_hi,&y_lo);                    // y = -sin^2
         ddf_inc_d(&y_hi,&y_lo,1.0);                // y = 1 - sin^2
         ddf_sqrt(y_hi,y_lo,&cosrnd_hi,&cosrnd_lo); // cos is sqrt(1-sin^2)

         Arehi[i][j] = cosrnd_hi; Arelo[i][j] = cosrnd_lo;
         Aimhi[i][j] = sinrnd_hi; Aimlo[i][j] = sinrnd_lo;
      }
   }
}
