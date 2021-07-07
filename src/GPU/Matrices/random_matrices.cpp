// The file random_matrices.cpp defines the functions with prototypes
// in random_matrices.h.

#include <cmath>
#include "random_numbers.h"
#include "random_series.h"
#include "random_matrices.h"

void random_dbl_series_vector
 ( int dim, int deg, double *x, double **v, bool expform )
{
   if(expform)
      for(int k=0; k<dim; k++) random_dbl_exponential(deg,&x[k],v[k]);
   else
      for(int k=0; k<dim; k++) random_dbl_logarithm(deg,&x[k],v[k]);
}

void random_cmplx_series_vector
 ( int dim, int deg, double *xre, double *xim,
   double **vre, double **vim, bool expform )
{
   if(expform)
      for(int k=0; k<dim; k++)
         random_cmplx_exponential(deg,&xre[k],&xim[k],vre[k],vim[k]);
   else
      for(int k=0; k<dim; k++)
         random_cmplx_logarithm(deg,&xre[k],&xim[k],vre[k],vim[k]);
}

void random_dbl_series_vectors
 ( int dim, int deg, double *x, double **plux, double **minx )
{
   for(int k=0; k<dim; k++)
      random_dbl_exponentials(deg,&x[k],plux[k],minx[k]);
}

void random_cmplx_series_vectors
 ( int dim, int deg, double *xre, double *xim,
   double **pluxre, double **pluxim, double **minxre, double **minxim )
{
   for(int k=0; k<dim; k++)
      random_cmplx_exponentials
         (deg,&xre[k],&xim[k],pluxre[k],pluxim[k],minxre[k],minxim[k]);
}

void random_dbl_series_matrix
 ( int rows, int cols, int deg, double **x, double ***A, bool expform )
{
   if(expform)
      for(int i=0; i<rows; i++)
         for(int j=0; j<cols; j++)
            random_dbl_exponential(deg,&x[i][j],A[i][j]);
   else
      for(int i=0; i<rows; i++)
         for(int j=0; j<cols; j++)
            random_dbl_logarithm(deg,&x[i][j],A[i][j]);
}

void random_dbl_upper_matrix ( int rows, int cols, double **A )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++) A[i][j] = 0.0;
      for(int j=i; j<cols; j++) A[i][j] = random_double();
   }
}

void random_dbl_matrix ( int rows, int cols, double **A )
{
   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++) A[i][j] = random_double();
}

void random_cmplx_upper_matrix
 ( int rows, int cols, double **Are, double **Aim )
{
   double rnd;

   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++) 
      {
         Are[i][j] = 0.0;
         Aim[i][j] = 0.0;
      }
      for(int j=i; j<cols; j++)
      {
         rnd = random_angle();
         Are[i][j] = cos(rnd);
         Aim[i][j] = sin(rnd);
      }
   }
}

void random_cmplx_matrix
 ( int rows, int cols, double **Are, double **Aim )
{
   double rnd;

   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
      {
         rnd = random_angle();
         Are[i][j] = cos(rnd);
         Aim[i][j] = sin(rnd);
      }
}

void random_dbl_upper_series_matrix
 ( int rows, int cols, int deg, double **x, double ***A, bool expform )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         for(int k=0; k<=deg; k++) A[i][j][k] = 0.0;
      }
      if(expform)
         for(int j=i; j<cols; j++)
            random_dbl_exponential(deg,&x[i][j],A[i][j]);
      else
         for(int j=i; j<cols; j++)
            random_dbl_logarithm(deg,&x[i][j],A[i][j]);
   }
}

void random_dbl_lower_series_matrix
 ( int rows, int cols, int deg, double **x, double ***A, bool expform )
{
   for(int i=0; i<rows; i++)
   {
      if(expform)
         for(int j=0; j<i; j++)
            random_dbl_exponential(deg,&x[i][j],A[i][j]);
      else
         for(int j=0; j<i; j++)
            random_dbl_logarithm(deg,&x[i][j],A[i][j]);
     
      x[i][i] = 0.0;
      A[i][i][0] = 1.0;
      for(int k=1; k<=deg; k++) A[i][i][k] = 0.0;

      for(int j=i+1; j<cols; j++)
      {
         for(int k=0; k<=deg; k++) A[i][j][k] = 0.0;
      }
   }
}

void random_cmplx_series_matrix
 ( int rows, int cols, int deg, double **xre, double **xim,
   double ***Are, double ***Aim, bool expform )
{
   for(int i=0; i<rows; i++)
      for(int j=0; j<cols; j++)
         if(expform)
            random_cmplx_exponential
               (deg,&xre[i][j],&xim[i][j],Are[i][j],Aim[i][j]);
         else
            random_cmplx_logarithm
               (deg,&xre[i][j],&xim[i][j],Are[i][j],Aim[i][j]);
}

void random_cmplx_upper_series_matrix
 ( int rows, int cols, int deg, double **xre, double **xim,
   double ***Are, double ***Aim, bool expform )
{
   for(int i=0; i<rows; i++)
   {
      for(int j=0; j<i; j++)
      {
         for(int k=0; k<=deg; k++)
         {
            Are[i][j][k] = 0.0;
            Aim[i][j][k] = 0.0;
         }
      }
      for(int j=i; j<cols; j++)
      {
         if(expform)
            random_cmplx_exponential
               (deg,&xre[i][j],&xim[i][j],Are[i][j],Aim[i][j]);
         else
            random_cmplx_logarithm
               (deg,&xre[i][j],&xim[i][j],Are[i][j],Aim[i][j]);
      }
   }
}

void random_cmplx_lower_series_matrix
 ( int rows, int cols, int deg, double **xre, double **xim,
   double ***Are, double ***Aim, bool expform )
{
   for(int i=0; i<rows; i++)
   {
      if(expform)
         for(int j=0; j<i; j++)
            random_cmplx_exponential
               (deg,&xre[i][j],&xim[i][j],Are[i][j],Aim[i][j]);

      else
         for(int j=0; j<i; j++)
            random_cmplx_logarithm
               (deg,&xre[i][j],&xim[i][j],Are[i][j],Aim[i][j]);

      xre[i][i] = 0.0; xim[i][i] = 0.0;
      Are[i][i][0] = 1.0; Aim[i][i][0] = 0.0;
      for(int k=1; k<=deg; k++)
      {
         Are[i][i][k] = 0.0;
         Aim[i][i][k] = 0.0;
      }
      for(int j=i+1; j<cols; j++)
      {
         for(int k=0; k<=deg; k++)
         {
            Are[i][j][k] = 0.0;
            Aim[i][j][k] = 0.0;
         }
      }
   }
}
