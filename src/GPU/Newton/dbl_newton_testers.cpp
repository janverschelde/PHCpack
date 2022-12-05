// The file dbl_newton_testers.cpp defines the functions with prototypes in
// the file dbl_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void real_start_series_vector ( int dim, int deg, double *r0, double **cff )
{
   double angle;

   for(int i=0; i<dim; i++)
   {
      cff[i][0] = r0[i] + 1.0e-8; // random_double(); => no convergence ...

      for(int j=1; j<=deg; j++) cff[i][j] = 0.0;
   }
}

void cmplx_start_series_vector
 ( int dim, int deg, double *r0re, double *r0im,
   double **cffre, double **cffim )
{
   double angle;

   for(int i=0; i<dim; i++)
   {
      // angle = random_angle(); 
      cffre[i][0] = r0re[i] + 1.0e-8; // cos(angle); => no convergence ...
      cffim[i][0] = r0im[i] + 1.0e-8; // sin(angle);

      for(int j=1; j<=deg; j++)
      {
         cffre[i][j] = 0.0;
         cffim[i][j] = 0.0;
      }
   }
}

void dbl_unit_series_vector ( int dim, int deg, double **cff )
{
   double angle;

   for(int i=0; i<dim; i++)
   {
      cff[i][0] = 1.0;
      for(int j=1; j<=deg; j++) cff[i][j] = 0.0;
   }
}

void cmplx_unit_series_vector
 ( int dim, int deg, double **cffre, double **cffim )
{
   for(int i=0; i<dim; i++)
   {
      cffre[i][0] = 1.0;
      cffim[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffre[i][j] = 0.0;
         cffim[i][j] = 0.0;
      }
   }
}

void dbl_unit_series_vectors ( int nbr, int dim, int deg, double ***cff )
{
   for(int i=0; i<nbr; i++)
      dbl_unit_series_vector(dim,deg,cff[i]);
}

void cmplx_unit_series_vectors
 ( int nbr, int dim, int deg, double ***cffre, double ***cffim )
{
   for(int i=0; i<nbr; i++)
      cmplx_unit_series_vector(dim,deg,cffre[i],cffim[i]);
}

void dbl_update_series
 ( int dim, int degp1, int startidx, double **x, double **dx, int vrblvl )
{
   if(vrblvl > 1) cout << scientific << setprecision(16);

   if(vrblvl > 1)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++) cout << x[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++) x[i][j] = x[i][j] + dx[j][i];

   if(vrblvl > 1)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++) cout << x[i][j] << endl;
      }
   }
}

void cmplx_update_series
 ( int dim, int degp1, int startidx,
   double **xre, double **xim, double **dxre, double **dxim, int vrblvl )
{
   if(vrblvl > 1) cout << scientific << setprecision(16);

   if(vrblvl > 1)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xre[i][j] << "  " << xim[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++)
      {
         xre[i][j] = xre[i][j] + dxre[j][i];
         xim[i][j] = xim[i][j] + dxim[j][i];
      }

   if(vrblvl > 1)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xre[i][j] << "  " << xim[i][j] << endl;
      }
   }
}

double dbl_error3sum
 ( int dim1, int dim2, int dim3,
   double ***data_h, double ***data_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int k=0; k<dim1; k++) // monomial k
      for(int i=0; i<dim2; i++)
         for(int j=0; j<dim3; j++)
         {
            if(vrblvl > 1)
            {
               cout << banner << "_h[" // "output_h["
                    << k << "][" << i << "][" << j << "] : "
                    << data_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << data_d[k][i][j] << endl;
            }
            errsum += abs(data_h[k][i][j] - data_d[k][i][j]);
         }

   return errsum;
}

double cmplx_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datare_h, double ***dataim_h,
   double ***datare_d, double ***dataim_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int k=0; k<dim1; k++) // monomial k
      for(int i=0; i<dim2; i++)
         for(int j=0; j<dim3; j++)
         {
            if(vrblvl > 1)
            {
               cout << banner << "_h[" // "output_h["
                    << k << "][" << i << "][" << j << "] : "
                    << datare_h[k][i][j] << "  "
                    << dataim_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datare_d[k][i][j] << "  "
                    << dataim_d[k][i][j] << endl;
            }
            errsum += abs(datare_h[k][i][j] - datare_d[k][i][j])
                    + abs(dataim_h[k][i][j] - dataim_d[k][i][j]);
         }

   return errsum;
}

double dbl_error2sum
 ( int nrows, int ncols, double **data_h, double **data_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << data_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << data_d[i][j] << endl;
         }
         errsum += abs(data_h[i][j] - data_d[i][j]);
      }

   return errsum;
}

double cmplx_error2sum
 ( int nrows, int ncols,
   double **datare_h, double **dataim_h,
   double **datare_d, double **dataim_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datare_h[i][j] << "  "
                 << dataim_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datare_d[i][j] << "  "
                 << dataim_d[i][j] << endl;
         }
         errsum += abs(datare_h[i][j] - datare_d[i][j])
                 + abs(dataim_h[i][j] - dataim_d[i][j]);
      }

   return errsum;
}
