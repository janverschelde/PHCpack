// The file dbl4_newton_testers.cpp defines the functions with prototypes in
// the file dbl4_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "quad_double_functions.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"
#include "dbl4_factorizations.h"
#include "dbl4_bals_host.h"
#include "dbl4_bals_kernels.h"
#include "dbl4_tail_kernels.h"
#include "dbl4_systems_host.h"
#include "dbl4_systems_kernels.h"

using namespace std;

void real4_start_series_vector
 ( int dim, int deg,
   double *r0hihi, double *r0lohi, double *r0hilo, double *r0lolo,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihi[i][0] = r0hihi[i];
      cfflohi[i][0] = r0lohi[i];
      cffhilo[i][0] = r0hilo[i];
      cfflolo[i][0] = r0lolo[i];
      qdf_inc(&cffhihi[i][0],&cfflohi[i][0],&cffhilo[i][0],&cfflolo[i][0],
              1.0e-32,0.0,0.0,0.0);

      for(int j=1; j<=deg; j++)
      {
         cffhihi[i][j] = 0.0; cfflohi[i][j] = 0.0;
         cffhilo[i][j] = 0.0; cfflolo[i][j] = 0.0;
      }
   }
}

void cmplx4_start_series_vector
 ( int dim, int deg,
   double *r0rehihi, double *r0relohi, double *r0rehilo, double *r0relolo,
   double *r0imhihi, double *r0imlohi, double *r0imhilo, double *r0imlolo,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehihi[i][0] = r0rehihi[i];
      cffrelohi[i][0] = r0relohi[i];
      cffrehilo[i][0] = r0rehilo[i];
      cffrelolo[i][0] = r0relolo[i];
      qdf_inc(&cffrehihi[i][0],&cffrelohi[i][0],
              &cffrehilo[i][0],&cffrelolo[i][0],1.0e-32,0.0,0.0,0.0);

      cffimhihi[i][0] = r0imhihi[i];
      cffimlohi[i][0] = r0imlohi[i];
      cffimhilo[i][0] = r0imhilo[i];
      cffimlolo[i][0] = r0imlolo[i];
      qdf_inc(&cffimhihi[i][0],&cffimlohi[i][0],
              &cffimhilo[i][0],&cffimlolo[i][0],1.0e-32,0.0,0.0,0.0);

      for(int j=1; j<=deg; j++)
      {
         cffrehihi[i][j] = 0.0; cffrelohi[i][j] = 0.0;
         cffrehilo[i][j] = 0.0; cffrelolo[i][j] = 0.0;
         cffimhihi[i][j] = 0.0; cffimlohi[i][j] = 0.0;
         cffimhilo[i][j] = 0.0; cffimlolo[i][j] = 0.0;
      }
   }
}

void dbl4_unit_series_vector
 ( int dim, int deg,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihi[i][0] = 1.0;
      cfflohi[i][0] = 0.0;
      cffhilo[i][0] = 0.0;
      cfflolo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhihi[i][j] = 0.0;
         cfflohi[i][j] = 0.0;
         cffhilo[i][j] = 0.0;
         cfflolo[i][j] = 0.0;
      }
   }
}

void cmplx4_unit_series_vector
 ( int dim, int deg,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehihi[i][0] = 1.0; cffrelohi[i][0] = 0.0;
      cffrehilo[i][0] = 0.0; cffrelolo[i][0] = 0.0;
      cffimhihi[i][0] = 0.0; cffimlohi[i][0] = 0.0;
      cffimhilo[i][0] = 0.0; cffimlolo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffrehihi[i][j] = 0.0; cffrelohi[i][j] = 0.0;
         cffrehilo[i][j] = 0.0; cffrelolo[i][j] = 0.0;
         cffimhihi[i][j] = 0.0; cffimlohi[i][j] = 0.0;
         cffimhilo[i][j] = 0.0; cffimlolo[i][j] = 0.0;
      }
   }
}

void dbl4_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo )
{
   for(int i=0; i<nbr; i++)
      dbl4_unit_series_vector
         (dim,deg,cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i]);
}

void cmplx4_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo )
{
   for(int i=0; i<nbr; i++)
      cmplx4_unit_series_vector
         (dim,deg,cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
                  cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i]);
}

void dbl4_update_series
 ( int dim, int degp1, int startidx,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **dxhihi, double **dxlohi, double **dxhilo, double **dxlolo,
   int vrblvl )
{
   if(vrblvl > 1) cout << scientific << setprecision(16);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << xhihi[i][j] << "  " << xlohi[i][j] << endl;
            cout << xhilo[i][j] << "  " << xlolo[i][j] << endl;
         }
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         qdf_inc(&xhihi[i][j],&xlohi[i][j],&xhilo[i][j],&xlolo[i][j],
                 dxhihi[j][i],dxlohi[j][i],dxhilo[j][i],dxlolo[j][i]);
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << xhihi[i][j] << "  " << xlohi[i][j] << endl;
            cout << xhilo[i][j] << "  " << xlolo[i][j] << endl;
         }
      }
   }
}

void cmplx4_update_series
 ( int dim, int degp1, int startidx,
   double **xrehihi, double **xrelohi, double **xrehilo, double **xrelolo,
   double **ximhihi, double **ximlohi, double **ximhilo, double **ximlolo,
   double **dxrehihi, double **dxrelohi, double **dxrehilo, double **dxrelolo,
   double **dximhihi, double **dximlohi, double **dximhilo, double **dximlolo,
   int vrblvl )
{
   if(vrblvl > 1) cout << scientific << setprecision(16);

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehihi[i][j] << "  " << xrelohi[i][j] << endl << "  "
                 << xrehilo[i][j] << "  " << xrelolo[i][j] << endl << "  "
                 << ximhihi[i][j] << "  " << ximlohi[i][j] << endl << "  "
                 << ximhilo[i][j] << "  " << ximlolo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++)
      {
         // xre[i][j] = xre[i][j] + dxre[j][i];
         qdf_inc(&xrehihi[i][j],&xrelohi[i][j],&xrehilo[i][j],&xrelolo[i][j],
                 dxrehihi[j][i],dxrelohi[j][i],dxrehilo[j][i],dxrelolo[j][i]);
         // xim[i][j] = xim[i][j] + dxim[j][i];
         qdf_inc(&ximhihi[i][j],&ximlohi[i][j],&ximhilo[i][j],&ximlolo[i][j],
                 dximhihi[j][i],dximlohi[j][i],dximhilo[j][i],dximlolo[j][i]);
      }
 
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehihi[i][j] << "  " << xrelohi[i][j] << endl << "  "
                 << xrehilo[i][j] << "  " << xrelolo[i][j] << endl << "  "
                 << ximhihi[i][j] << "  " << ximlohi[i][j] << endl << "  "
                 << ximhilo[i][j] << "  " << ximlolo[i][j] << endl;
      }
   }
}

double dbl4_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datahihi_h, double ***datalohi_h,
   double ***datahilo_h, double ***datalolo_h,
   double ***datahihi_d, double ***datalohi_d,
   double ***datahilo_d, double ***datalolo_d,
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
                    << datahihi_h[k][i][j] << "  "
                    << datalohi_h[k][i][j] << endl << "  "
                    << datahilo_h[k][i][j] << "  "
                    << datalolo_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datahihi_d[k][i][j] << "  "
                    << datalohi_d[k][i][j] << endl << "  "
                    << datahilo_d[k][i][j] << "  "
                    << datalolo_d[k][i][j] << endl;
            }
            errsum += abs(datahihi_h[k][i][j] - datahihi_d[k][i][j])
                    + abs(datalohi_h[k][i][j] - datalohi_d[k][i][j])
                    + abs(datahilo_h[k][i][j] - datahilo_d[k][i][j])
                    + abs(datalolo_h[k][i][j] - datalolo_d[k][i][j]);
         }

   return errsum;
}

double cmplx4_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datarehihi_h, double ***datarelohi_h,
   double ***datarehilo_h, double ***datarelolo_h,
   double ***dataimhihi_h, double ***dataimlohi_h,
   double ***dataimhilo_h, double ***dataimlolo_h,
   double ***datarehihi_d, double ***datarelohi_d,
   double ***datarehilo_d, double ***datarelolo_d,
   double ***dataimhihi_d, double ***dataimlohi_d,
   double ***dataimhilo_d, double ***dataimlolo_d,
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
                    << datarehihi_h[k][i][j] << "  "
                    << datarelohi_h[k][i][j] << endl << "  "
                    << datarehilo_h[k][i][j] << "  "
                    << datarelolo_h[k][i][j] << endl << "  "
                    << dataimhihi_h[k][i][j] << "  "
                    << dataimlohi_h[k][i][j] << endl << "  "
                    << dataimhilo_h[k][i][j] << "  "
                    << dataimlolo_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datarehihi_d[k][i][j] << "  "
                    << datarelohi_d[k][i][j] << endl << "  "
                    << datarehilo_d[k][i][j] << "  "
                    << datarelolo_d[k][i][j] << endl << "  "
                    << dataimhihi_d[k][i][j] << "  "
                    << dataimlohi_d[k][i][j] << endl << "  "
                    << dataimhilo_d[k][i][j] << "  "
                    << dataimlolo_d[k][i][j] << endl;
            }
            errsum += abs(datarehihi_h[k][i][j] - datarehihi_d[k][i][j])
                    + abs(datarelohi_h[k][i][j] - datarelohi_d[k][i][j])
                    + abs(datarehilo_h[k][i][j] - datarehilo_d[k][i][j])
                    + abs(datarelolo_h[k][i][j] - datarelolo_d[k][i][j])
                    + abs(dataimhihi_h[k][i][j] - dataimhihi_d[k][i][j])
                    + abs(dataimlohi_h[k][i][j] - dataimlohi_d[k][i][j])
                    + abs(dataimhilo_h[k][i][j] - dataimhilo_d[k][i][j])
                    + abs(dataimlolo_h[k][i][j] - dataimlolo_d[k][i][j]);
         }

   return errsum;
}

double dbl4_error2sum
 ( int nrows, int ncols,
   double **datahihi_h, double **datalohi_h,
   double **datahilo_h, double **datalolo_h,
   double **datahihi_d, double **datalohi_d,
   double **datahilo_d, double **datalolo_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datahihi_h[i][j] << "  "
                 << datalohi_h[i][j] << endl << "  "
                 << datahilo_h[i][j] << "  "
                 << datalolo_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datahihi_d[i][j] << "  "
                 << datalohi_d[i][j] << endl << "  "
                 << datahilo_d[i][j] << "  "
                 << datalolo_d[i][j] << endl;
         }
         errsum += abs(datahihi_h[i][j] - datahihi_d[i][j])
                 + abs(datalohi_h[i][j] - datalohi_d[i][j])
                 + abs(datahilo_h[i][j] - datahilo_d[i][j])
                 + abs(datalolo_h[i][j] - datalolo_d[i][j]);
      }

   return errsum;
}

double cmplx4_error2sum
 ( int nrows, int ncols,
   double **datarehihi_h, double **datarelohi_h,
   double **datarehilo_h, double **datarelolo_h,
   double **dataimhihi_h, double **dataimlohi_h,
   double **dataimhilo_h, double **dataimlolo_h,
   double **datarehihi_d, double **datarelohi_d,
   double **datarehilo_d, double **datarelolo_d,
   double **dataimhihi_d, double **dataimlohi_d,
   double **dataimhilo_d, double **dataimlolo_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datarehihi_h[i][j] << "  "
                 << datarelohi_h[i][j] << endl << "  "
                 << datarehilo_h[i][j] << "  "
                 << datarelolo_h[i][j] << endl << "  "
                 << dataimhihi_h[i][j] << "  "
                 << dataimlohi_h[i][j] << endl << "  "
                 << dataimhilo_h[i][j] << "  "
                 << dataimlolo_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datarehihi_d[i][j] << "  "
                 << datarelohi_d[i][j] << endl << "  "
                 << datarehilo_d[i][j] << "  "
                 << datarelolo_d[i][j] << endl << "  "
                 << dataimhihi_d[i][j] << "  "
                 << dataimlohi_d[i][j] << endl << "  "
                 << dataimhilo_d[i][j] << "  "
                 << dataimlolo_d[i][j] << endl;
         }
         errsum += abs(datarehihi_h[i][j] - datarehihi_d[i][j])
                 + abs(datarelohi_h[i][j] - datarelohi_d[i][j])
                 + abs(datarehilo_h[i][j] - datarehilo_d[i][j])
                 + abs(datarelolo_h[i][j] - datarelolo_d[i][j])
                 + abs(dataimhihi_h[i][j] - dataimhihi_d[i][j])
                 + abs(dataimlohi_h[i][j] - dataimlohi_d[i][j])
                 + abs(dataimhilo_h[i][j] - dataimhilo_d[i][j])
                 + abs(dataimlolo_h[i][j] - dataimlolo_d[i][j]);
      }

   return errsum;
}
