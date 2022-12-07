// The file dbl8_newton_testers.cpp defines the functions with prototypes in
// the file dbl8_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "octo_double_functions.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"
#include "dbl8_factorizations.h"
#include "dbl8_bals_host.h"
#include "dbl8_bals_kernels.h"
#include "dbl8_tail_kernels.h"
#include "dbl8_systems_host.h"
#include "dbl8_systems_kernels.h"

using namespace std;

void real8_start_series_vector
 ( int dim, int deg,
   double *r0hihihi, double *r0lohihi, double *r0hilohi, double *r0lolohi,
   double *r0hihilo, double *r0lohilo, double *r0hilolo, double *r0lololo,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihihi[i][0] = r0hihihi[i]; cfflohihi[i][0] = r0lohihi[i];
      cffhilohi[i][0] = r0hilohi[i]; cfflolohi[i][0] = r0lolohi[i];
      cffhihilo[i][0] = r0hihilo[i]; cfflohilo[i][0] = r0lohilo[i];
      cffhilolo[i][0] = r0hilolo[i]; cfflololo[i][0] = r0lololo[i];

      odf_inc(&cffhihihi[i][0],&cfflohihi[i][0],
              &cffhilohi[i][0],&cfflolohi[i][0],
              &cffhihilo[i][0],&cfflohilo[i][0],
              &cffhilolo[i][0],&cfflololo[i][0],
              1.0e-64,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

      for(int j=1; j<=deg; j++)
      {
         cffhihihi[i][j] = 0.0; cfflohihi[i][j] = 0.0;
         cffhilohi[i][j] = 0.0; cfflolohi[i][j] = 0.0;
         cffhihilo[i][j] = 0.0; cfflohilo[i][j] = 0.0;
         cffhilolo[i][j] = 0.0; cfflololo[i][j] = 0.0;
      }
   }
}

void cmplx8_start_series_vector
 ( int dim, int deg,
   double *r0rehihihi, double *r0relohihi,
   double *r0rehilohi, double *r0relolohi,
   double *r0rehihilo, double *r0relohilo,
   double *r0rehilolo, double *r0relololo,
   double *r0imhihihi, double *r0imlohihi,
   double *r0imhilohi, double *r0imlolohi,
   double *r0imhihilo, double *r0imlohilo,
   double *r0imhilolo, double *r0imlololo,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehihihi[i][0] = r0rehihihi[i];
      cffrelohihi[i][0] = r0relohihi[i];
      cffrehilohi[i][0] = r0rehilohi[i];
      cffrelolohi[i][0] = r0relolohi[i];
      cffrehihilo[i][0] = r0rehihilo[i];
      cffrelohilo[i][0] = r0relohilo[i];
      cffrehilolo[i][0] = r0rehilolo[i];
      cffrelololo[i][0] = r0relololo[i];
      odf_inc(&cffrehihihi[i][0],&cffrelohihi[i][0],
              &cffrehilohi[i][0],&cffrelolohi[i][0],
              &cffrehihilo[i][0],&cffrelohilo[i][0],
              &cffrehilolo[i][0],&cffrelololo[i][0],
              1.0e-64,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

      cffimhihihi[i][0] = r0imhihihi[i];
      cffimlohihi[i][0] = r0imlohihi[i];
      cffimhilohi[i][0] = r0imhilohi[i];
      cffimlolohi[i][0] = r0imlolohi[i];
      cffimhihilo[i][0] = r0imhihilo[i];
      cffimlohilo[i][0] = r0imlohilo[i];
      cffimhilolo[i][0] = r0imhilolo[i];
      cffimlololo[i][0] = r0imlololo[i];
      odf_inc(&cffimhihihi[i][0],&cffimlohihi[i][0],
              &cffimhilohi[i][0],&cffimlolohi[i][0],
              &cffimhihilo[i][0],&cffimlohilo[i][0],
              &cffimhilolo[i][0],&cffimlololo[i][0],
              1.0e-64,0.0,0.0,0.0,0.0,0.0,0.0,0.0);

      for(int j=1; j<=deg; j++)
      {
         cffrehihihi[i][j] = 0.0; cffrelohihi[i][j] = 0.0;
         cffrehilohi[i][j] = 0.0; cffrelolohi[i][j] = 0.0;
         cffrehihilo[i][j] = 0.0; cffrelohilo[i][j] = 0.0;
         cffrehilolo[i][j] = 0.0; cffrelololo[i][j] = 0.0;
         cffimhihihi[i][j] = 0.0; cffimlohihi[i][j] = 0.0;
         cffimhilohi[i][j] = 0.0; cffimlolohi[i][j] = 0.0;
         cffimhihilo[i][j] = 0.0; cffimlohilo[i][j] = 0.0;
         cffimhilolo[i][j] = 0.0; cffimlololo[i][j] = 0.0;
      }
   }
}

void dbl8_unit_series_vector
 ( int dim, int deg,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihihi[i][0] = 1.0; cfflohihi[i][0] = 0.0;
      cffhilohi[i][0] = 0.0; cfflolohi[i][0] = 0.0;
      cffhihilo[i][0] = 0.0; cfflohilo[i][0] = 0.0;
      cffhilolo[i][0] = 0.0; cfflololo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhihihi[i][j] = 0.0; cfflohihi[i][j] = 0.0;
         cffhilohi[i][j] = 0.0; cfflolohi[i][j] = 0.0;
         cffhihilo[i][j] = 0.0; cfflohilo[i][j] = 0.0;
         cffhilolo[i][j] = 0.0; cfflololo[i][j] = 0.0;
      }
   }
}

void cmplx8_unit_series_vector
 ( int dim, int deg,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehihihi[i][0] = 1.0; cffrelohihi[i][0] = 0.0;
      cffrehilohi[i][0] = 0.0; cffrelolohi[i][0] = 0.0;
      cffrehihilo[i][0] = 0.0; cffrelohilo[i][0] = 0.0;
      cffrehilolo[i][0] = 0.0; cffrelololo[i][0] = 0.0;
      cffimhihihi[i][0] = 0.0; cffimlohihi[i][0] = 0.0;
      cffimhilohi[i][0] = 0.0; cffimlolohi[i][0] = 0.0;
      cffimhihilo[i][0] = 0.0; cffimlohilo[i][0] = 0.0;
      cffimhilolo[i][0] = 0.0; cffimlololo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffrehihihi[i][j] = 0.0; cffrelohihi[i][j] = 0.0;
         cffrehilohi[i][j] = 0.0; cffrelolohi[i][j] = 0.0;
         cffrehihilo[i][j] = 0.0; cffrelohilo[i][j] = 0.0;
         cffrehilolo[i][j] = 0.0; cffrelololo[i][j] = 0.0;
         cffimhihihi[i][j] = 0.0; cffimlohihi[i][j] = 0.0;
         cffimhilohi[i][j] = 0.0; cffimlolohi[i][j] = 0.0;
         cffimhihilo[i][j] = 0.0; cffimlohilo[i][j] = 0.0;
         cffimhilolo[i][j] = 0.0; cffimlololo[i][j] = 0.0;
      }
   }
}

void dbl8_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo )
{
   for(int i=0; i<nbr; i++)
      dbl8_unit_series_vector
         (dim,deg,cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
                  cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i]);
}

void cmplx8_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo )
{
   for(int i=0; i<nbr; i++)
      cmplx8_unit_series_vector
         (dim,deg,
          cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
          cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
          cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
          cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i]);
}

void dbl8_update_series
 ( int dim, int degp1, int startidx,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **dxhihihi, double **dxlohihi, double **dxhilohi, double **dxlolohi,
   double **dxhihilo, double **dxlohilo, double **dxhilolo, double **dxlololo,
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
            cout << xhihihi[i][j] << "  " << xlohihi[i][j] << endl;
            cout << xhilohi[i][j] << "  " << xlolohi[i][j] << endl;
            cout << xhihilo[i][j] << "  " << xlohilo[i][j] << endl;
            cout << xhilolo[i][j] << "  " << xlololo[i][j] << endl;
         }
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         odf_inc(&xhihihi[i][j],&xlohihi[i][j],&xhilohi[i][j],&xlolohi[i][j],
                 &xhihilo[i][j],&xlohilo[i][j],&xhilolo[i][j],&xlololo[i][j],
                 dxhihihi[j][i],dxlohihi[j][i],dxhilohi[j][i],dxlolohi[j][i],
                 dxhihilo[j][i],dxlohilo[j][i],dxhilolo[j][i],dxlololo[j][i]);
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
            cout << xhihihi[i][j] << "  " << xlohihi[i][j] << endl;
            cout << xhilohi[i][j] << "  " << xlolohi[i][j] << endl;
            cout << xhihilo[i][j] << "  " << xlohilo[i][j] << endl;
            cout << xhilolo[i][j] << "  " << xlololo[i][j] << endl;
         }
      }
   }
}

void cmplx8_update_series
 ( int dim, int degp1, int startidx,
   double **xrehihihi, double **xrelohihi,
   double **xrehilohi, double **xrelolohi,
   double **xrehihilo, double **xrelohilo,
   double **xrehilolo, double **xrelololo,
   double **ximhihihi, double **ximlohihi,
   double **ximhilohi, double **ximlolohi,
   double **ximhihilo, double **ximlohilo,
   double **ximhilolo, double **ximlololo,
   double **dxrehihihi, double **dxrelohihi,
   double **dxrehilohi, double **dxrelolohi,
   double **dxrehihilo, double **dxrelohilo,
   double **dxrehilolo, double **dxrelololo,
   double **dximhihihi, double **dximlohihi,
   double **dximhilohi, double **dximlolohi,
   double **dximhihilo, double **dximlohilo,
   double **dximhilolo, double **dximlololo, int vrblvl )
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
            cout << xrehihihi[i][j] << "  " << xrelohihi[i][j] << endl
                 << "  "
                 << xrehilohi[i][j] << "  " << xrelolohi[i][j] << endl
                 << "  "
                 << xrehihilo[i][j] << "  " << xrelohilo[i][j] << endl
                 << "  "
                 << xrehilolo[i][j] << "  " << xrelololo[i][j] << endl
                 << "  "
                 << ximhihihi[i][j] << "  " << ximlohihi[i][j] << endl
                 << "  "
                 << ximhilohi[i][j] << "  " << ximlolohi[i][j] << endl
                 << "  "
                 << ximhihilo[i][j] << "  " << ximlohilo[i][j] << endl
                 << "  "
                 << ximhilolo[i][j] << "  " << ximlololo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++)
      {
         // xre[i][j] = xre[i][j] + dxre[j][i];
         odf_inc(&xrehihihi[i][j],&xrelohihi[i][j],
                 &xrehilohi[i][j],&xrelolohi[i][j],
                 &xrehihilo[i][j],&xrelohilo[i][j],
                 &xrehilolo[i][j],&xrelololo[i][j],
                 dxrehihihi[j][i],dxrelohihi[j][i],
                 dxrehilohi[j][i],dxrelolohi[j][i],
                 dxrehihilo[j][i],dxrelohilo[j][i],
                 dxrehilolo[j][i],dxrelololo[j][i]);
         // xim[i][j] = xim[i][j] + dxim[j][i];
         odf_inc(&ximhihihi[i][j],&ximlohihi[i][j],
                 &ximhilohi[i][j],&ximlolohi[i][j],
                 &ximhihilo[i][j],&ximlohilo[i][j],
                 &ximhilolo[i][j],&ximlololo[i][j],
                 dximhihihi[j][i],dximlohihi[j][i],
                 dximhilohi[j][i],dximlolohi[j][i],
                 dximhihilo[j][i],dximlohilo[j][i],
                 dximhilolo[j][i],dximlololo[j][i]);
      }
 
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehihihi[i][j] << "  " << xrelohihi[i][j] << endl << "  "
                 << xrehilohi[i][j] << "  " << xrelolohi[i][j] << endl << "  "
                 << xrehihilo[i][j] << "  " << xrelohilo[i][j] << endl << "  "
                 << xrehilolo[i][j] << "  " << xrelololo[i][j] << endl << "  "
                 << ximhihihi[i][j] << "  " << ximlohihi[i][j] << endl << "  "
                 << ximhilohi[i][j] << "  " << ximlolohi[i][j] << endl << "  "
                 << ximhihilo[i][j] << "  " << ximlohilo[i][j] << endl << "  "
                 << ximhilolo[i][j] << "  " << ximlololo[i][j] << endl;
      }
   }
}

double dbl8_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datahihihi_h, double ***datalohihi_h,
   double ***datahilohi_h, double ***datalolohi_h,
   double ***datahihilo_h, double ***datalohilo_h,
   double ***datahilolo_h, double ***datalololo_h,
   double ***datahihihi_d, double ***datalohihi_d,
   double ***datahilohi_d, double ***datalolohi_d,
   double ***datahihilo_d, double ***datalohilo_d,
   double ***datahilolo_d, double ***datalololo_d,
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
                    << datahihihi_h[k][i][j] << "  "
                    << datalohihi_h[k][i][j] << endl << "  "
                    << datahilohi_h[k][i][j] << "  "
                    << datalolohi_h[k][i][j] << endl << "  "
                    << datahihilo_h[k][i][j] << "  "
                    << datalohilo_h[k][i][j] << endl << "  "
                    << datahilolo_h[k][i][j] << "  "
                    << datalololo_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datahihihi_d[k][i][j] << "  "
                    << datalohihi_d[k][i][j] << endl << "  "
                    << datahilohi_d[k][i][j] << "  "
                    << datalolohi_d[k][i][j] << endl << "  "
                    << datahihilo_d[k][i][j] << "  "
                    << datalohilo_d[k][i][j] << endl << "  "
                    << datahilolo_d[k][i][j] << "  "
                    << datalololo_d[k][i][j] << endl;
            }
            errsum += abs(datahihihi_h[k][i][j] - datahihihi_d[k][i][j])
                    + abs(datalohihi_h[k][i][j] - datalohihi_d[k][i][j])
                    + abs(datahilohi_h[k][i][j] - datahilohi_d[k][i][j])
                    + abs(datalolohi_h[k][i][j] - datalolohi_d[k][i][j])
                    + abs(datahihilo_h[k][i][j] - datahihilo_d[k][i][j])
                    + abs(datalohilo_h[k][i][j] - datalohilo_d[k][i][j])
                    + abs(datahilolo_h[k][i][j] - datahilolo_d[k][i][j])
                    + abs(datalololo_h[k][i][j] - datalololo_d[k][i][j]);
         }

   return errsum;
}

double cmplx8_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datarehihihi_h, double ***datarelohihi_h,
   double ***datarehilohi_h, double ***datarelolohi_h,
   double ***datarehihilo_h, double ***datarelohilo_h,
   double ***datarehilolo_h, double ***datarelololo_h,
   double ***dataimhihihi_h, double ***dataimlohihi_h,
   double ***dataimhilohi_h, double ***dataimlolohi_h,
   double ***dataimhihilo_h, double ***dataimlohilo_h,
   double ***dataimhilolo_h, double ***dataimlololo_h,
   double ***datarehihihi_d, double ***datarelohihi_d,
   double ***datarehilohi_d, double ***datarelolohi_d,
   double ***datarehihilo_d, double ***datarelohilo_d,
   double ***datarehilolo_d, double ***datarelololo_d,
   double ***dataimhihihi_d, double ***dataimlohihi_d,
   double ***dataimhilohi_d, double ***dataimlolohi_d,
   double ***dataimhihilo_d, double ***dataimlohilo_d,
   double ***dataimhilolo_d, double ***dataimlololo_d,
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
                    << datarehihihi_h[k][i][j] << "  "
                    << datarelohihi_h[k][i][j] << endl << "  "
                    << datarehilohi_h[k][i][j] << "  "
                    << datarelolohi_h[k][i][j] << endl << "  "
                    << datarehihilo_h[k][i][j] << "  "
                    << datarelohilo_h[k][i][j] << endl << "  "
                    << datarehilolo_h[k][i][j] << "  "
                    << datarelololo_h[k][i][j] << endl << "  "
                    << dataimhihihi_h[k][i][j] << "  "
                    << dataimlohihi_h[k][i][j] << endl << "  "
                    << dataimhilohi_h[k][i][j] << "  "
                    << dataimlolohi_h[k][i][j] << endl << "  "
                    << dataimhihilo_h[k][i][j] << "  "
                    << dataimlohilo_h[k][i][j] << endl << "  "
                    << dataimhilolo_h[k][i][j] << "  "
                    << dataimlololo_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datarehihihi_d[k][i][j] << "  "
                    << datarelohihi_d[k][i][j] << endl << "  "
                    << datarehilohi_d[k][i][j] << "  "
                    << datarelolohi_d[k][i][j] << endl << "  "
                    << datarehihilo_d[k][i][j] << "  "
                    << datarelohilo_d[k][i][j] << endl << "  "
                    << datarehilolo_d[k][i][j] << "  "
                    << datarelololo_d[k][i][j] << endl << "  "
                    << dataimhihihi_d[k][i][j] << "  "
                    << dataimlohihi_d[k][i][j] << endl << "  "
                    << dataimhilohi_d[k][i][j] << "  "
                    << dataimlolohi_d[k][i][j] << endl << "  "
                    << dataimhihilo_d[k][i][j] << "  "
                    << dataimlohilo_d[k][i][j] << endl << "  "
                    << dataimhilolo_d[k][i][j] << "  "
                    << dataimlololo_d[k][i][j] << endl;
            }
            errsum += abs(datarehihihi_h[k][i][j] - datarehihihi_d[k][i][j])
                    + abs(datarelohihi_h[k][i][j] - datarelohihi_d[k][i][j])
                    + abs(datarehilohi_h[k][i][j] - datarehilohi_d[k][i][j])
                    + abs(datarelolohi_h[k][i][j] - datarelolohi_d[k][i][j])
                    + abs(datarehihilo_h[k][i][j] - datarehihilo_d[k][i][j])
                    + abs(datarelohilo_h[k][i][j] - datarelohilo_d[k][i][j])
                    + abs(datarehilolo_h[k][i][j] - datarehilolo_d[k][i][j])
                    + abs(datarelololo_h[k][i][j] - datarelololo_d[k][i][j])
                    + abs(dataimhihihi_h[k][i][j] - dataimhihihi_d[k][i][j])
                    + abs(dataimlohihi_h[k][i][j] - dataimlohihi_d[k][i][j])
                    + abs(dataimhilohi_h[k][i][j] - dataimhilohi_d[k][i][j])
                    + abs(dataimlolohi_h[k][i][j] - dataimlolohi_d[k][i][j])
                    + abs(dataimhihilo_h[k][i][j] - dataimhihilo_d[k][i][j])
                    + abs(dataimlohilo_h[k][i][j] - dataimlohilo_d[k][i][j])
                    + abs(dataimhilolo_h[k][i][j] - dataimhilolo_d[k][i][j])
                    + abs(dataimlololo_h[k][i][j] - dataimlololo_d[k][i][j]);
         }

   return errsum;
}

double dbl8_error2sum
 ( int nrows, int ncols,
   double **datahihihi_h, double **datalohihi_h,
   double **datahilohi_h, double **datalolohi_h,
   double **datahihilo_h, double **datalohilo_h,
   double **datahilolo_h, double **datalololo_h,
   double **datahihihi_d, double **datalohihi_d,
   double **datahilohi_d, double **datalolohi_d,
   double **datahihilo_d, double **datalohilo_d,
   double **datahilolo_d, double **datalololo_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datahihihi_h[i][j] << "  "
                 << datalohihi_h[i][j] << endl << "  "
                 << datahilohi_h[i][j] << "  "
                 << datalolohi_h[i][j] << endl << "  "
                 << datahihilo_h[i][j] << "  "
                 << datalohilo_h[i][j] << endl << "  "
                 << datahilolo_h[i][j] << "  "
                 << datalololo_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datahihihi_d[i][j] << "  "
                 << datalohihi_d[i][j] << endl << "  "
                 << datahilohi_d[i][j] << "  "
                 << datalolohi_d[i][j] << endl << "  "
                 << datahihilo_d[i][j] << "  "
                 << datalohilo_d[i][j] << endl << "  "
                 << datahilolo_d[i][j] << "  "
                 << datalololo_d[i][j] << endl;
         }
         errsum += abs(datahihihi_h[i][j] - datahihihi_d[i][j])
                 + abs(datalohihi_h[i][j] - datalohihi_d[i][j])
                 + abs(datahilohi_h[i][j] - datahilohi_d[i][j])
                 + abs(datalolohi_h[i][j] - datalolohi_d[i][j])
                 + abs(datahihilo_h[i][j] - datahihilo_d[i][j])
                 + abs(datalohilo_h[i][j] - datalohilo_d[i][j])
                 + abs(datahilolo_h[i][j] - datahilolo_d[i][j])
                 + abs(datalololo_h[i][j] - datalololo_d[i][j]);
      }

   return errsum;
}

double cmplx8_error2sum
 ( int nrows, int ncols,
   double **datarehihihi_h, double **datarelohihi_h,
   double **datarehilohi_h, double **datarelolohi_h,
   double **datarehihilo_h, double **datarelohilo_h,
   double **datarehilolo_h, double **datarelololo_h,
   double **dataimhihihi_h, double **dataimlohihi_h,
   double **dataimhilohi_h, double **dataimlolohi_h,
   double **dataimhihilo_h, double **dataimlohilo_h,
   double **dataimhilolo_h, double **dataimlololo_h,
   double **datarehihihi_d, double **datarelohihi_d,
   double **datarehilohi_d, double **datarelolohi_d,
   double **datarehihilo_d, double **datarelohilo_d,
   double **datarehilolo_d, double **datarelololo_d,
   double **dataimhihihi_d, double **dataimlohihi_d,
   double **dataimhilohi_d, double **dataimlolohi_d,
   double **dataimhihilo_d, double **dataimlohilo_d,
   double **dataimhilolo_d, double **dataimlololo_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datarehihihi_h[i][j] << "  "
                 << datarelohihi_h[i][j] << endl << "  "
                 << datarehilohi_h[i][j] << "  "
                 << datarelolohi_h[i][j] << endl << "  "
                 << datarehihilo_h[i][j] << "  "
                 << datarelohilo_h[i][j] << endl << "  "
                 << datarehilolo_h[i][j] << "  "
                 << datarelololo_h[i][j] << endl << "  "
                 << dataimhihihi_h[i][j] << "  "
                 << dataimlohihi_h[i][j] << endl << "  "
                 << dataimhilohi_h[i][j] << "  "
                 << dataimlolohi_h[i][j] << endl << "  "
                 << dataimhihilo_h[i][j] << "  "
                 << dataimlohilo_h[i][j] << endl << "  "
                 << dataimhilolo_h[i][j] << "  "
                 << dataimlololo_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datarehihihi_d[i][j] << "  "
                 << datarelohihi_d[i][j] << endl << "  "
                 << datarehilohi_d[i][j] << "  "
                 << datarelolohi_d[i][j] << endl << "  "
                 << datarehihilo_d[i][j] << "  "
                 << datarelohilo_d[i][j] << endl << "  "
                 << datarehilolo_d[i][j] << "  "
                 << datarelololo_d[i][j] << endl << "  "
                 << dataimhihihi_d[i][j] << "  "
                 << dataimlohihi_d[i][j] << endl << "  "
                 << dataimhilohi_d[i][j] << "  "
                 << dataimlolohi_d[i][j] << endl << "  "
                 << dataimhihilo_d[i][j] << "  "
                 << dataimlohilo_d[i][j] << endl << "  "
                 << dataimhilolo_d[i][j] << "  "
                 << dataimlololo_d[i][j] << endl;
         }
         errsum += abs(datarehihihi_h[i][j] - datarehihihi_d[i][j])
                 + abs(datarelohihi_h[i][j] - datarelohihi_d[i][j])
                 + abs(datarehilohi_h[i][j] - datarehilohi_d[i][j])
                 + abs(datarelolohi_h[i][j] - datarelolohi_d[i][j])
                 + abs(datarehihilo_h[i][j] - datarehihilo_d[i][j])
                 + abs(datarelohilo_h[i][j] - datarelohilo_d[i][j])
                 + abs(datarehilolo_h[i][j] - datarehilolo_d[i][j])
                 + abs(datarelololo_h[i][j] - datarelololo_d[i][j])
                 + abs(dataimhihihi_h[i][j] - dataimhihihi_d[i][j])
                 + abs(dataimlohihi_h[i][j] - dataimlohihi_d[i][j])
                 + abs(dataimhilohi_h[i][j] - dataimhilohi_d[i][j])
                 + abs(dataimlolohi_h[i][j] - dataimlolohi_d[i][j])
                 + abs(dataimhihilo_h[i][j] - dataimhihilo_d[i][j])
                 + abs(dataimlohilo_h[i][j] - dataimlohilo_d[i][j])
                 + abs(dataimhilolo_h[i][j] - dataimhilolo_d[i][j])
                 + abs(dataimlololo_h[i][j] - dataimlololo_d[i][j]);
      }

   return errsum;
}
