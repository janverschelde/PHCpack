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

void dbl8_start_series_vector
 ( int dim, int deg,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihihi[i][0] = 1.000005; cfflohihi[i][0] = 0.0;
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

void cmplx8_start_series_vector
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
      cffrehihihi[i][0] = 1.000005; cffrelohihi[i][0] = 0.0;
      cffrehilohi[i][0] = 0.0; cffrelolohi[i][0] = 0.0;
      cffrehihilo[i][0] = 0.0; cffrelohilo[i][0] = 0.0;
      cffrehilolo[i][0] = 0.0; cffrelololo[i][0] = 0.0;
      cffimhihihi[i][0] = 0.000005; cffimlohihi[i][0] = 0.0;
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

void dbl8_update_series
 ( int dim, int degp1,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **dxhihihi, double **dxlohihi, double **dxhilohi, double **dxlolohi,
   double **dxhihilo, double **dxlohilo, double **dxhilolo, double **dxlololo,
   int vrblvl )
{
   if(vrblvl > 0) cout << scientific << setprecision(16);

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
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         odf_inc(&xhihihi[i][j],&xlohihi[i][j],&xhilohi[i][j],&xlolohi[i][j],
                 &xhihilo[i][j],&xlohilo[i][j],&xhilolo[i][j],&xlololo[i][j],
                 dxhihihi[j][i],dxlohihi[j][i],dxhilohi[j][i],dxlolohi[j][i],
                 dxhihilo[j][i],dxlohilo[j][i],dxhilolo[j][i],dxlololo[j][i]);
      }

   if(vrblvl > 0)
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
 ( int dim, int degp1,
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
   if(vrblvl > 0) cout << scientific << setprecision(16);

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
   for(int j=0; j<degp1; j++) 
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
 
   if(vrblvl > 0)
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

void dbl8_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **workrhshihihi, double **workrhslohihi,
   double **workrhshilohi, double **workrhslolohi,
   double **workrhshihilo, double **workrhslohilo,
   double **workrhshilolo, double **workrhslololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int *ipvt, int vrblvl )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   dbl8_unit_series_vector
      (dim,deg,cffhihihi,cfflohihi,cffhilohi,cfflolohi,
               cffhihilo,cfflohilo,cffhilolo,cfflololo);

   CPU_dbl8_evaluate_monomials
      (dim,deg,nvr,idx,exp,nbrfac,expfac,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       acchihihi,acclohihi,acchilohi,acclolohi,
       acchihilo,acclohilo,acchilolo,acclololo,
       inputhihihi, inputlohihi, inputhilohi, inputlolohi,
       inputhihilo, inputlohilo, inputhilolo, inputlololo,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,vrblvl);

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhihihi[i][j][k] = 0.0; jacvallohihi[i][j][k] = 0.0;
            jacvalhilohi[i][j][k] = 0.0; jacvallolohi[i][j][k] = 0.0;
            jacvalhihilo[i][j][k] = 0.0; jacvallohilo[i][j][k] = 0.0;
            jacvalhilolo[i][j][k] = 0.0; jacvallololo[i][j][k] = 0.0;
         }

   dbl8_linearize_evaldiff_output
      (dim,degp1,nvr,idx,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,
       funvalhihihi,funvallohihi,funvalhilohi,funvallolohi,
       funvalhihilo,funvallohilo,funvalhilolo,funvallololo,
       rhshihihi,rhslohihi,rhshilohi,rhslolohi,
       rhshihilo,rhslohilo,rhshilolo,rhslololo,
       jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
       jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,vrblvl);

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         workrhshihihi[i][j] = rhshihihi[i][j];
         workrhslohihi[i][j] = rhslohihi[i][j];
         workrhshilohi[i][j] = rhshilohi[i][j];
         workrhslolohi[i][j] = rhslolohi[i][j];
         workrhshihilo[i][j] = rhshihilo[i][j];
         workrhslohilo[i][j] = rhslohilo[i][j];
         workrhshilolo[i][j] = rhshilolo[i][j];
         workrhslololo[i][j] = rhslololo[i][j];
      }
   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhihihi[i][j] = 0.0; sollohihi[i][j] = 0.0;
         solhilohi[i][j] = 0.0; sollolohi[i][j] = 0.0;
         solhihilo[i][j] = 0.0; sollohilo[i][j] = 0.0;
         solhilolo[i][j] = 0.0; sollololo[i][j] = 0.0;
      }
 
   CPU_dbl8_lusb_solve
      (dim,degp1,
       jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
       jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
       workrhshihihi,workrhslohihi,workrhshilohi,workrhslolohi,
       workrhshihilo,workrhslohilo,workrhshilolo,workrhslololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,
       workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
       workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
       workvechihihi,workveclohihi,workvechilohi,workveclolohi,
       workvechihilo,workveclohilo,workvechilolo,workveclololo,
       ipvt,0); // vrblvl);

   CPU_dbl8_linear_residue
      (dim,degp1,
       jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
       jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
       rhshihihi,rhslohihi,rhshilohi,rhslolohi,
       rhshihilo,rhslohilo,rhshilolo,rhslololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,
       resvechihihi,resveclohihi,resvechilohi,resveclolohi,
       resvechihilo,resveclohilo,resvechilolo,resveclololo,
       resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
       resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);

   if(vrblvl > 0)
      cout << "maximum residual : " << *resmaxhihihi << endl;

   dbl8_update_series
      (dim,degp1,
       inputhihihi,inputlohihi,inputhilohi,inputlolohi,
       inputhihilo,inputlohilo,inputhilolo,inputlololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,vrblvl);
}

void dbl8_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
   double ***outputhihihi_h, double ***outputlohihi_h,
   double ***outputhilohi_h, double ***outputlolohi_h,
   double ***outputhihilo_h, double ***outputlohilo_h,
   double ***outputhilolo_h, double ***outputlololo_h,
   double ***outputhihihi_d, double ***outputlohihi_d,
   double ***outputhilohi_d, double ***outputlolohi_d,
   double ***outputhihilo_d, double ***outputlohilo_d,
   double ***outputhilolo_d, double ***outputlololo_d,
   double **funvalhihihi_h, double **funvallohihi_h,
   double **funvalhilohi_h, double **funvallolohi_h,
   double **funvalhihilo_h, double **funvallohilo_h,
   double **funvalhilolo_h, double **funvallololo_h,
   double **funvalhihihi_d, double **funvallohihi_d,
   double **funvalhilohi_d, double **funvallolohi_d,
   double **funvalhihilo_d, double **funvallohilo_d,
   double **funvalhilolo_d, double **funvallololo_d,
   double ***jacvalhihihi_h, double ***jacvallohihi_h,
   double ***jacvalhilohi_h, double ***jacvallolohi_h,
   double ***jacvalhihilo_h, double ***jacvallohilo_h,
   double ***jacvalhilolo_h, double ***jacvallololo_h,
   double ***jacvalhihihi_d, double ***jacvallohihi_d,
   double ***jacvalhilohi_d, double ***jacvallolohi_d,
   double ***jacvalhihilo_d, double ***jacvallohilo_d,
   double ***jacvalhilolo_d, double ***jacvallololo_d,
   double **rhshihihi_h, double **rhslohihi_h,
   double **rhshilohi_h, double **rhslolohi_h,
   double **rhshihilo_h, double **rhslohilo_h,
   double **rhshilolo_h, double **rhslololo_h,
   double **rhshihihi_d, double **rhslohihi_d,
   double **rhshilohi_d, double **rhslolohi_d,
   double **rhshihilo_d, double **rhslohilo_d,
   double **rhshilolo_d, double **rhslololo_d,
   double **urhshihihi_h, double **urhslohihi_h,
   double **urhshilohi_h, double **urhslolohi_h,
   double **urhshihilo_h, double **urhslohilo_h,
   double **urhshilolo_h, double **urhslololo_h,
   double **urhshihihi_d, double **urhslohihi_d,
   double **urhshilohi_d, double **urhslolohi_d,
   double **urhshihilo_d, double **urhslohilo_d,
   double **urhshilolo_d, double **urhslololo_d,
   double **solhihihi_h, double **sollohihi_h,
   double **solhilohi_h, double **sollolohi_h,
   double **solhihilo_h, double **sollohilo_h,
   double **solhilolo_h, double **sollololo_h,
   double **solhihihi_d, double **sollohihi_d,
   double **solhilohi_d, double **sollolohi_d,
   double **solhihilo_d, double **sollohilo_d,
   double **solhilolo_d, double **sollololo_d,
   double **Qhihihi_h, double **Qlohihi_h,
   double **Qhilohi_h, double **Qlolohi_h,
   double **Qhihilo_h, double **Qlohilo_h,
   double **Qhilolo_h, double **Qlololo_h,
   double **Qhihihi_d, double **Qlohihi_d,
   double **Qhilohi_d, double **Qlolohi_d,
   double **Qhihilo_d, double **Qlohilo_d,
   double **Qhilolo_d, double **Qlololo_d,
   double **Rhihihi_h, double **Rlohihi_h,
   double **Rhilohi_h, double **Rlolohi_h,
   double **Rhihilo_h, double **Rlohilo_h,
   double **Rhilolo_h, double **Rlololo_h,
   double **Rhihihi_d, double **Rlohihi_d,
   double **Rhilohi_d, double **Rlolohi_d,
   double **Rhihilo_d, double **Rlohilo_d,
   double **Rhilolo_d, double **Rlololo_d,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **resvechihihi, double **resveclohihi, 
   double **resvechilohi, double **resveclolohi, 
   double **resvechihilo, double **resveclohilo, 
   double **resvechilolo, double **resveclololo, 
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      // The series coefficients accumulate common factors,
      // initially the coefficients are set to one.
      dbl8_unit_series_vector
         (dim,deg,cffhihihi,cfflohihi,cffhilohi,cfflolohi,
                  cffhihilo,cfflohilo,cffhilolo,cfflololo);

      if(vrblvl > 0)
         cout << "calling CPU_dbl8_evaluate_monomials ..." << endl;

      CPU_dbl8_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      // reset the coefficients
      dbl8_unit_series_vector
         (dim,deg,cffhihihi,cfflohihi,cffhilohi,cfflolohi,
                  cffhihilo,cfflohilo,cffhilolo,cfflololo);

      if(vrblvl > 0)
         cout << "calling GPU_dbl8_evaluate_monomials ..." << endl;

      GPU_dbl8_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
          vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;

      double errsum = 0.0;

      errsum = dbl8_error3sum(dim,dim+1,degp1,
                  outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
                  outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
                  outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
                  outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
                  "output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << "sum of errors : " << errsum << endl;
   }

   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhihihi_h[i][j][k] = 0.0; jacvallohihi_h[i][j][k] = 0.0;
            jacvalhilohi_h[i][j][k] = 0.0; jacvallolohi_h[i][j][k] = 0.0;
            jacvalhihilo_h[i][j][k] = 0.0; jacvallohilo_h[i][j][k] = 0.0;
            jacvalhilolo_h[i][j][k] = 0.0; jacvallololo_h[i][j][k] = 0.0;
            jacvalhihihi_d[i][j][k] = 0.0; jacvallohihi_d[i][j][k] = 0.0;
            jacvalhilohi_d[i][j][k] = 0.0; jacvallolohi_d[i][j][k] = 0.0;
            jacvalhihilo_d[i][j][k] = 0.0; jacvallohilo_d[i][j][k] = 0.0;
            jacvalhilolo_d[i][j][k] = 0.0; jacvallololo_d[i][j][k] = 0.0;
         }

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      dbl8_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
          funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
          rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
          rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,vrblvl);
   }
   if((mode == 1) || (mode == 2))
   {
      dbl8_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
          funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
          funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
          rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
          rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      errsum = dbl8_error2sum(dim,degp1,
                  funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
                  funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
                  funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
                  funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
                  "funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = dbl8_error3sum(degp1,dim,dim,
                  jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
                  jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
                  jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
                  jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
                  "jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = dbl8_error2sum(degp1,dim,
                  rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
                  rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
                  rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
                  rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
                  "rhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhshihihi_h[i][j] = rhshihihi_h[i][j];
         urhslohihi_h[i][j] = rhslohihi_h[i][j];
         urhshilohi_h[i][j] = rhshilohi_h[i][j];
         urhslolohi_h[i][j] = rhslolohi_h[i][j];
         urhshihilo_h[i][j] = rhshihilo_h[i][j];
         urhslohilo_h[i][j] = rhslohilo_h[i][j];
         urhshilolo_h[i][j] = rhshilolo_h[i][j];
         urhslololo_h[i][j] = rhslololo_h[i][j];
         urhshihihi_d[i][j] = rhshihihi_d[i][j];
         urhslohihi_d[i][j] = rhslohihi_d[i][j];
         urhshilohi_d[i][j] = rhshilohi_d[i][j];
         urhslolohi_d[i][j] = rhslolohi_d[i][j];
         urhshihilo_d[i][j] = rhshihilo_d[i][j];
         urhslohilo_d[i][j] = rhslohilo_d[i][j];
         urhshilolo_d[i][j] = rhshilolo_d[i][j];
         urhslololo_d[i][j] = rhslololo_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhihihi_h[i][j] = 0.0; sollohihi_h[i][j] = 0.0;
         solhilohi_h[i][j] = 0.0; sollolohi_h[i][j] = 0.0;
         solhihilo_h[i][j] = 0.0; sollohilo_h[i][j] = 0.0;
         solhilolo_h[i][j] = 0.0; sollololo_h[i][j] = 0.0;
         solhihihi_d[i][j] = 0.0; sollohihi_d[i][j] = 0.0;
         solhilohi_d[i][j] = 0.0; sollolohi_d[i][j] = 0.0;
         solhihilo_d[i][j] = 0.0; sollohilo_d[i][j] = 0.0;
         solhilolo_d[i][j] = 0.0; sollololo_d[i][j] = 0.0;
      }
 
   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling CPU_dbl8_qrbs_solve ..." << endl;
      CPU_dbl8_qrbs_solve
         (dim,degp1,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
          urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
          urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl8_linear_residue ..." << endl;

         CPU_dbl8_linear_residue
            (dim,degp1,
             jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
             jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
             rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
             rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
             solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
             solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
             resvechihihi,resveclohihi,resvechilohi,resveclolohi,
             resvechihilo,resveclohilo,resvechilolo,resveclololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);

         cout << "maximum residual : " << *resmaxhihihi << endl;
      }
      dbl8_update_series
         (dim,degp1,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling GPU_dbl8_bals_solve ..." << endl;

      GPU_dbl8_bals_solve
         (dim,degp1,szt,nbt,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
          urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl8_linear_residue ..." << endl;

         GPU_dbl8_linear_residue
            (dim,degp1,szt,nbt,
             jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
             jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
             rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
             rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
             solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
             solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
             resvechihihi,resveclohihi,resvechilohi,resveclolohi,
             resvechihilo,resveclohilo,resvechilolo,resveclololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);

         cout << "maximum residual : " << *resmaxhihihi << endl;
      }
      dbl8_update_series
         (dim,degp1,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;

      cout << "comparing CPU with GPU matrices Q ... " << endl;
      errsum = dbl8_error2sum(dim,dim,
                  Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
                  Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
                  Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
                  Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = dbl8_error2sum(dim,dim,
                  Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
                  Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
                  Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
                  Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = dbl8_error2sum(degp1,dim,
                  urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
                  urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
                  urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
                  urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
                  "urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = dbl8_error2sum(degp1,dim,
                  solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
                  solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
                  solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
                  solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
                  "sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = dbl8_error2sum(dim,degp1,
                  inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
                  inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
                  inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
                  inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
                  "input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

void cmplx8_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihihi, double **cffrelohihi,
   double **cffrehilohi, double **cffrelolohi,
   double **cffrehihilo, double **cffrelohilo,
   double **cffrehilolo, double **cffrelololo,
   double **cffimhihihi, double **cffimlohihi,
   double **cffimhilohi, double **cffimlolohi,
   double **cffimhihilo, double **cffimlohilo,
   double **cffimhilolo, double **cffimlololo,
   double *accrehihihi, double *accrelohihi,
   double *accrehilohi, double *accrelolohi,
   double *accrehihilo, double *accrelohilo,
   double *accrehilolo, double *accrelololo,
   double *accimhihihi, double *accimlohihi,
   double *accimhilohi, double *accimlolohi,
   double *accimhihilo, double *accimlohilo,
   double *accimhilolo, double *accimlololo,
   double **inputrehihihi_h, double **inputrelohihi_h,
   double **inputrehilohi_h, double **inputrelolohi_h,
   double **inputrehihilo_h, double **inputrelohilo_h,
   double **inputrehilolo_h, double **inputrelololo_h,
   double **inputimhihihi_h, double **inputimlohihi_h,
   double **inputimhilohi_h, double **inputimlolohi_h,
   double **inputimhihilo_h, double **inputimlohilo_h,
   double **inputimhilolo_h, double **inputimlololo_h,
   double **inputrehihihi_d, double **inputrelohihi_d,
   double **inputrehilohi_d, double **inputrelolohi_d,
   double **inputrehihilo_d, double **inputrelohilo_d,
   double **inputrehilolo_d, double **inputrelololo_d,
   double **inputimhihihi_d, double **inputimlohihi_d,
   double **inputimhilohi_d, double **inputimlolohi_d,
   double **inputimhihilo_d, double **inputimlohilo_d,
   double **inputimhilolo_d, double **inputimlololo_d,
   double ***outputrehihihi_h, double ***outputrelohihi_h,
   double ***outputrehilohi_h, double ***outputrelolohi_h,
   double ***outputrehihilo_h, double ***outputrelohilo_h,
   double ***outputrehilolo_h, double ***outputrelololo_h,
   double ***outputimhihihi_h, double ***outputimlohihi_h,
   double ***outputimhilohi_h, double ***outputimlolohi_h,
   double ***outputimhihilo_h, double ***outputimlohilo_h,
   double ***outputimhilolo_h, double ***outputimlololo_h,
   double ***outputrehihihi_d, double ***outputrelohihi_d,
   double ***outputrehilohi_d, double ***outputrelolohi_d,
   double ***outputrehihilo_d, double ***outputrelohilo_d,
   double ***outputrehilolo_d, double ***outputrelololo_d,
   double ***outputimhihihi_d, double ***outputimlohihi_d,
   double ***outputimhilohi_d, double ***outputimlolohi_d,
   double ***outputimhihilo_d, double ***outputimlohilo_d,
   double ***outputimhilolo_d, double ***outputimlololo_d,
   double **funvalrehihihi_h, double **funvalrelohihi_h,
   double **funvalrehilohi_h, double **funvalrelolohi_h,
   double **funvalrehihilo_h, double **funvalrelohilo_h,
   double **funvalrehilolo_h, double **funvalrelololo_h,
   double **funvalimhihihi_h, double **funvalimlohihi_h,
   double **funvalimhilohi_h, double **funvalimlolohi_h,
   double **funvalimhihilo_h, double **funvalimlohilo_h,
   double **funvalimhilolo_h, double **funvalimlololo_h,
   double **funvalrehihihi_d, double **funvalrelohihi_d,
   double **funvalrehilohi_d, double **funvalrelolohi_d,
   double **funvalrehihilo_d, double **funvalrelohilo_d,
   double **funvalrehilolo_d, double **funvalrelololo_d,
   double **funvalimhihihi_d, double **funvalimlohihi_d,
   double **funvalimhilohi_d, double **funvalimlolohi_d,
   double **funvalimhihilo_d, double **funvalimlohilo_d,
   double **funvalimhilolo_d, double **funvalimlololo_d,
   double ***jacvalrehihihi_h, double ***jacvalrelohihi_h,
   double ***jacvalrehilohi_h, double ***jacvalrelolohi_h,
   double ***jacvalrehihilo_h, double ***jacvalrelohilo_h,
   double ***jacvalrehilolo_h, double ***jacvalrelololo_h,
   double ***jacvalimhihihi_h, double ***jacvalimlohihi_h,
   double ***jacvalimhilohi_h, double ***jacvalimlolohi_h,
   double ***jacvalimhihilo_h, double ***jacvalimlohilo_h,
   double ***jacvalimhilolo_h, double ***jacvalimlololo_h,
   double ***jacvalrehihihi_d, double ***jacvalrelohihi_d,
   double ***jacvalrehilohi_d, double ***jacvalrelolohi_d,
   double ***jacvalrehihilo_d, double ***jacvalrelohilo_d,
   double ***jacvalrehilolo_d, double ***jacvalrelololo_d,
   double ***jacvalimhihihi_d, double ***jacvalimlohihi_d,
   double ***jacvalimhilohi_d, double ***jacvalimlolohi_d,
   double ***jacvalimhihilo_d, double ***jacvalimlohilo_d,
   double ***jacvalimhilolo_d, double ***jacvalimlololo_d,
   double **rhsrehihihi_h, double **rhsrelohihi_h,
   double **rhsrehilohi_h, double **rhsrelolohi_h,
   double **rhsrehihilo_h, double **rhsrelohilo_h,
   double **rhsrehilolo_h, double **rhsrelololo_h,
   double **rhsimhihihi_h, double **rhsimlohihi_h,
   double **rhsimhilohi_h, double **rhsimlolohi_h,
   double **rhsimhihilo_h, double **rhsimlohilo_h,
   double **rhsimhilolo_h, double **rhsimlololo_h,
   double **rhsrehihihi_d, double **rhsrelohihi_d, 
   double **rhsrehilohi_d, double **rhsrelolohi_d, 
   double **rhsrehihilo_d, double **rhsrelohilo_d, 
   double **rhsrehilolo_d, double **rhsrelololo_d, 
   double **rhsimhihihi_d, double **rhsimlohihi_d,
   double **rhsimhilohi_d, double **rhsimlolohi_d,
   double **rhsimhihilo_d, double **rhsimlohilo_d,
   double **rhsimhilolo_d, double **rhsimlololo_d,
   double **urhsrehihihi_h, double **urhsrelohihi_h,
   double **urhsrehilohi_h, double **urhsrelolohi_h,
   double **urhsrehihilo_h, double **urhsrelohilo_h,
   double **urhsrehilolo_h, double **urhsrelololo_h,
   double **urhsimhihihi_h, double **urhsimlohihi_h,
   double **urhsimhilohi_h, double **urhsimlolohi_h,
   double **urhsimhihilo_h, double **urhsimlohilo_h,
   double **urhsimhilolo_h, double **urhsimlololo_h,
   double **urhsrehihihi_d, double **urhsrelohihi_d,
   double **urhsrehilohi_d, double **urhsrelolohi_d,
   double **urhsrehihilo_d, double **urhsrelohilo_d,
   double **urhsrehilolo_d, double **urhsrelololo_d,
   double **urhsimhihihi_d, double **urhsimlohihi_d,
   double **urhsimhilohi_d, double **urhsimlolohi_d,
   double **urhsimhihilo_d, double **urhsimlohilo_d,
   double **urhsimhilolo_d, double **urhsimlololo_d,
   double **solrehihihi_h, double **solrelohihi_h,
   double **solrehilohi_h, double **solrelolohi_h,
   double **solrehihilo_h, double **solrelohilo_h,
   double **solrehilolo_h, double **solrelololo_h,
   double **solimhihihi_h, double **solimlohihi_h, 
   double **solimhilohi_h, double **solimlolohi_h, 
   double **solimhihilo_h, double **solimlohilo_h, 
   double **solimhilolo_h, double **solimlololo_h, 
   double **solrehihihi_d, double **solrelohihi_d, 
   double **solrehilohi_d, double **solrelolohi_d, 
   double **solrehihilo_d, double **solrelohilo_d, 
   double **solrehilolo_d, double **solrelololo_d, 
   double **solimhihihi_d, double **solimlohihi_d,
   double **solimhilohi_d, double **solimlolohi_d,
   double **solimhihilo_d, double **solimlohilo_d,
   double **solimhilolo_d, double **solimlololo_d,
   double **Qrehihihi_h, double **Qrelohihi_h,
   double **Qrehilohi_h, double **Qrelolohi_h,
   double **Qrehihilo_h, double **Qrelohilo_h,
   double **Qrehilolo_h, double **Qrelololo_h,
   double **Qimhihihi_h, double **Qimlohihi_h,
   double **Qimhilohi_h, double **Qimlolohi_h,
   double **Qimhihilo_h, double **Qimlohilo_h,
   double **Qimhilolo_h, double **Qimlololo_h,
   double **Qrehihihi_d, double **Qrelohihi_d,
   double **Qrehilohi_d, double **Qrelolohi_d,
   double **Qrehihilo_d, double **Qrelohilo_d,
   double **Qrehilolo_d, double **Qrelololo_d,
   double **Qimhihihi_d, double **Qimlohihi_d, 
   double **Qimhilohi_d, double **Qimlolohi_d, 
   double **Qimhihilo_d, double **Qimlohilo_d, 
   double **Qimhilolo_d, double **Qimlololo_d, 
   double **Rrehihihi_h, double **Rrelohihi_h,
   double **Rrehilohi_h, double **Rrelolohi_h,
   double **Rrehihilo_h, double **Rrelohilo_h,
   double **Rrehilolo_h, double **Rrelololo_h,
   double **Rimhihihi_h, double **Rimlohihi_h, 
   double **Rimhilohi_h, double **Rimlolohi_h, 
   double **Rimhihilo_h, double **Rimlohilo_h, 
   double **Rimhilolo_h, double **Rimlololo_h, 
   double **Rrehihihi_d, double **Rrelohihi_d,
   double **Rrehilohi_d, double **Rrelolohi_d,
   double **Rrehihilo_d, double **Rrelohilo_d,
   double **Rrehilolo_d, double **Rrelololo_d,
   double **Rimhihihi_d, double **Rimlohihi_d,
   double **Rimhilohi_d, double **Rimlolohi_d,
   double **Rimhihilo_d, double **Rimlohilo_d,
   double **Rimhilolo_d, double **Rimlololo_d,
   double **workmatrehihihi, double **workmatrelohihi,
   double **workmatrehilohi, double **workmatrelolohi,
   double **workmatrehihilo, double **workmatrelohilo,
   double **workmatrehilolo, double **workmatrelololo,
   double **workmatimhihihi, double **workmatimlohihi,
   double **workmatimhilohi, double **workmatimlolohi,
   double **workmatimhihilo, double **workmatimlohilo,
   double **workmatimhilolo, double **workmatimlololo,
   double *workvecrehihihi, double *workvecrelohihi,
   double *workvecrehilohi, double *workvecrelolohi,
   double *workvecrehihilo, double *workvecrelohilo,
   double *workvecrehilolo, double *workvecrelololo,
   double *workvecimhihihi, double *workvecimlohihi,
   double *workvecimhilohi, double *workvecimlolohi,
   double *workvecimhihilo, double *workvecimlohilo,
   double *workvecimhilolo, double *workvecimlololo,
   double **resvecrehihihi, double **resvecrelohihi,
   double **resvecrehilohi, double **resvecrelolohi,
   double **resvecrehihilo, double **resvecrelohilo,
   double **resvecrehilolo, double **resvecrelololo,
   double **resvecimhihihi, double **resvecimlohihi,
   double **resvecimhilohi, double **resvecimlolohi,
   double **resvecimhihilo, double **resvecimlohilo,
   double **resvecimhilolo, double **resvecimlololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      // The series coefficients accumulate common factors,
      // initially the coefficients are set to one.
      cmplx8_unit_series_vector
         (dim,deg,cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
                  cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
                  cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
                  cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);

      if(vrblvl > 0)
         cout << "calling CPU_cmplx8_evaluate_monomials ..." << endl;

      CPU_cmplx8_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
          accrehihihi,accrelohihi,accrehilohi,accrelolohi,
          accrehihilo,accrelohilo,accrehilolo,accrelololo,
          accimhihihi,accimlohihi,accimhilohi,accimlolohi,
          accimhihilo,accimlohilo,accimhilolo,accimlololo,
          inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
          inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
          inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
          inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h,
          outputrehihihi_h,outputrelohihi_h,outputrehilohi_h,outputrelolohi_h,
          outputrehihilo_h,outputrelohilo_h,outputrehilolo_h,outputrelololo_h,
          outputimhihihi_h,outputimlohihi_h,outputimhilohi_h,outputimlolohi_h,
          outputimhihilo_h,outputimlohilo_h,outputimhilolo_h,outputimlololo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      // reset the coefficients
      cmplx8_unit_series_vector
         (dim,deg,cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
                  cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
                  cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
                  cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);

      if(vrblvl > 0)
         cout << "calling GPU_cmplx8_evaluate_monomials ..." << endl;

      GPU_cmplx8_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
          accrehihihi,accrelohihi,accrehilohi,accrelolohi,
          accrehihilo,accrelohilo,accrehilolo,accrelololo,
          accimhihihi,accimlohihi,accimhilohi,accimlolohi,
          accimhihilo,accimlohilo,accimhilolo,accimlololo,
          inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
          inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
          inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
          inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
          outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
          outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
          outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
          outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d,
          vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;

      cout << "comparing CPU with GPU evaluations ... " << endl;

      errsum = cmplx8_error3sum(dim,dim+1,degp1,
                  outputrehihihi_h,outputrelohihi_h,
                  outputrehilohi_h,outputrelolohi_h,
                  outputrehihilo_h,outputrelohilo_h,
                  outputrehilolo_h,outputrelololo_h,
                  outputimhihihi_h,outputimlohihi_h,
                  outputimhilohi_h,outputimlolohi_h,
                  outputimhihilo_h,outputimlohilo_h,
                  outputimhilolo_h,outputimlololo_h,
                  outputrehihihi_d,outputrelohihi_d,
                  outputrehilohi_d,outputrelolohi_d,
                  outputrehihilo_d,outputrelohilo_d,
                  outputrehilolo_d,outputrelololo_d,
                  outputimhihihi_d,outputimlohihi_d,
                  outputimhilohi_d,outputimlolohi_d,
                  outputimhihilo_d,outputimlohilo_d,
                  outputimhilolo_d,outputimlololo_d,"output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalrehihihi_h[i][j][k] = 0.0; jacvalrelohihi_h[i][j][k] = 0.0;
            jacvalrehilohi_h[i][j][k] = 0.0; jacvalrelolohi_h[i][j][k] = 0.0;
            jacvalimhihilo_h[i][j][k] = 0.0; jacvalimlohilo_h[i][j][k] = 0.0;
            jacvalimhilolo_h[i][j][k] = 0.0; jacvalimlololo_h[i][j][k] = 0.0;
            jacvalrehihihi_d[i][j][k] = 0.0; jacvalrelohihi_d[i][j][k] = 0.0;
            jacvalrehilohi_d[i][j][k] = 0.0; jacvalrelolohi_d[i][j][k] = 0.0;
            jacvalimhihilo_d[i][j][k] = 0.0; jacvalimlohilo_d[i][j][k] = 0.0;
            jacvalimhilolo_d[i][j][k] = 0.0; jacvalimlololo_d[i][j][k] = 0.0;
         }

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      cmplx8_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputrehihihi_h,outputrelohihi_h,outputrehilohi_h,outputrelolohi_h,
          outputrehihilo_h,outputrelohilo_h,outputrehilolo_h,outputrelololo_h,
          outputimhihihi_h,outputimlohihi_h,outputimhilohi_h,outputimlolohi_h,
          outputimhihilo_h,outputimlohilo_h,outputimhilolo_h,outputimlololo_h,
          funvalrehihihi_h,funvalrelohihi_h,funvalrehilohi_h,funvalrelolohi_h,
          funvalrehihilo_h,funvalrelohilo_h,funvalrehilolo_h,funvalrelololo_h,
          funvalimhihihi_h,funvalimlohihi_h,funvalimhilohi_h,funvalimlolohi_h,
          funvalimhihilo_h,funvalimlohilo_h,funvalimhilolo_h,funvalimlololo_h,
          rhsrehihihi_h,rhsrelohihi_h,rhsrehilohi_h,rhsrelolohi_h,
          rhsrehihilo_h,rhsrelohilo_h,rhsrehilolo_h,rhsrelololo_h,
          rhsimhihihi_h,rhsimlohihi_h,rhsimhilohi_h,rhsimlolohi_h,
          rhsimhihilo_h,rhsimlohilo_h,rhsimhilolo_h,rhsimlololo_h,
          jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
          jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
          jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
          jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx8_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
          outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
          outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
          outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d,
          funvalrehihihi_d,funvalrelohihi_d,funvalrehilohi_d,funvalrelolohi_d,
          funvalrehihilo_d,funvalrelohilo_d,funvalrehilolo_d,funvalrelololo_d,
          funvalimhihihi_d,funvalimlohihi_d,funvalimhilohi_d,funvalimlolohi_d,
          funvalimhihilo_d,funvalimlohilo_d,funvalimhilolo_d,funvalimlololo_d,
          rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
          rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
          rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
          rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
          jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
          jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
          jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
          jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
          vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU function values ... " << endl;
      errsum = cmplx8_error2sum(dim,degp1,
                  funvalrehihihi_h,funvalrelohihi_h,
                  funvalrehilohi_h,funvalrelolohi_h,
                  funvalrehihilo_h,funvalrelohilo_h,
                  funvalrehilolo_h,funvalrelololo_h,
                  funvalimhihihi_h,funvalimlohihi_h,
                  funvalimhilohi_h,funvalimlolohi_h,
                  funvalimhihilo_h,funvalimlohilo_h,
                  funvalimhilolo_h,funvalimlololo_h,
                  funvalrehihihi_d,funvalrelohihi_d,
                  funvalrehilohi_d,funvalrelolohi_d,
                  funvalrehihilo_d,funvalrelohilo_d,
                  funvalrehilolo_d,funvalrelololo_d,
                  funvalimhihihi_d,funvalimlohihi_d,
                  funvalimhilohi_d,funvalimlolohi_d,
                  funvalimhihilo_d,funvalimlohilo_d,
                  funvalimhilolo_d,funvalimlololo_d,"funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = cmplx8_error3sum(degp1,dim,dim,
                  jacvalrehihihi_h,jacvalrelohihi_h,
                  jacvalrehilohi_h,jacvalrelolohi_h,
                  jacvalrehihilo_h,jacvalrelohilo_h,
                  jacvalrehilolo_h,jacvalrelololo_h,
                  jacvalimhihihi_h,jacvalimlohihi_h,
                  jacvalimhilohi_h,jacvalimlolohi_h,
                  jacvalimhihilo_h,jacvalimlohilo_h,
                  jacvalimhilolo_h,jacvalimlololo_h,
                  jacvalrehihihi_d,jacvalrelohihi_d,
                  jacvalrehilohi_d,jacvalrelolohi_d,
                  jacvalrehihilo_d,jacvalrelohilo_d,
                  jacvalrehilolo_d,jacvalrelololo_d,
                  jacvalimhihihi_d,jacvalimlohihi_d,
                  jacvalimhilohi_d,jacvalimlolohi_d,
                  jacvalimhihilo_d,jacvalimlohilo_d,
                  jacvalimhilolo_d,jacvalimlololo_d,"jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = cmplx8_error2sum(degp1,dim,
                  rhsrehihihi_h,rhsrelohihi_h,rhsrehilohi_h,rhsrelolohi_h,
                  rhsrehihilo_h,rhsrelohilo_h,rhsrehilolo_h,rhsrelololo_h,
                  rhsimhihihi_h,rhsimlohihi_h,rhsimhilohi_h,rhsimlolohi_h,
                  rhsimhihilo_h,rhsimlohilo_h,rhsimhilolo_h,rhsimlololo_h,
                  rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
                  rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
                  rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
                  rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
                  "rhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhsrehihihi_h[i][j] = rhsrehihihi_h[i][j];
         urhsrehilohi_h[i][j] = rhsrehilohi_h[i][j];
         urhsrelohihi_h[i][j] = rhsrelohihi_h[i][j];
         urhsrelolohi_h[i][j] = rhsrelolohi_h[i][j];
         urhsrehihilo_h[i][j] = rhsrehihilo_h[i][j];
         urhsrehilolo_h[i][j] = rhsrehilolo_h[i][j];
         urhsrelohilo_h[i][j] = rhsrelohilo_h[i][j];
         urhsrelololo_h[i][j] = rhsrelololo_h[i][j];
         urhsimhihihi_h[i][j] = rhsimhihihi_h[i][j];
         urhsimhilohi_h[i][j] = rhsimhilohi_h[i][j];
         urhsimlohihi_h[i][j] = rhsimlohihi_h[i][j];
         urhsimlolohi_h[i][j] = rhsimlolohi_h[i][j];
         urhsimhihilo_h[i][j] = rhsimhihilo_h[i][j];
         urhsimhilolo_h[i][j] = rhsimhilolo_h[i][j];
         urhsimlohilo_h[i][j] = rhsimlohilo_h[i][j];
         urhsimlololo_h[i][j] = rhsimlololo_h[i][j];
         urhsrehihihi_d[i][j] = rhsrehihihi_d[i][j];
         urhsrehilohi_d[i][j] = rhsrehilohi_d[i][j];
         urhsrelohihi_d[i][j] = rhsrelohihi_d[i][j];
         urhsrelolohi_d[i][j] = rhsrelolohi_d[i][j];
         urhsrehihilo_d[i][j] = rhsrehihilo_d[i][j];
         urhsrehilolo_d[i][j] = rhsrehilolo_d[i][j];
         urhsrelohilo_d[i][j] = rhsrelohilo_d[i][j];
         urhsrelololo_d[i][j] = rhsrelololo_d[i][j];
         urhsimhihihi_d[i][j] = rhsimhihihi_d[i][j];
         urhsimhilohi_d[i][j] = rhsimhilohi_d[i][j];
         urhsimlohihi_d[i][j] = rhsimlohihi_d[i][j];
         urhsimlolohi_d[i][j] = rhsimlolohi_d[i][j];
         urhsimhihilo_d[i][j] = rhsimhihilo_d[i][j];
         urhsimhilolo_d[i][j] = rhsimhilolo_d[i][j];
         urhsimlohilo_d[i][j] = rhsimlohilo_d[i][j];
         urhsimlololo_d[i][j] = rhsimlololo_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solrehihihi_h[i][j] = 0.0; solimhihihi_h[i][j] = 0.0;
         solrelohihi_h[i][j] = 0.0; solimlohihi_h[i][j] = 0.0;
         solrehilohi_h[i][j] = 0.0; solimhilohi_h[i][j] = 0.0;
         solrelolohi_h[i][j] = 0.0; solimlolohi_h[i][j] = 0.0;
         solrehihilo_h[i][j] = 0.0; solimhihilo_h[i][j] = 0.0;
         solrelohilo_h[i][j] = 0.0; solimlohilo_h[i][j] = 0.0;
         solrehilolo_h[i][j] = 0.0; solimhilolo_h[i][j] = 0.0;
         solrelololo_h[i][j] = 0.0; solimlololo_h[i][j] = 0.0;
         solrehihihi_d[i][j] = 0.0; solimhihihi_d[i][j] = 0.0;
         solrelohihi_d[i][j] = 0.0; solimlohihi_d[i][j] = 0.0;
         solrehilohi_d[i][j] = 0.0; solimhilohi_d[i][j] = 0.0;
         solrelolohi_d[i][j] = 0.0; solimlolohi_d[i][j] = 0.0;
         solrehihilo_d[i][j] = 0.0; solimhihilo_d[i][j] = 0.0;
         solrelohilo_d[i][j] = 0.0; solimlohilo_d[i][j] = 0.0;
         solrehilolo_d[i][j] = 0.0; solimhilolo_d[i][j] = 0.0;
         solrelololo_d[i][j] = 0.0; solimlololo_d[i][j] = 0.0;
      }

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_cmplx8_qrbs_solve ..." << endl;

      CPU_cmplx8_qrbs_solve
         (dim,degp1,
          jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
          jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
          jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
          jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,
          urhsrehihihi_h,urhsrelohihi_h,urhsrehilohi_h,urhsrelolohi_h,
          urhsrehihilo_h,urhsrelohilo_h,urhsrehilolo_h,urhsrelololo_h,
          urhsimhihihi_h,urhsimlohihi_h,urhsimhilohi_h,urhsimlolohi_h,
          urhsimhihilo_h,urhsimlohilo_h,urhsimhilolo_h,urhsimlololo_h,
          solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
          solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
          solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
          solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,
          workmatrehihihi,workmatrelohihi,workmatrehilohi,workmatrelolohi,
          workmatrehihilo,workmatrelohilo,workmatrehilolo,workmatrelololo,
          workmatimhihihi,workmatimlohihi,workmatimhilohi,workmatimlolohi,
          workmatimhihilo,workmatimlohilo,workmatimhilolo,workmatimlololo,
          Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
          Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
          Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
          Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
          Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
          Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
          Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
          Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,
          workvecrehihihi,workvecrelohihi,workvecrehilohi,workvecrelolohi,
          workvecrehihilo,workvecrelohilo,workvecrehilolo,workvecrelololo,
          workvecimhihihi,workvecimlohihi,workvecimhilohi,workvecimlolohi,
          workvecimhihilo,workvecimlohilo,workvecimhilolo,workvecimlololo,
          vrblvl);
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_cmplx8_linear_residue ..." << endl;

         CPU_cmplx8_linear_residue
            (dim,degp1,
             jacvalrehihihi_h,jacvalrelohihi_h,
             jacvalrehilohi_h,jacvalrelolohi_h,
             jacvalrehihilo_h,jacvalrelohilo_h,
             jacvalrehilolo_h,jacvalrelololo_h,
             jacvalimhihihi_h,jacvalimlohihi_h,
             jacvalimhilohi_h,jacvalimlolohi_h,
             jacvalimhihilo_h,jacvalimlohilo_h,
             jacvalimhilolo_h,jacvalimlololo_h,
             rhsrehihihi_h,rhsrelohihi_h,rhsrehilohi_h,rhsrelolohi_h,
             rhsrehihilo_h,rhsrelohilo_h,rhsrehilolo_h,rhsrelololo_h,
             rhsimhihihi_h,rhsimlohihi_h,rhsimhilohi_h,rhsimlolohi_h,
             rhsimhihilo_h,rhsimlohilo_h,rhsimhilolo_h,rhsimlololo_h,
             solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
             solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
             solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
             solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,
             resvecrehihihi,resvecrelohihi,resvecrehilohi,resvecrelolohi,
             resvecrehihilo,resvecrelohilo,resvecrehilolo,resvecrelololo,
             resvecimhihihi,resvecimlohihi,resvecimhilohi,resvecimlolohi,
             resvecimhihilo,resvecimlohilo,resvecimhilolo,resvecimlololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,
             vrblvl);
         cout << "maximum residual : " << *resmaxhihihi << endl;
      }
      cmplx8_update_series
         (dim,degp1,
          inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
          inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
          inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
          inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h,
          solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
          solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
          solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
          solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_cmplx8_bals_solve ..." << endl;

      GPU_cmplx8_bals_solve
         (dim,degp1,szt,nbt,
          jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
          jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
          jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
          jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
          Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
          Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
          Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
          Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,
          Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
          Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
          Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
          Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,
          urhsrehihihi_d,urhsrelohihi_d,urhsrehilohi_d,urhsrelolohi_d,
          urhsrehihilo_d,urhsrelohilo_d,urhsrehilolo_d,urhsrelololo_d,
          urhsimhihihi_d,urhsimlohihi_d,urhsimhilohi_d,urhsimlolohi_d,
          urhsimhihilo_d,urhsimlohilo_d,urhsimhilolo_d,urhsimlololo_d,
          solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
          solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
          solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
          solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling GPU_cmplx8_linear_residue ..." << endl;

         GPU_cmplx8_linear_residue
            (dim,degp1,szt,nbt,
             jacvalrehihihi_d,jacvalrelohihi_d,
             jacvalrehilohi_d,jacvalrelolohi_d,
             jacvalrehihilo_d,jacvalrelohilo_d,
             jacvalrehilolo_d,jacvalrelololo_d,
             jacvalimhihihi_d,jacvalimlohihi_d,
             jacvalimhilohi_d,jacvalimlolohi_d,
             jacvalimhihilo_d,jacvalimlohilo_d,
             jacvalimhilolo_d,jacvalimlololo_d,
             rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
             rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
             rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
             rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
             solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
             solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
             solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
             solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,
             resvecrehihihi,resvecrelohihi,resvecrehilohi,resvecrelolohi,
             resvecrehihilo,resvecrelohilo,resvecrehilolo,resvecrelololo,
             resvecimhihihi,resvecimlohihi,resvecimhilohi,resvecimlolohi,
             resvecimhihilo,resvecimlohilo,resvecimhilolo,resvecimlololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);
         cout << "maximum residual : " << *resmaxhihihi << endl;
      }
      cmplx8_update_series
         (dim,degp1,
          inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
          inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
          inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
          inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
          solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
          solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
          solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
          solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;

      cout << "comparing CPU with GPU matrices Q ... " << endl;
      errsum = cmplx8_error2sum(dim,dim,
                  Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
                  Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
                  Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
                  Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
                  Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
                  Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
                  Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
                  Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = cmplx8_error2sum(dim,dim,
                  Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
                  Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
                  Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
                  Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,
                  Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
                  Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
                  Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
                  Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = cmplx8_error2sum(degp1,dim,
                  urhsrehihihi_h,urhsrelohihi_h,urhsrehilohi_h,urhsrelolohi_h,
                  urhsrehihilo_h,urhsrelohilo_h,urhsrehilolo_h,urhsrelololo_h,
                  urhsimhihihi_h,urhsimlohihi_h,urhsimhilohi_h,urhsimlolohi_h,
                  urhsimhihilo_h,urhsimlohilo_h,urhsimhilolo_h,urhsimlololo_h,
                  urhsrehihihi_d,urhsrelohihi_d,urhsrehilohi_d,urhsrelolohi_d,
                  urhsrehihilo_d,urhsrelohilo_d,urhsrehilolo_d,urhsrelololo_d,
                  urhsimhihihi_d,urhsimlohihi_d,urhsimhilohi_d,urhsimlolohi_d,
                  urhsimhihilo_d,urhsimlohilo_d,urhsimhilolo_d,urhsimlololo_d,
                  "urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = cmplx8_error2sum(degp1,dim,
                  solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
                  solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
                  solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
                  solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,
                  solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
                  solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
                  solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
                  solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,
                  "sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = cmplx8_error2sum(dim,degp1,
                  inputrehihihi_h,inputrelohihi_h,
                  inputrehilohi_h,inputrelolohi_h,
                  inputrehihilo_h,inputrelohilo_h,
                  inputrehilolo_h,inputrelololo_h,
                  inputimhihihi_h,inputimlohihi_h,
                  inputimhilohi_h,inputimlolohi_h,
                  inputimhihilo_h,inputimlohilo_h,
                  inputimhilolo_h,inputimlololo_h,
                  inputrehihihi_d,inputrelohihi_d,
                  inputrehilohi_d,inputrelolohi_d,
                  inputrehihilo_d,inputrelohilo_d,
                  inputrehilolo_d,inputrelololo_d,
                  inputimhihihi_d,inputimlohihi_d,
                  inputimhilohi_d,inputimlolohi_d,
                  inputimhihilo_d,inputimlohilo_d,
                  inputimhilolo_d,inputimlololo_d,"input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl8_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhihihi_h = new double*[dim];
   double **inputlohihi_h = new double*[dim];
   double **inputhilohi_h = new double*[dim];
   double **inputlolohi_h = new double*[dim];
   double **inputhihilo_h = new double*[dim];
   double **inputlohilo_h = new double*[dim];
   double **inputhilolo_h = new double*[dim];
   double **inputlololo_h = new double*[dim];
   double **inputhihihi_d = new double*[dim];
   double **inputlohihi_d = new double*[dim];
   double **inputhilohi_d = new double*[dim];
   double **inputlolohi_d = new double*[dim];
   double **inputhihilo_d = new double*[dim];
   double **inputlohilo_d = new double*[dim];
   double **inputhilolo_d = new double*[dim];
   double **inputlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihihi_h[i] = new double[degp1];
      inputlohihi_h[i] = new double[degp1];
      inputhilohi_h[i] = new double[degp1];
      inputlolohi_h[i] = new double[degp1];
      inputhihilo_h[i] = new double[degp1];
      inputlohilo_h[i] = new double[degp1];
      inputhilolo_h[i] = new double[degp1];
      inputlololo_h[i] = new double[degp1];
      inputhihihi_d[i] = new double[degp1];
      inputlohihi_d[i] = new double[degp1];
      inputhilohi_d[i] = new double[degp1];
      inputlolohi_d[i] = new double[degp1];
      inputhihilo_d[i] = new double[degp1];
      inputlohilo_d[i] = new double[degp1];
      inputhilolo_d[i] = new double[degp1];
      inputlololo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acchihihi = new double[degp1]; // accumulated power series
   double *acclohihi = new double[degp1];
   double *acchilohi = new double[degp1];
   double *acclolohi = new double[degp1];
   double *acchihilo = new double[degp1];
   double *acclohilo = new double[degp1];
   double *acchilolo = new double[degp1];
   double *acclololo = new double[degp1];
   double **cffhihihi = new double*[dim]; // the coefficients of monomials
   double **cfflohihi = new double*[dim];
   double **cffhilohi = new double*[dim];
   double **cfflolohi = new double*[dim];
   double **cffhihilo = new double*[dim];
   double **cfflohilo = new double*[dim];
   double **cffhilolo = new double*[dim];
   double **cfflololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double[degp1];
      cfflohihi[i] = new double[degp1];
      cffhilohi[i] = new double[degp1];
      cfflolohi[i] = new double[degp1];
      cffhihilo[i] = new double[degp1];
      cfflohilo[i] = new double[degp1];
      cffhilolo[i] = new double[degp1];
      cfflololo[i] = new double[degp1];
   }
   double ***outputhihihi_h = new double**[dim];
   double ***outputlohihi_h = new double**[dim];
   double ***outputhilohi_h = new double**[dim];
   double ***outputlolohi_h = new double**[dim];
   double ***outputhihilo_h = new double**[dim];
   double ***outputlohilo_h = new double**[dim];
   double ***outputhilolo_h = new double**[dim];
   double ***outputlololo_h = new double**[dim];
   double ***outputhihihi_d = new double**[dim];
   double ***outputlohihi_d = new double**[dim];
   double ***outputhilohi_d = new double**[dim];
   double ***outputlolohi_d = new double**[dim];
   double ***outputhihilo_d = new double**[dim];
   double ***outputlohilo_d = new double**[dim];
   double ***outputhilolo_d = new double**[dim];
   double ***outputlololo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihihi_h[i] = new double*[dim+1];
      outputlohihi_h[i] = new double*[dim+1];
      outputhilohi_h[i] = new double*[dim+1];
      outputlolohi_h[i] = new double*[dim+1];
      outputhihilo_h[i] = new double*[dim+1];
      outputlohilo_h[i] = new double*[dim+1];
      outputhilolo_h[i] = new double*[dim+1];
      outputlololo_h[i] = new double*[dim+1];
      outputhihihi_d[i] = new double*[dim+1];
      outputlohihi_d[i] = new double*[dim+1];
      outputhilohi_d[i] = new double*[dim+1];
      outputlolohi_d[i] = new double*[dim+1];
      outputhihilo_d[i] = new double*[dim+1];
      outputlohilo_d[i] = new double*[dim+1];
      outputhilolo_d[i] = new double*[dim+1];
      outputlololo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihihi_h[i][j] = new double[degp1];
         outputlohihi_h[i][j] = new double[degp1];
         outputhilohi_h[i][j] = new double[degp1];
         outputlolohi_h[i][j] = new double[degp1];
         outputhihilo_h[i][j] = new double[degp1];
         outputlohilo_h[i][j] = new double[degp1];
         outputhilolo_h[i][j] = new double[degp1];
         outputlololo_h[i][j] = new double[degp1];
         outputhihihi_d[i][j] = new double[degp1];
         outputlohihi_d[i][j] = new double[degp1];
         outputhilohi_d[i][j] = new double[degp1];
         outputlolohi_d[i][j] = new double[degp1];
         outputhihilo_d[i][j] = new double[degp1];
         outputlohilo_d[i][j] = new double[degp1];
         outputhilolo_d[i][j] = new double[degp1];
         outputlololo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihihi_h = new double*[dim];
   double **funvallohihi_h = new double*[dim];
   double **funvalhilohi_h = new double*[dim];
   double **funvallolohi_h = new double*[dim];
   double **funvalhihilo_h = new double*[dim];
   double **funvallohilo_h = new double*[dim];
   double **funvalhilolo_h = new double*[dim];
   double **funvallololo_h = new double*[dim];
   double **funvalhihihi_d = new double*[dim];
   double **funvallohihi_d = new double*[dim];
   double **funvalhilohi_d = new double*[dim];
   double **funvallolohi_d = new double*[dim];
   double **funvalhihilo_d = new double*[dim];
   double **funvallohilo_d = new double*[dim];
   double **funvalhilolo_d = new double*[dim];
   double **funvallololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhihihi_h[i] = new double[degp1];
      funvallohihi_h[i] = new double[degp1];
      funvalhilohi_h[i] = new double[degp1];
      funvallolohi_h[i] = new double[degp1];
      funvalhihilo_h[i] = new double[degp1];
      funvallohilo_h[i] = new double[degp1];
      funvalhilolo_h[i] = new double[degp1];
      funvallololo_h[i] = new double[degp1];
      funvalhihihi_d[i] = new double[degp1];
      funvallohihi_d[i] = new double[degp1];
      funvalhilohi_d[i] = new double[degp1];
      funvallolohi_d[i] = new double[degp1];
      funvalhihilo_d[i] = new double[degp1];
      funvallohilo_d[i] = new double[degp1];
      funvalhilolo_d[i] = new double[degp1];
      funvallololo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihihi_h = new double**[degp1];
   double ***jacvallohihi_h = new double**[degp1];
   double ***jacvalhilohi_h = new double**[degp1];
   double ***jacvallolohi_h = new double**[degp1];
   double ***jacvalhihilo_h = new double**[degp1];
   double ***jacvallohilo_h = new double**[degp1];
   double ***jacvalhilolo_h = new double**[degp1];
   double ***jacvallololo_h = new double**[degp1];
   double ***jacvalhihihi_d = new double**[degp1];
   double ***jacvallohihi_d = new double**[degp1];
   double ***jacvalhilohi_d = new double**[degp1];
   double ***jacvallolohi_d = new double**[degp1];
   double ***jacvalhihilo_d = new double**[degp1];
   double ***jacvallohilo_d = new double**[degp1];
   double ***jacvalhilolo_d = new double**[degp1];
   double ***jacvallololo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhihihi_h[i] = new double*[dim];
      jacvallohihi_h[i] = new double*[dim];
      jacvalhilohi_h[i] = new double*[dim];
      jacvallolohi_h[i] = new double*[dim];
      jacvalhihilo_h[i] = new double*[dim];
      jacvallohilo_h[i] = new double*[dim];
      jacvalhilolo_h[i] = new double*[dim];
      jacvallololo_h[i] = new double*[dim];
      jacvalhihihi_d[i] = new double*[dim];
      jacvallohihi_d[i] = new double*[dim];
      jacvalhilohi_d[i] = new double*[dim];
      jacvallolohi_d[i] = new double*[dim];
      jacvalhihilo_d[i] = new double*[dim];
      jacvallohilo_d[i] = new double*[dim];
      jacvalhilolo_d[i] = new double*[dim];
      jacvallololo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhihihi_h[i][j] = new double[dim];
         jacvallohihi_h[i][j] = new double[dim];
         jacvalhilohi_h[i][j] = new double[dim];
         jacvallolohi_h[i][j] = new double[dim];
         jacvalhihilo_h[i][j] = new double[dim];
         jacvallohilo_h[i][j] = new double[dim];
         jacvalhilolo_h[i][j] = new double[dim];
         jacvallololo_h[i][j] = new double[dim];
         jacvalhihihi_d[i][j] = new double[dim];
         jacvallohihi_d[i][j] = new double[dim];
         jacvalhilohi_d[i][j] = new double[dim];
         jacvallolohi_d[i][j] = new double[dim];
         jacvalhihilo_d[i][j] = new double[dim];
         jacvallohilo_d[i][j] = new double[dim];
         jacvalhilolo_d[i][j] = new double[dim];
         jacvallololo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihihi_h = new double*[degp1];
   double **sollohihi_h = new double*[degp1];
   double **solhilohi_h = new double*[degp1];
   double **sollolohi_h = new double*[degp1];
   double **solhihilo_h = new double*[degp1];
   double **sollohilo_h = new double*[degp1];
   double **solhilolo_h = new double*[degp1];
   double **sollololo_h = new double*[degp1];
   double **solhihihi_d = new double*[degp1];
   double **sollohihi_d = new double*[degp1];
   double **solhilohi_d = new double*[degp1];
   double **sollolohi_d = new double*[degp1];
   double **solhihilo_d = new double*[degp1];
   double **sollohilo_d = new double*[degp1];
   double **solhilolo_d = new double*[degp1];
   double **sollololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihihi_h[i] = new double[dim];
      sollohihi_h[i] = new double[dim];
      solhilohi_h[i] = new double[dim];
      sollolohi_h[i] = new double[dim];
      solhihilo_h[i] = new double[dim];
      sollohilo_h[i] = new double[dim];
      solhilolo_h[i] = new double[dim];
      sollololo_h[i] = new double[dim];
      solhihihi_d[i] = new double[dim];
      sollohihi_d[i] = new double[dim];
      solhilohi_d[i] = new double[dim];
      sollolohi_d[i] = new double[dim];
      solhihilo_d[i] = new double[dim];
      sollohilo_d[i] = new double[dim];
      solhilolo_d[i] = new double[dim];
      sollololo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihihi_h = new double*[degp1];
   double **rhslohihi_h = new double*[degp1];
   double **rhshilohi_h = new double*[degp1];
   double **rhslolohi_h = new double*[degp1];
   double **rhshihilo_h = new double*[degp1];
   double **rhslohilo_h = new double*[degp1];
   double **rhshilolo_h = new double*[degp1];
   double **rhslololo_h = new double*[degp1];
   double **rhshihihi_d = new double*[degp1];
   double **rhslohihi_d = new double*[degp1];
   double **rhshilohi_d = new double*[degp1];
   double **rhslolohi_d = new double*[degp1];
   double **rhshihilo_d = new double*[degp1];
   double **rhslohilo_d = new double*[degp1];
   double **rhshilolo_d = new double*[degp1];
   double **rhslololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihihi_h[i] = new double[dim];
      rhslohihi_h[i] = new double[dim];
      rhshilohi_h[i] = new double[dim];
      rhslolohi_h[i] = new double[dim];
      rhshihilo_h[i] = new double[dim];
      rhslohilo_h[i] = new double[dim];
      rhshilolo_h[i] = new double[dim];
      rhslololo_h[i] = new double[dim];
      rhshihihi_d[i] = new double[dim];
      rhslohihi_d[i] = new double[dim];
      rhshilohi_d[i] = new double[dim];
      rhslolohi_d[i] = new double[dim];
      rhshihilo_d[i] = new double[dim];
      rhslohilo_d[i] = new double[dim];
      rhshilolo_d[i] = new double[dim];
      rhslololo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihihi_h = new double*[degp1];
   double **urhslohihi_h = new double*[degp1];
   double **urhshilohi_h = new double*[degp1];
   double **urhslolohi_h = new double*[degp1];
   double **urhshihilo_h = new double*[degp1];
   double **urhslohilo_h = new double*[degp1];
   double **urhshilolo_h = new double*[degp1];
   double **urhslololo_h = new double*[degp1];
   double **urhshihihi_d = new double*[degp1];
   double **urhslohihi_d = new double*[degp1];
   double **urhshilohi_d = new double*[degp1];
   double **urhslolohi_d = new double*[degp1];
   double **urhshihilo_d = new double*[degp1];
   double **urhslohilo_d = new double*[degp1];
   double **urhshilolo_d = new double*[degp1];
   double **urhslololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhshihihi_h[i] = new double[dim];
      urhslohihi_h[i] = new double[dim];
      urhshilohi_h[i] = new double[dim];
      urhslolohi_h[i] = new double[dim];
      urhshihilo_h[i] = new double[dim];
      urhslohilo_h[i] = new double[dim];
      urhshilolo_h[i] = new double[dim];
      urhslololo_h[i] = new double[dim];
      urhshihihi_d[i] = new double[dim];
      urhslohihi_d[i] = new double[dim];
      urhshilohi_d[i] = new double[dim];
      urhslolohi_d[i] = new double[dim];
      urhshihilo_d[i] = new double[dim];
      urhslohilo_d[i] = new double[dim];
      urhshilolo_d[i] = new double[dim];
      urhslololo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihihi = new double*[dim];
   double **workmatlohihi = new double*[dim];
   double **workmathilohi = new double*[dim];
   double **workmatlolohi = new double*[dim];
   double **workmathihilo = new double*[dim];
   double **workmatlohilo = new double*[dim];
   double **workmathilolo = new double*[dim];
   double **workmatlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihihi[i] = new double[dim];
      workmatlohihi[i] = new double[dim];
      workmathilohi[i] = new double[dim];
      workmatlolohi[i] = new double[dim];
      workmathihilo[i] = new double[dim];
      workmatlohilo[i] = new double[dim];
      workmathilolo[i] = new double[dim];
      workmatlololo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechihihi = new double[dim];
   double *workveclohihi = new double[dim];
   double *workvechilohi = new double[dim];
   double *workveclolohi = new double[dim];
   double *workvechihilo = new double[dim];
   double *workveclohilo = new double[dim];
   double *workvechilolo = new double[dim];
   double *workveclololo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihihi = new double*[degp1];
   double **workrhslohihi = new double*[degp1];
   double **workrhshilohi = new double*[degp1];
   double **workrhslolohi = new double*[degp1];
   double **workrhshihilo = new double*[degp1];
   double **workrhslohilo = new double*[degp1];
   double **workrhshilolo = new double*[degp1];
   double **workrhslololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihihi[i] = new double[dim];
      workrhslohihi[i] = new double[dim];
      workrhshilohi[i] = new double[dim];
      workrhslolohi[i] = new double[dim];
      workrhshihilo[i] = new double[dim];
      workrhslohilo[i] = new double[dim];
      workrhshilolo[i] = new double[dim];
      workrhslololo[i] = new double[dim];
   }
   double **resvechihihi = new double*[degp1];
   double **resveclohihi = new double*[degp1];
   double **resvechilohi = new double*[degp1];
   double **resveclolohi = new double*[degp1];
   double **resvechihilo = new double*[degp1];
   double **resveclohilo = new double*[degp1];
   double **resvechilolo = new double*[degp1];
   double **resveclololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihihi[i] = new double[dim];
      resveclohihi[i] = new double[dim];
      resvechilohi[i] = new double[dim];
      resveclolohi[i] = new double[dim];
      resvechihilo[i] = new double[dim];
      resveclohilo[i] = new double[dim];
      resvechilolo[i] = new double[dim];
      resveclololo[i] = new double[dim];
   }
   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

   double **Qhihihi_h = new double*[dim];
   double **Qlohihi_h = new double*[dim];
   double **Qhilohi_h = new double*[dim];
   double **Qlolohi_h = new double*[dim];
   double **Qhihilo_h = new double*[dim];
   double **Qlohilo_h = new double*[dim];
   double **Qhilolo_h = new double*[dim];
   double **Qlololo_h = new double*[dim];
   double **Qhihihi_d = new double*[dim];
   double **Qlohihi_d = new double*[dim];
   double **Qhilohi_d = new double*[dim];
   double **Qlolohi_d = new double*[dim];
   double **Qhihilo_d = new double*[dim];
   double **Qlohilo_d = new double*[dim];
   double **Qhilolo_d = new double*[dim];
   double **Qlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qhihihi_h[i] = new double[dim];
      Qlohihi_h[i] = new double[dim];
      Qhilohi_h[i] = new double[dim];
      Qlolohi_h[i] = new double[dim];
      Qhihilo_h[i] = new double[dim];
      Qlohilo_h[i] = new double[dim];
      Qhilolo_h[i] = new double[dim];
      Qlololo_h[i] = new double[dim];
      Qhihihi_d[i] = new double[dim];
      Qlohihi_d[i] = new double[dim];
      Qhilohi_d[i] = new double[dim];
      Qlolohi_d[i] = new double[dim];
      Qhihilo_d[i] = new double[dim];
      Qlohilo_d[i] = new double[dim];
      Qhilolo_d[i] = new double[dim];
      Qlololo_d[i] = new double[dim];
   }
   double **Rhihihi_h = new double*[dim];
   double **Rlohihi_h = new double*[dim];
   double **Rhilohi_h = new double*[dim];
   double **Rlolohi_h = new double*[dim];
   double **Rhihilo_h = new double*[dim];
   double **Rlohilo_h = new double*[dim];
   double **Rhilolo_h = new double*[dim];
   double **Rlololo_h = new double*[dim];
   double **Rhihihi_d = new double*[dim];
   double **Rlohihi_d = new double*[dim];
   double **Rhilohi_d = new double*[dim];
   double **Rlolohi_d = new double*[dim];
   double **Rhihilo_d = new double*[dim];
   double **Rlohilo_d = new double*[dim];
   double **Rhilolo_d = new double*[dim];
   double **Rlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhihihi_h[i] = new double[dim];
      Rlohihi_h[i] = new double[dim];
      Rhilohi_h[i] = new double[dim];
      Rlolohi_h[i] = new double[dim];
      Rhihilo_h[i] = new double[dim];
      Rlohilo_h[i] = new double[dim];
      Rhilolo_h[i] = new double[dim];
      Rlololo_h[i] = new double[dim];
      Rhihihi_d[i] = new double[dim];
      Rlohihi_d[i] = new double[dim];
      Rhilohi_d[i] = new double[dim];
      Rlolohi_d[i] = new double[dim];
      Rhihilo_d[i] = new double[dim];
      Rlohilo_d[i] = new double[dim];
      Rhilolo_d[i] = new double[dim];
      Rlololo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl8_start_series_vector
      (dim,deg,inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
               inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h);
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihihi_d[i][j] = inputhihihi_h[i][j];
         inputlohihi_d[i][j] = inputlohihi_h[i][j];
         inputhilohi_d[i][j] = inputhilohi_h[i][j];
         inputlolohi_d[i][j] = inputlolohi_h[i][j];
         inputhihilo_d[i][j] = inputhihilo_h[i][j];
         inputlohilo_d[i][j] = inputlohilo_h[i][j];
         inputhilolo_d[i][j] = inputhilolo_h[i][j];
         inputlololo_d[i][j] = inputlololo_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihihi_h[i][0] << "  "
                            << inputlohihi_h[i][0] << endl;
         cout << "     " << inputhilohi_h[i][0] << "  "
                         << inputlolohi_h[i][0] << endl;
         cout << "     " << inputhihilo_h[i][0] << "  "
                         << inputlohilo_h[i][0] << endl;
         cout << "     " << inputhilolo_h[i][0] << "  "
                         << inputlololo_h[i][0] << endl;
      }
   }
   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;
/*
      dbl8_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi, inputlohihi, inputhilohi, inputlolohi,
          inputhihilo, inputlohilo, inputhilolo, inputlololo,
          outputhihihi,outputlohihi,outputhilohi,outputlolohi,
          outputhihilo,outputlohilo,outputhilolo,outputlololo,
          funvalhihihi,funvallohihi,funvalhilohi,funvallolohi,
          funvalhihilo,funvallohilo,funvalhilolo,funvallololo,
          jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
          jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          workrhshihihi,workrhslohihi,workrhshilohi,workrhslolohi,
          workrhshihilo,workrhslohilo,workrhshilolo,workrhslololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          ipvt,vrblvl);
 */
      dbl8_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
          funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
          funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
          funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
          funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
          rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
          rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
          rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
          rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
          urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
          urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
          urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
          urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          vrblvl,mode);
   }
   return 0;
}

int test_dbl8_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputrehihihi_h = new double*[dim];
   double **inputrelohihi_h = new double*[dim];
   double **inputrehilohi_h = new double*[dim];
   double **inputrelolohi_h = new double*[dim];
   double **inputrehihilo_h = new double*[dim];
   double **inputrelohilo_h = new double*[dim];
   double **inputrehilolo_h = new double*[dim];
   double **inputrelololo_h = new double*[dim];
   double **inputimhihihi_h = new double*[dim];
   double **inputimlohihi_h = new double*[dim];
   double **inputimhilohi_h = new double*[dim];
   double **inputimlolohi_h = new double*[dim];
   double **inputimhihilo_h = new double*[dim];
   double **inputimlohilo_h = new double*[dim];
   double **inputimhilolo_h = new double*[dim];
   double **inputimlololo_h = new double*[dim];
   double **inputrehihihi_d = new double*[dim];
   double **inputrelohihi_d = new double*[dim];
   double **inputrehilohi_d = new double*[dim];
   double **inputrelolohi_d = new double*[dim];
   double **inputrehihilo_d = new double*[dim];
   double **inputrelohilo_d = new double*[dim];
   double **inputrehilolo_d = new double*[dim];
   double **inputrelololo_d = new double*[dim];
   double **inputimhihihi_d = new double*[dim];
   double **inputimlohihi_d = new double*[dim];
   double **inputimhilohi_d = new double*[dim];
   double **inputimlolohi_d = new double*[dim];
   double **inputimhihilo_d = new double*[dim];
   double **inputimlohilo_d = new double*[dim];
   double **inputimhilolo_d = new double*[dim];
   double **inputimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputrehihihi_h[i] = new double[degp1];
       inputrelohihi_h[i] = new double[degp1];
       inputrehilohi_h[i] = new double[degp1];
       inputrelolohi_h[i] = new double[degp1];
       inputrehihilo_h[i] = new double[degp1];
       inputrelohilo_h[i] = new double[degp1];
       inputrehilolo_h[i] = new double[degp1];
       inputrelololo_h[i] = new double[degp1];
       inputimhihihi_h[i] = new double[degp1];
       inputimlohihi_h[i] = new double[degp1];
       inputimhilohi_h[i] = new double[degp1];
       inputimlolohi_h[i] = new double[degp1];
       inputimhihilo_h[i] = new double[degp1];
       inputimlohilo_h[i] = new double[degp1];
       inputimhilolo_h[i] = new double[degp1];
       inputimlololo_h[i] = new double[degp1];
       inputrehihihi_d[i] = new double[degp1];
       inputrelohihi_d[i] = new double[degp1];
       inputrehilohi_d[i] = new double[degp1];
       inputrelolohi_d[i] = new double[degp1];
       inputrehihilo_d[i] = new double[degp1];
       inputrelohilo_d[i] = new double[degp1];
       inputrehilolo_d[i] = new double[degp1];
       inputrelololo_d[i] = new double[degp1];
       inputimhihihi_d[i] = new double[degp1];
       inputimlohihi_d[i] = new double[degp1];
       inputimhilohi_d[i] = new double[degp1];
       inputimlolohi_d[i] = new double[degp1];
       inputimhihilo_d[i] = new double[degp1];
       inputimlohilo_d[i] = new double[degp1];
       inputimhilolo_d[i] = new double[degp1];
       inputimlololo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *accrehihihi = new double[degp1]; // accumulated power series
   double *accrelohihi = new double[degp1];
   double *accrehilohi = new double[degp1];
   double *accrelolohi = new double[degp1];
   double *accrehihilo = new double[degp1];
   double *accrelohilo = new double[degp1];
   double *accrehilolo = new double[degp1];
   double *accrelololo = new double[degp1];
   double *accimhihihi = new double[degp1];
   double *accimlohihi = new double[degp1];
   double *accimhilohi = new double[degp1];
   double *accimlolohi = new double[degp1];
   double *accimhihilo = new double[degp1];
   double *accimlohilo = new double[degp1];
   double *accimhilolo = new double[degp1];
   double *accimlololo = new double[degp1];
   double **cffrehihihi = new double*[dim]; // the coefficients of monomials
   double **cffrelohihi = new double*[dim];
   double **cffrehilohi = new double*[dim]; 
   double **cffrelolohi = new double*[dim];
   double **cffrehihilo = new double*[dim];
   double **cffrelohilo = new double*[dim];
   double **cffrehilolo = new double*[dim]; 
   double **cffrelololo = new double*[dim];
   double **cffimhihihi = new double*[dim]; 
   double **cffimlohihi = new double*[dim]; 
   double **cffimhilohi = new double*[dim]; 
   double **cffimlolohi = new double*[dim]; 
   double **cffimhihilo = new double*[dim]; 
   double **cffimlohilo = new double*[dim]; 
   double **cffimhilolo = new double*[dim]; 
   double **cffimlololo = new double*[dim]; 

   for(int i=0; i<dim; i++)
   {
      cffrehihihi[i] = new double[degp1];
      cffrelohihi[i] = new double[degp1];
      cffrehilohi[i] = new double[degp1];
      cffrelolohi[i] = new double[degp1];
      cffrehihilo[i] = new double[degp1];
      cffrelohilo[i] = new double[degp1];
      cffrehilolo[i] = new double[degp1];
      cffrelololo[i] = new double[degp1];
      cffimhihihi[i] = new double[degp1];
      cffimlohihi[i] = new double[degp1];
      cffimhilohi[i] = new double[degp1];
      cffimlolohi[i] = new double[degp1];
      cffimhihilo[i] = new double[degp1];
      cffimlohilo[i] = new double[degp1];
      cffimhilolo[i] = new double[degp1];
      cffimlololo[i] = new double[degp1];
   }
   double ***outputrehihihi_h = new double**[dim];
   double ***outputrelohihi_h = new double**[dim];
   double ***outputrehilohi_h = new double**[dim];
   double ***outputrelolohi_h = new double**[dim];
   double ***outputrehihilo_h = new double**[dim];
   double ***outputrelohilo_h = new double**[dim];
   double ***outputrehilolo_h = new double**[dim];
   double ***outputrelololo_h = new double**[dim];
   double ***outputimhihihi_h = new double**[dim];
   double ***outputimlohihi_h = new double**[dim];
   double ***outputimhilohi_h = new double**[dim];
   double ***outputimlolohi_h = new double**[dim];
   double ***outputimhihilo_h = new double**[dim];
   double ***outputimlohilo_h = new double**[dim];
   double ***outputimhilolo_h = new double**[dim];
   double ***outputimlololo_h = new double**[dim];
   double ***outputrehihihi_d = new double**[dim];
   double ***outputrelohihi_d = new double**[dim];
   double ***outputrehilohi_d = new double**[dim];
   double ***outputrelolohi_d = new double**[dim];
   double ***outputrehihilo_d = new double**[dim];
   double ***outputrelohilo_d = new double**[dim];
   double ***outputrehilolo_d = new double**[dim];
   double ***outputrelololo_d = new double**[dim];
   double ***outputimhihihi_d = new double**[dim];
   double ***outputimlohihi_d = new double**[dim];
   double ***outputimhilohi_d = new double**[dim];
   double ***outputimlolohi_d = new double**[dim];
   double ***outputimhihilo_d = new double**[dim];
   double ***outputimlohilo_d = new double**[dim];
   double ***outputimhilolo_d = new double**[dim];
   double ***outputimlololo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehihihi_h[i] = new double*[dim+1];
      outputrelohihi_h[i] = new double*[dim+1];
      outputrehilohi_h[i] = new double*[dim+1];
      outputrelolohi_h[i] = new double*[dim+1];
      outputrehihilo_h[i] = new double*[dim+1];
      outputrelohilo_h[i] = new double*[dim+1];
      outputrehilolo_h[i] = new double*[dim+1];
      outputrelololo_h[i] = new double*[dim+1];
      outputimhihihi_h[i] = new double*[dim+1];
      outputimlohihi_h[i] = new double*[dim+1];
      outputimhilohi_h[i] = new double*[dim+1];
      outputimlolohi_h[i] = new double*[dim+1];
      outputimhihilo_h[i] = new double*[dim+1];
      outputimlohilo_h[i] = new double*[dim+1];
      outputimhilolo_h[i] = new double*[dim+1];
      outputimlololo_h[i] = new double*[dim+1];
      outputrehihihi_d[i] = new double*[dim+1];
      outputrelohihi_d[i] = new double*[dim+1];
      outputrehilohi_d[i] = new double*[dim+1];
      outputrelolohi_d[i] = new double*[dim+1];
      outputrehihilo_d[i] = new double*[dim+1];
      outputrelohilo_d[i] = new double*[dim+1];
      outputrehilolo_d[i] = new double*[dim+1];
      outputrelololo_d[i] = new double*[dim+1];
      outputimhihihi_d[i] = new double*[dim+1];
      outputimlohihi_d[i] = new double*[dim+1];
      outputimhilohi_d[i] = new double*[dim+1];
      outputimlolohi_d[i] = new double*[dim+1];
      outputimhihilo_d[i] = new double*[dim+1];
      outputimlohilo_d[i] = new double*[dim+1];
      outputimhilolo_d[i] = new double*[dim+1];
      outputimlololo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehihihi_h[i][j] = new double[degp1];
         outputrelohihi_h[i][j] = new double[degp1];
         outputrehilohi_h[i][j] = new double[degp1];
         outputrelolohi_h[i][j] = new double[degp1];
         outputrehihilo_h[i][j] = new double[degp1];
         outputrelohilo_h[i][j] = new double[degp1];
         outputrehilolo_h[i][j] = new double[degp1];
         outputrelololo_h[i][j] = new double[degp1];
         outputimhihihi_h[i][j] = new double[degp1];
         outputimlohihi_h[i][j] = new double[degp1];
         outputimhilohi_h[i][j] = new double[degp1];
         outputimlolohi_h[i][j] = new double[degp1];
         outputimhihilo_h[i][j] = new double[degp1];
         outputimlohilo_h[i][j] = new double[degp1];
         outputimhilolo_h[i][j] = new double[degp1];
         outputimlololo_h[i][j] = new double[degp1];
         outputrehihihi_d[i][j] = new double[degp1];
         outputrelohihi_d[i][j] = new double[degp1];
         outputrehilohi_d[i][j] = new double[degp1];
         outputrelolohi_d[i][j] = new double[degp1];
         outputrehihilo_d[i][j] = new double[degp1];
         outputrelohilo_d[i][j] = new double[degp1];
         outputrehilolo_d[i][j] = new double[degp1];
         outputrelololo_d[i][j] = new double[degp1];
         outputimhihihi_d[i][j] = new double[degp1];
         outputimlohihi_d[i][j] = new double[degp1];
         outputimhilohi_d[i][j] = new double[degp1];
         outputimlolohi_d[i][j] = new double[degp1];
         outputimhihilo_d[i][j] = new double[degp1];
         outputimlohilo_d[i][j] = new double[degp1];
         outputimhilolo_d[i][j] = new double[degp1];
         outputimlololo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalrehihihi_h = new double*[dim];
   double **funvalrelohihi_h = new double*[dim];
   double **funvalrehilohi_h = new double*[dim];
   double **funvalrelolohi_h = new double*[dim];
   double **funvalrehihilo_h = new double*[dim];
   double **funvalrelohilo_h = new double*[dim];
   double **funvalrehilolo_h = new double*[dim];
   double **funvalrelololo_h = new double*[dim];
   double **funvalimhihihi_h = new double*[dim];
   double **funvalimlohihi_h = new double*[dim];
   double **funvalimhilohi_h = new double*[dim];
   double **funvalimlolohi_h = new double*[dim];
   double **funvalimhihilo_h = new double*[dim];
   double **funvalimlohilo_h = new double*[dim];
   double **funvalimhilolo_h = new double*[dim];
   double **funvalimlololo_h = new double*[dim];
   double **funvalrehihihi_d = new double*[dim];
   double **funvalrelohihi_d = new double*[dim];
   double **funvalrehilohi_d = new double*[dim];
   double **funvalrelolohi_d = new double*[dim];
   double **funvalrehihilo_d = new double*[dim];
   double **funvalrelohilo_d = new double*[dim];
   double **funvalrehilolo_d = new double*[dim];
   double **funvalrelololo_d = new double*[dim];
   double **funvalimhihihi_d = new double*[dim];
   double **funvalimlohihi_d = new double*[dim];
   double **funvalimhilohi_d = new double*[dim];
   double **funvalimlolohi_d = new double*[dim];
   double **funvalimhihilo_d = new double*[dim];
   double **funvalimlohilo_d = new double*[dim];
   double **funvalimhilolo_d = new double*[dim];
   double **funvalimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalrehihihi_h[i] = new double[degp1];
      funvalrelohihi_h[i] = new double[degp1];
      funvalrehilohi_h[i] = new double[degp1];
      funvalrelolohi_h[i] = new double[degp1];
      funvalrehihilo_h[i] = new double[degp1];
      funvalrelohilo_h[i] = new double[degp1];
      funvalrehilolo_h[i] = new double[degp1];
      funvalrelololo_h[i] = new double[degp1];
      funvalimhihihi_h[i] = new double[degp1];
      funvalimlohihi_h[i] = new double[degp1];
      funvalimhilohi_h[i] = new double[degp1];
      funvalimlolohi_h[i] = new double[degp1];
      funvalimhihilo_h[i] = new double[degp1];
      funvalimlohilo_h[i] = new double[degp1];
      funvalimhilolo_h[i] = new double[degp1];
      funvalimlololo_h[i] = new double[degp1];
      funvalrehihihi_d[i] = new double[degp1];
      funvalrelohihi_d[i] = new double[degp1];
      funvalrehilohi_d[i] = new double[degp1];
      funvalrelolohi_d[i] = new double[degp1];
      funvalrehihilo_d[i] = new double[degp1];
      funvalrelohilo_d[i] = new double[degp1];
      funvalrehilolo_d[i] = new double[degp1];
      funvalrelololo_d[i] = new double[degp1];
      funvalimhihihi_d[i] = new double[degp1];
      funvalimlohihi_d[i] = new double[degp1];
      funvalimhilohi_d[i] = new double[degp1];
      funvalimlolohi_d[i] = new double[degp1];
      funvalimhihilo_d[i] = new double[degp1];
      funvalimlohilo_d[i] = new double[degp1];
      funvalimhilolo_d[i] = new double[degp1];
      funvalimlololo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalrehihihi_h = new double**[degp1];
   double ***jacvalrelohihi_h = new double**[degp1];
   double ***jacvalrehilohi_h = new double**[degp1];
   double ***jacvalrelolohi_h = new double**[degp1];
   double ***jacvalrehihilo_h = new double**[degp1];
   double ***jacvalrelohilo_h = new double**[degp1];
   double ***jacvalrehilolo_h = new double**[degp1];
   double ***jacvalrelololo_h = new double**[degp1];
   double ***jacvalimhihihi_h = new double**[degp1];
   double ***jacvalimlohihi_h = new double**[degp1];
   double ***jacvalimhilohi_h = new double**[degp1];
   double ***jacvalimlolohi_h = new double**[degp1];
   double ***jacvalimhihilo_h = new double**[degp1];
   double ***jacvalimlohilo_h = new double**[degp1];
   double ***jacvalimhilolo_h = new double**[degp1];
   double ***jacvalimlololo_h = new double**[degp1];
   double ***jacvalrehihihi_d = new double**[degp1];
   double ***jacvalrelohihi_d = new double**[degp1];
   double ***jacvalrehilohi_d = new double**[degp1];
   double ***jacvalrelolohi_d = new double**[degp1];
   double ***jacvalrehihilo_d = new double**[degp1];
   double ***jacvalrelohilo_d = new double**[degp1];
   double ***jacvalrehilolo_d = new double**[degp1];
   double ***jacvalrelololo_d = new double**[degp1];
   double ***jacvalimhihihi_d = new double**[degp1];
   double ***jacvalimlohihi_d = new double**[degp1];
   double ***jacvalimhilohi_d = new double**[degp1];
   double ***jacvalimlolohi_d = new double**[degp1];
   double ***jacvalimhihilo_d = new double**[degp1];
   double ***jacvalimlohilo_d = new double**[degp1];
   double ***jacvalimhilolo_d = new double**[degp1];
   double ***jacvalimlololo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalrehihihi_h[i] = new double*[dim];
      jacvalrelohihi_h[i] = new double*[dim];
      jacvalrehilohi_h[i] = new double*[dim];
      jacvalrelolohi_h[i] = new double*[dim];
      jacvalrehihilo_h[i] = new double*[dim];
      jacvalrelohilo_h[i] = new double*[dim];
      jacvalrehilolo_h[i] = new double*[dim];
      jacvalrelololo_h[i] = new double*[dim];
      jacvalimhihihi_h[i] = new double*[dim];
      jacvalimlohihi_h[i] = new double*[dim];
      jacvalimhilohi_h[i] = new double*[dim];
      jacvalimlolohi_h[i] = new double*[dim];
      jacvalimhihilo_h[i] = new double*[dim];
      jacvalimlohilo_h[i] = new double*[dim];
      jacvalimhilolo_h[i] = new double*[dim];
      jacvalimlololo_h[i] = new double*[dim];
      jacvalrehihihi_d[i] = new double*[dim];
      jacvalrelohihi_d[i] = new double*[dim];
      jacvalrehilohi_d[i] = new double*[dim];
      jacvalrelolohi_d[i] = new double*[dim];
      jacvalrehihilo_d[i] = new double*[dim];
      jacvalrelohilo_d[i] = new double*[dim];
      jacvalrehilolo_d[i] = new double*[dim];
      jacvalrelololo_d[i] = new double*[dim];
      jacvalimhihihi_d[i] = new double*[dim];
      jacvalimlohihi_d[i] = new double*[dim];
      jacvalimhilohi_d[i] = new double*[dim];
      jacvalimlolohi_d[i] = new double*[dim];
      jacvalimhihilo_d[i] = new double*[dim];
      jacvalimlohilo_d[i] = new double*[dim];
      jacvalimhilolo_d[i] = new double*[dim];
      jacvalimlololo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalrehihihi_h[i][j] = new double[dim];
         jacvalrelohihi_h[i][j] = new double[dim];
         jacvalrehilohi_h[i][j] = new double[dim];
         jacvalrelolohi_h[i][j] = new double[dim];
         jacvalrehihilo_h[i][j] = new double[dim];
         jacvalrelohilo_h[i][j] = new double[dim];
         jacvalrehilolo_h[i][j] = new double[dim];
         jacvalrelololo_h[i][j] = new double[dim];
         jacvalimhihihi_h[i][j] = new double[dim];
         jacvalimlohihi_h[i][j] = new double[dim];
         jacvalimhilohi_h[i][j] = new double[dim];
         jacvalimlolohi_h[i][j] = new double[dim];
         jacvalimhihilo_h[i][j] = new double[dim];
         jacvalimlohilo_h[i][j] = new double[dim];
         jacvalimhilolo_h[i][j] = new double[dim];
         jacvalimlololo_h[i][j] = new double[dim];
         jacvalrehihihi_d[i][j] = new double[dim];
         jacvalrelohihi_d[i][j] = new double[dim];
         jacvalrehilohi_d[i][j] = new double[dim];
         jacvalrelolohi_d[i][j] = new double[dim];
         jacvalrehihilo_d[i][j] = new double[dim];
         jacvalrelohilo_d[i][j] = new double[dim];
         jacvalrehilolo_d[i][j] = new double[dim];
         jacvalrelololo_d[i][j] = new double[dim];
         jacvalimhihihi_d[i][j] = new double[dim];
         jacvalimlohihi_d[i][j] = new double[dim];
         jacvalimhilohi_d[i][j] = new double[dim];
         jacvalimlolohi_d[i][j] = new double[dim];
         jacvalimhihilo_d[i][j] = new double[dim];
         jacvalimlohilo_d[i][j] = new double[dim];
         jacvalimhilolo_d[i][j] = new double[dim];
         jacvalimlololo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solrehihihi_h = new double*[degp1];
   double **solrelohihi_h = new double*[degp1];
   double **solrehilohi_h = new double*[degp1];
   double **solrelolohi_h = new double*[degp1];
   double **solrehihilo_h = new double*[degp1];
   double **solrelohilo_h = new double*[degp1];
   double **solrehilolo_h = new double*[degp1];
   double **solrelololo_h = new double*[degp1];
   double **solimhihihi_h = new double*[degp1];
   double **solimlohihi_h = new double*[degp1];
   double **solimhilohi_h = new double*[degp1];
   double **solimlolohi_h = new double*[degp1];
   double **solimhihilo_h = new double*[degp1];
   double **solimlohilo_h = new double*[degp1];
   double **solimhilolo_h = new double*[degp1];
   double **solimlololo_h = new double*[degp1];
   double **solrehihihi_d = new double*[degp1];
   double **solrelohihi_d = new double*[degp1];
   double **solrehilohi_d = new double*[degp1];
   double **solrelolohi_d = new double*[degp1];
   double **solrehihilo_d = new double*[degp1];
   double **solrelohilo_d = new double*[degp1];
   double **solrehilolo_d = new double*[degp1];
   double **solrelololo_d = new double*[degp1];
   double **solimhihihi_d = new double*[degp1];
   double **solimlohihi_d = new double*[degp1];
   double **solimhilohi_d = new double*[degp1];
   double **solimlolohi_d = new double*[degp1];
   double **solimhihilo_d = new double*[degp1];
   double **solimlohilo_d = new double*[degp1];
   double **solimhilolo_d = new double*[degp1];
   double **solimlololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++) 
   {
      solrehihihi_h[i] = new double[dim];
      solrelohihi_h[i] = new double[dim];
      solrehilohi_h[i] = new double[dim];
      solrelolohi_h[i] = new double[dim];
      solrehihilo_h[i] = new double[dim];
      solrelohilo_h[i] = new double[dim];
      solrehilolo_h[i] = new double[dim];
      solrelololo_h[i] = new double[dim];
      solimhihihi_h[i] = new double[dim];
      solimlohihi_h[i] = new double[dim];
      solimhilohi_h[i] = new double[dim];
      solimlolohi_h[i] = new double[dim];
      solimhihilo_h[i] = new double[dim];
      solimlohilo_h[i] = new double[dim];
      solimhilolo_h[i] = new double[dim];
      solimlololo_h[i] = new double[dim];
      solrehihihi_d[i] = new double[dim];
      solrelohihi_d[i] = new double[dim];
      solrehilohi_d[i] = new double[dim];
      solrelolohi_d[i] = new double[dim];
      solrehihilo_d[i] = new double[dim];
      solrelohilo_d[i] = new double[dim];
      solrehilolo_d[i] = new double[dim];
      solrelololo_d[i] = new double[dim];
      solimhihihi_d[i] = new double[dim];
      solimlohihi_d[i] = new double[dim];
      solimhilohi_d[i] = new double[dim];
      solimlolohi_d[i] = new double[dim];
      solimhihilo_d[i] = new double[dim];
      solimlohilo_d[i] = new double[dim];
      solimhilolo_d[i] = new double[dim];
      solimlololo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsrehihihi_h = new double*[degp1];
   double **rhsrelohihi_h = new double*[degp1];
   double **rhsrehilohi_h = new double*[degp1];
   double **rhsrelolohi_h = new double*[degp1];
   double **rhsrehihilo_h = new double*[degp1];
   double **rhsrelohilo_h = new double*[degp1];
   double **rhsrehilolo_h = new double*[degp1];
   double **rhsrelololo_h = new double*[degp1];
   double **rhsimhihihi_h = new double*[degp1];
   double **rhsimlohihi_h = new double*[degp1];
   double **rhsimhilohi_h = new double*[degp1];
   double **rhsimlolohi_h = new double*[degp1];
   double **rhsimhihilo_h = new double*[degp1];
   double **rhsimlohilo_h = new double*[degp1];
   double **rhsimhilolo_h = new double*[degp1];
   double **rhsimlololo_h = new double*[degp1];
   double **rhsrehihihi_d = new double*[degp1];
   double **rhsrelohihi_d = new double*[degp1];
   double **rhsrehilohi_d = new double*[degp1];
   double **rhsrelolohi_d = new double*[degp1];
   double **rhsrehihilo_d = new double*[degp1];
   double **rhsrelohilo_d = new double*[degp1];
   double **rhsrehilolo_d = new double*[degp1];
   double **rhsrelololo_d = new double*[degp1];
   double **rhsimhihihi_d = new double*[degp1];
   double **rhsimlohihi_d = new double*[degp1];
   double **rhsimhilohi_d = new double*[degp1];
   double **rhsimlolohi_d = new double*[degp1];
   double **rhsimhihilo_d = new double*[degp1];
   double **rhsimlohilo_d = new double*[degp1];
   double **rhsimhilolo_d = new double*[degp1];
   double **rhsimlololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsrehihihi_h[i] = new double[dim];
      rhsrelohihi_h[i] = new double[dim];
      rhsrehilohi_h[i] = new double[dim];
      rhsrelolohi_h[i] = new double[dim];
      rhsrehihilo_h[i] = new double[dim];
      rhsrelohilo_h[i] = new double[dim];
      rhsrehilolo_h[i] = new double[dim];
      rhsrelololo_h[i] = new double[dim];
      rhsimhihihi_h[i] = new double[dim];
      rhsimlohihi_h[i] = new double[dim];
      rhsimhilohi_h[i] = new double[dim];
      rhsimlolohi_h[i] = new double[dim];
      rhsimhihilo_h[i] = new double[dim];
      rhsimlohilo_h[i] = new double[dim];
      rhsimhilolo_h[i] = new double[dim];
      rhsimlololo_h[i] = new double[dim];
      rhsrehihihi_d[i] = new double[dim];
      rhsrelohihi_d[i] = new double[dim];
      rhsrehilohi_d[i] = new double[dim];
      rhsrelolohi_d[i] = new double[dim];
      rhsrehihilo_d[i] = new double[dim];
      rhsrelohilo_d[i] = new double[dim];
      rhsrehilolo_d[i] = new double[dim];
      rhsrelololo_d[i] = new double[dim];
      rhsimhihihi_d[i] = new double[dim];
      rhsimlohihi_d[i] = new double[dim];
      rhsimhilohi_d[i] = new double[dim];
      rhsimlolohi_d[i] = new double[dim];
      rhsimhihilo_d[i] = new double[dim];
      rhsimlohilo_d[i] = new double[dim];
      rhsimhilolo_d[i] = new double[dim];
      rhsimlololo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmatrehihihi = new double*[dim];
   double **workmatrelohihi = new double*[dim];
   double **workmatrehilohi = new double*[dim];
   double **workmatrelolohi = new double*[dim];
   double **workmatrehihilo = new double*[dim];
   double **workmatrelohilo = new double*[dim];
   double **workmatrehilolo = new double*[dim];
   double **workmatrelololo = new double*[dim];
   double **workmatimhihihi = new double*[dim];
   double **workmatimlohihi = new double*[dim];
   double **workmatimhilohi = new double*[dim];
   double **workmatimlolohi = new double*[dim];
   double **workmatimhihilo = new double*[dim];
   double **workmatimlohilo = new double*[dim];
   double **workmatimhilolo = new double*[dim];
   double **workmatimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmatrehihihi[i] = new double[dim];
      workmatrelohihi[i] = new double[dim];
      workmatrehilohi[i] = new double[dim];
      workmatrelolohi[i] = new double[dim];
      workmatrehihilo[i] = new double[dim];
      workmatrelohilo[i] = new double[dim];
      workmatrehilolo[i] = new double[dim];
      workmatrelololo[i] = new double[dim];
      workmatimhihihi[i] = new double[dim];
      workmatimlohihi[i] = new double[dim];
      workmatimhilohi[i] = new double[dim];
      workmatimlolohi[i] = new double[dim];
      workmatimhihilo[i] = new double[dim];
      workmatimlohilo[i] = new double[dim];
      workmatimhilolo[i] = new double[dim];
      workmatimlololo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvecrehihihi = new double[dim];
   double *workvecrelohihi = new double[dim];
   double *workvecrehilohi = new double[dim];
   double *workvecrelolohi = new double[dim];
   double *workvecrehihilo = new double[dim];
   double *workvecrelohilo = new double[dim];
   double *workvecrehilolo = new double[dim];
   double *workvecrelololo = new double[dim];
   double *workvecimhihihi = new double[dim];
   double *workvecimlohihi = new double[dim];
   double *workvecimhilohi = new double[dim];
   double *workvecimlolohi = new double[dim];
   double *workvecimhihilo = new double[dim];
   double *workvecimlohilo = new double[dim];
   double *workvecimhilolo = new double[dim];
   double *workvecimlololo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **urhsrehihihi_h = new double*[degp1];
   double **urhsrelohihi_h = new double*[degp1];
   double **urhsrehilohi_h = new double*[degp1];
   double **urhsrelolohi_h = new double*[degp1];
   double **urhsrehihilo_h = new double*[degp1];
   double **urhsrelohilo_h = new double*[degp1];
   double **urhsrehilolo_h = new double*[degp1];
   double **urhsrelololo_h = new double*[degp1];
   double **urhsimhihihi_h = new double*[degp1];
   double **urhsimlohihi_h = new double*[degp1];
   double **urhsimhilohi_h = new double*[degp1];
   double **urhsimlolohi_h = new double*[degp1];
   double **urhsimhihilo_h = new double*[degp1];
   double **urhsimlohilo_h = new double*[degp1];
   double **urhsimhilolo_h = new double*[degp1];
   double **urhsimlololo_h = new double*[degp1];
   double **urhsrehihihi_d = new double*[degp1];
   double **urhsrelohihi_d = new double*[degp1];
   double **urhsrehilohi_d = new double*[degp1];
   double **urhsrelolohi_d = new double*[degp1];
   double **urhsrehihilo_d = new double*[degp1];
   double **urhsrelohilo_d = new double*[degp1];
   double **urhsrehilolo_d = new double*[degp1];
   double **urhsrelololo_d = new double*[degp1];
   double **urhsimhihihi_d = new double*[degp1];
   double **urhsimlohihi_d = new double*[degp1];
   double **urhsimhilohi_d = new double*[degp1];
   double **urhsimlolohi_d = new double*[degp1];
   double **urhsimhihilo_d = new double*[degp1];
   double **urhsimlohilo_d = new double*[degp1];
   double **urhsimhilolo_d = new double*[degp1];
   double **urhsimlololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhsrehihihi_h[i] = new double[dim];
      urhsrelohihi_h[i] = new double[dim];
      urhsrehilohi_h[i] = new double[dim];
      urhsrelolohi_h[i] = new double[dim];
      urhsrehihilo_h[i] = new double[dim];
      urhsrelohilo_h[i] = new double[dim];
      urhsrehilolo_h[i] = new double[dim];
      urhsrelololo_h[i] = new double[dim];
      urhsimhihihi_h[i] = new double[dim];
      urhsimlohihi_h[i] = new double[dim];
      urhsimhilohi_h[i] = new double[dim];
      urhsimlolohi_h[i] = new double[dim];
      urhsimhihilo_h[i] = new double[dim];
      urhsimlohilo_h[i] = new double[dim];
      urhsimhilolo_h[i] = new double[dim];
      urhsimlololo_h[i] = new double[dim];
      urhsrehihihi_d[i] = new double[dim];
      urhsrelohihi_d[i] = new double[dim];
      urhsrehilohi_d[i] = new double[dim];
      urhsrelolohi_d[i] = new double[dim];
      urhsrehihilo_d[i] = new double[dim];
      urhsrelohilo_d[i] = new double[dim];
      urhsrehilolo_d[i] = new double[dim];
      urhsrelololo_d[i] = new double[dim];
      urhsimhihihi_d[i] = new double[dim];
      urhsimlohihi_d[i] = new double[dim];
      urhsimhilohi_d[i] = new double[dim];
      urhsimlolohi_d[i] = new double[dim];
      urhsimhihilo_d[i] = new double[dim];
      urhsimlohilo_d[i] = new double[dim];
      urhsimhilolo_d[i] = new double[dim];
      urhsimlololo_d[i] = new double[dim];
   }
   double **resvecrehihihi = new double*[degp1];
   double **resvecrelohihi = new double*[degp1];
   double **resvecrehilohi = new double*[degp1];
   double **resvecrelolohi = new double*[degp1];
   double **resvecrehihilo = new double*[degp1];
   double **resvecrelohilo = new double*[degp1];
   double **resvecrehilolo = new double*[degp1];
   double **resvecrelololo = new double*[degp1];
   double **resvecimhihihi = new double*[degp1];
   double **resvecimlohihi = new double*[degp1];
   double **resvecimhilohi = new double*[degp1];
   double **resvecimlolohi = new double*[degp1];
   double **resvecimhihilo = new double*[degp1];
   double **resvecimlohilo = new double*[degp1];
   double **resvecimhilolo = new double*[degp1];
   double **resvecimlololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecrehihihi[i] = new double[dim];
      resvecrelohihi[i] = new double[dim];
      resvecrehilohi[i] = new double[dim];
      resvecrelolohi[i] = new double[dim];
      resvecrehihilo[i] = new double[dim];
      resvecrelohilo[i] = new double[dim];
      resvecrehilolo[i] = new double[dim];
      resvecrelololo[i] = new double[dim];
      resvecimhihihi[i] = new double[dim];
      resvecimlohihi[i] = new double[dim];
      resvecimhilohi[i] = new double[dim];
      resvecimlolohi[i] = new double[dim];
      resvecimhihilo[i] = new double[dim];
      resvecimlohilo[i] = new double[dim];
      resvecimhilolo[i] = new double[dim];
      resvecimlololo[i] = new double[dim];
   }
   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;
   double **Qrehihihi_h = new double*[dim];
   double **Qrelohihi_h = new double*[dim];
   double **Qrehilohi_h = new double*[dim];
   double **Qrelolohi_h = new double*[dim];
   double **Qrehihilo_h = new double*[dim];
   double **Qrelohilo_h = new double*[dim];
   double **Qrehilolo_h = new double*[dim];
   double **Qrelololo_h = new double*[dim];
   double **Qimhihihi_h = new double*[dim];
   double **Qimlohihi_h = new double*[dim];
   double **Qimhilohi_h = new double*[dim];
   double **Qimlolohi_h = new double*[dim];
   double **Qimhihilo_h = new double*[dim];
   double **Qimlohilo_h = new double*[dim];
   double **Qimhilolo_h = new double*[dim];
   double **Qimlololo_h = new double*[dim];
   double **Qrehihihi_d = new double*[dim];
   double **Qrelohihi_d = new double*[dim];
   double **Qrehilohi_d = new double*[dim];
   double **Qrelolohi_d = new double*[dim];
   double **Qrehihilo_d = new double*[dim];
   double **Qrelohilo_d = new double*[dim];
   double **Qrehilolo_d = new double*[dim];
   double **Qrelololo_d = new double*[dim];
   double **Qimhihihi_d = new double*[dim];
   double **Qimlohihi_d = new double*[dim];
   double **Qimhilohi_d = new double*[dim];
   double **Qimlolohi_d = new double*[dim];
   double **Qimhihilo_d = new double*[dim];
   double **Qimlohilo_d = new double*[dim];
   double **Qimhilolo_d = new double*[dim];
   double **Qimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qrehihihi_h[i] = new double[dim];
      Qrelohihi_h[i] = new double[dim];
      Qrehilohi_h[i] = new double[dim];
      Qrelolohi_h[i] = new double[dim];
      Qrehihilo_h[i] = new double[dim];
      Qrelohilo_h[i] = new double[dim];
      Qrehilolo_h[i] = new double[dim];
      Qrelololo_h[i] = new double[dim];
      Qimhihihi_h[i] = new double[dim];
      Qimlohihi_h[i] = new double[dim];
      Qimhilohi_h[i] = new double[dim];
      Qimlolohi_h[i] = new double[dim];
      Qimhihilo_h[i] = new double[dim];
      Qimlohilo_h[i] = new double[dim];
      Qimhilolo_h[i] = new double[dim];
      Qimlololo_h[i] = new double[dim];
      Qrehihihi_d[i] = new double[dim];
      Qrelohihi_d[i] = new double[dim];
      Qrehilohi_d[i] = new double[dim];
      Qrelolohi_d[i] = new double[dim];
      Qrehihilo_d[i] = new double[dim];
      Qrelohilo_d[i] = new double[dim];
      Qrehilolo_d[i] = new double[dim];
      Qrelololo_d[i] = new double[dim];
      Qimhihihi_d[i] = new double[dim];
      Qimlohihi_d[i] = new double[dim];
      Qimhilohi_d[i] = new double[dim];
      Qimlolohi_d[i] = new double[dim];
      Qimhihilo_d[i] = new double[dim];
      Qimlohilo_d[i] = new double[dim];
      Qimhilolo_d[i] = new double[dim];
      Qimlololo_d[i] = new double[dim];
   }
   double **Rrehihihi_h = new double*[dim];
   double **Rrelohihi_h = new double*[dim];
   double **Rrehilohi_h = new double*[dim];
   double **Rrelolohi_h = new double*[dim];
   double **Rrehihilo_h = new double*[dim];
   double **Rrelohilo_h = new double*[dim];
   double **Rrehilolo_h = new double*[dim];
   double **Rrelololo_h = new double*[dim];
   double **Rimhihihi_h = new double*[dim];
   double **Rimlohihi_h = new double*[dim];
   double **Rimhilohi_h = new double*[dim];
   double **Rimlolohi_h = new double*[dim];
   double **Rimhihilo_h = new double*[dim];
   double **Rimlohilo_h = new double*[dim];
   double **Rimhilolo_h = new double*[dim];
   double **Rimlololo_h = new double*[dim];
   double **Rrehihihi_d = new double*[dim];
   double **Rrelohihi_d = new double*[dim];
   double **Rrehilohi_d = new double*[dim];
   double **Rrelolohi_d = new double*[dim];
   double **Rrehihilo_d = new double*[dim];
   double **Rrelohilo_d = new double*[dim];
   double **Rrehilolo_d = new double*[dim];
   double **Rrelololo_d = new double*[dim];
   double **Rimhihihi_d = new double*[dim];
   double **Rimlohihi_d = new double*[dim];
   double **Rimhilohi_d = new double*[dim];
   double **Rimlolohi_d = new double*[dim];
   double **Rimhihilo_d = new double*[dim];
   double **Rimlohilo_d = new double*[dim];
   double **Rimhilolo_d = new double*[dim];
   double **Rimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rrehihihi_h[i] = new double[dim];
      Rrelohihi_h[i] = new double[dim];
      Rrehilohi_h[i] = new double[dim];
      Rrelolohi_h[i] = new double[dim];
      Rrehihilo_h[i] = new double[dim];
      Rrelohilo_h[i] = new double[dim];
      Rrehilolo_h[i] = new double[dim];
      Rrelololo_h[i] = new double[dim];
      Rimhihihi_h[i] = new double[dim];
      Rimlohihi_h[i] = new double[dim];
      Rimhilohi_h[i] = new double[dim];
      Rimlolohi_h[i] = new double[dim];
      Rimhihilo_h[i] = new double[dim];
      Rimlohilo_h[i] = new double[dim];
      Rimhilolo_h[i] = new double[dim];
      Rimlololo_h[i] = new double[dim];
      Rrehihihi_d[i] = new double[dim];
      Rrelohihi_d[i] = new double[dim];
      Rrehilohi_d[i] = new double[dim];
      Rrelolohi_d[i] = new double[dim];
      Rrehihilo_d[i] = new double[dim];
      Rrelohilo_d[i] = new double[dim];
      Rrehilolo_d[i] = new double[dim];
      Rrelololo_d[i] = new double[dim];
      Rimhihihi_d[i] = new double[dim];
      Rimlohihi_d[i] = new double[dim];
      Rimhilohi_d[i] = new double[dim];
      Rimlolohi_d[i] = new double[dim];
      Rimhihilo_d[i] = new double[dim];
      Rimlohilo_d[i] = new double[dim];
      Rimhilolo_d[i] = new double[dim];
      Rimlololo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   cmplx8_start_series_vector
      (dim,deg,
       inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
       inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
       inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
       inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehihihi_d[i][j] = inputrehihihi_h[i][j];
         inputrelohihi_d[i][j] = inputrelohihi_h[i][j];
         inputrehilohi_d[i][j] = inputrehilohi_h[i][j];
         inputrelolohi_d[i][j] = inputrelolohi_h[i][j];
         inputrehihilo_d[i][j] = inputrehihilo_h[i][j];
         inputrelohilo_d[i][j] = inputrelohilo_h[i][j];
         inputrehilolo_d[i][j] = inputrehilolo_h[i][j];
         inputrelololo_d[i][j] = inputrelololo_h[i][j];
         inputimhihihi_d[i][j] = inputimhihihi_h[i][j];
         inputimlohihi_d[i][j] = inputimlohihi_h[i][j];
         inputimhilohi_d[i][j] = inputimhilohi_h[i][j];
         inputimlolohi_d[i][j] = inputimlolohi_h[i][j];
         inputimhihilo_d[i][j] = inputimhihilo_h[i][j];
         inputimlohilo_d[i][j] = inputimlohilo_h[i][j];
         inputimhilolo_d[i][j] = inputimhilolo_h[i][j];
         inputimlololo_d[i][j] = inputimlololo_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << inputrehihihi_h[i][0] << "  "
              << inputrelohihi_h[i][0] << endl << "  "
              << inputrehilohi_h[i][0] << "  "
              << inputrelolohi_h[i][0] << endl << "  "
              << inputrehihilo_h[i][0] << "  "
              << inputrelohilo_h[i][0] << endl << "  "
              << inputrehilolo_h[i][0] << "  "
              << inputrelololo_h[i][0] << endl << "  "
              << inputimhihihi_h[i][0] << "  "
              << inputimlohihi_h[i][0] << endl << "  "
              << inputimhilohi_h[i][0] << "  "
              << inputimlolohi_h[i][0] << endl << "  "
              << inputimhihilo_h[i][0] << "  "
              << inputimlohilo_h[i][0] << endl << "  "
              << inputimhilolo_h[i][0] << "  "
              << inputimlololo_h[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;

      cmplx8_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
          cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
          cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
          cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
          accrehihihi,accrelohihi,accrehilohi,accrelolohi,
          accrehihilo,accrelohilo,accrehilolo,accrelololo,
          accimhihihi,accimlohihi,accimhilohi,accimlolohi,
          accimhihilo,accimlohilo,accimhilolo,accimlololo,
          inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
          inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
          inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
          inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h,
          inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
          inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
          inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
          inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
          outputrehihihi_h,outputrelohihi_h,outputrehilohi_h,outputrelolohi_h,
          outputrehihilo_h,outputrelohilo_h,outputrehilolo_h,outputrelololo_h,
          outputimhihihi_h,outputimlohihi_h,outputimhilohi_h,outputimlolohi_h,
          outputimhihilo_h,outputimlohilo_h,outputimhilolo_h,outputimlololo_h,
          outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
          outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
          outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
          outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d,
          funvalrehihihi_h,funvalrelohihi_h,funvalrehilohi_h,funvalrelolohi_h,
          funvalrehihilo_h,funvalrelohilo_h,funvalrehilolo_h,funvalrelololo_h,
          funvalimhihihi_h,funvalimlohihi_h,funvalimhilohi_h,funvalimlolohi_h,
          funvalimhihilo_h,funvalimlohilo_h,funvalimhilolo_h,funvalimlololo_h,
          funvalrehihihi_d,funvalrelohihi_d,funvalrehilohi_d,funvalrelolohi_d,
          funvalrehihilo_d,funvalrelohilo_d,funvalrehilolo_d,funvalrelololo_d,
          funvalimhihihi_d,funvalimlohihi_d,funvalimhilohi_d,funvalimlolohi_d,
          funvalimhihilo_d,funvalimlohilo_d,funvalimhilolo_d,funvalimlololo_d,
          jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
          jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
          jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
          jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,
          jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
          jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
          jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
          jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
          rhsrehihihi_h,rhsrelohihi_h,rhsrehilohi_h,rhsrelolohi_h,
          rhsrehihilo_h,rhsrelohilo_h,rhsrehilolo_h,rhsrelololo_h,
          rhsimhihihi_h,rhsimlohihi_h,rhsimhilohi_h,rhsimlolohi_h,
          rhsimhihilo_h,rhsimlohilo_h,rhsimhilolo_h,rhsimlololo_h,
          rhsrehihihi_d,rhsrelohihi_d,rhsrehilohi_d,rhsrelolohi_d,
          rhsrehihilo_d,rhsrelohilo_d,rhsrehilolo_d,rhsrelololo_d,
          rhsimhihihi_d,rhsimlohihi_d,rhsimhilohi_d,rhsimlolohi_d,
          rhsimhihilo_d,rhsimlohilo_d,rhsimhilolo_d,rhsimlololo_d,
          urhsrehihihi_h,urhsrelohihi_h,urhsrehilohi_h,urhsrelolohi_h,
          urhsrehihilo_h,urhsrelohilo_h,urhsrehilolo_h,urhsrelololo_h,
          urhsimhihihi_h,urhsimlohihi_h,urhsimhilohi_h,urhsimlolohi_h,
          urhsimhihilo_h,urhsimlohilo_h,urhsimhilolo_h,urhsimlololo_h,
          urhsrehihihi_d,urhsrelohihi_d,urhsrehilohi_d,urhsrelolohi_d,
          urhsrehihilo_d,urhsrelohilo_d,urhsrehilolo_d,urhsrelololo_d,
          urhsimhihihi_d,urhsimlohihi_d,urhsimhilohi_d,urhsimlolohi_d,
          urhsimhihilo_d,urhsimlohilo_d,urhsimhilolo_d,urhsimlololo_d,
          solrehihihi_h,solrelohihi_h,solrehilohi_h,solrelolohi_h,
          solrehihilo_h,solrelohilo_h,solrehilolo_h,solrelololo_h,
          solimhihihi_h,solimlohihi_h,solimhilohi_h,solimlolohi_h,
          solimhihilo_h,solimlohilo_h,solimhilolo_h,solimlololo_h,
          solrehihihi_d,solrelohihi_d,solrehilohi_d,solrelolohi_d,
          solrehihilo_d,solrelohilo_d,solrehilolo_d,solrelololo_d,
          solimhihihi_d,solimlohihi_d,solimhilohi_d,solimlolohi_d,
          solimhihilo_d,solimlohilo_d,solimhilolo_d,solimlololo_d,
          Qrehihihi_h,Qrelohihi_h,Qrehilohi_h,Qrelolohi_h,
          Qrehihilo_h,Qrelohilo_h,Qrehilolo_h,Qrelololo_h,
          Qimhihihi_h,Qimlohihi_h,Qimhilohi_h,Qimlolohi_h,
          Qimhihilo_h,Qimlohilo_h,Qimhilolo_h,Qimlololo_h,
          Qrehihihi_d,Qrelohihi_d,Qrehilohi_d,Qrelolohi_d,
          Qrehihilo_d,Qrelohilo_d,Qrehilolo_d,Qrelololo_d,
          Qimhihihi_d,Qimlohihi_d,Qimhilohi_d,Qimlolohi_d,
          Qimhihilo_d,Qimlohilo_d,Qimhilolo_d,Qimlololo_d,
          Rrehihihi_h,Rrelohihi_h,Rrehilohi_h,Rrelolohi_h,
          Rrehihilo_h,Rrelohilo_h,Rrehilolo_h,Rrelololo_h,
          Rimhihihi_h,Rimlohihi_h,Rimhilohi_h,Rimlolohi_h,
          Rimhihilo_h,Rimlohilo_h,Rimhilolo_h,Rimlololo_h,
          Rrehihihi_d,Rrelohihi_d,Rrehilohi_d,Rrelolohi_d,
          Rrehihilo_d,Rrelohilo_d,Rrehilolo_d,Rrelololo_d,
          Rimhihihi_d,Rimlohihi_d,Rimhilohi_d,Rimlolohi_d,
          Rimhihilo_d,Rimlohilo_d,Rimhilolo_d,Rimlololo_d,
          workmatrehihihi,workmatrelohihi,workmatrehilohi,workmatrelolohi,
          workmatrehihilo,workmatrelohilo,workmatrehilolo,workmatrelololo,
          workmatimhihihi,workmatimlohihi,workmatimhilohi,workmatimlolohi,
          workmatimhihilo,workmatimlohilo,workmatimhilolo,workmatimlololo,
          workvecrehihihi,workvecrelohihi,workvecrehilohi,workvecrelolohi,
          workvecrehihilo,workvecrelohilo,workvecrehilolo,workvecrelololo,
          workvecimhihihi,workvecimlohihi,workvecimhilohi,workvecimlolohi,
          workvecimhihilo,workvecimlohilo,workvecimhilolo,workvecimlololo,
          resvecrehihihi,resvecrelohihi,resvecrehilohi,resvecrelolohi,
          resvecrehihilo,resvecrelohilo,resvecrehilolo,resvecrelololo,
          resvecimhihihi,resvecimlohihi,resvecimhilohi,resvecimlolohi,
          resvecimhihilo,resvecimlohilo,resvecimhilolo,resvecimlololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,vrblvl,mode);
   }
   return 0;
}
