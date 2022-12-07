// The file dbl2_newton_testers.cpp defines the functions with prototypes in
// the file dbl2_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_factorizations.h"
#include "dbl2_bals_host.h"
#include "dbl2_bals_kernels.h"
#include "dbl2_tail_kernels.h"
#include "dbl2_systems_host.h"
#include "dbl2_systems_kernels.h"
#include "dbl_bals_flopcounts.h"

using namespace std;

void real2_start_series_vector
 ( int dim, int deg, double *r0hi, double *r0lo,
   double **cffhi, double **cfflo )
{
   for(int i=0; i<dim; i++)
   {
      cffhi[i][0] = r0hi[i];
      cfflo[i][0] = r0lo[i];
      ddf_inc(&cffhi[i][0],&cfflo[i][0],1.0e-16,0.0);

      for(int j=1; j<=deg; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
      }
   }
}

void cmplx2_start_series_vector
 ( int dim, int deg,
   double *r0rehi, double *r0relo, double *r0imhi, double *r0imlo,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehi[i][0] = r0rehi[i];
      cffrelo[i][0] = r0relo[i];
      ddf_inc(&cffrehi[i][0],&cffrelo[i][0],1.0e-16,0.0);

      cffimhi[i][0] = r0imhi[i];
      cffimlo[i][0] = r0imlo[i];
      ddf_inc(&cffimhi[i][0],&cffimlo[i][0],1.0e-16,0.0);

      for(int j=1; j<=deg; j++)
      {
         cffrehi[i][j] = 0.0; cffrelo[i][j] = 0.0;
         cffimhi[i][j] = 0.0; cffimlo[i][j] = 0.0;
      }
   }
}

void dbl2_unit_series_vector
 ( int dim, int deg, double **cffhi, double **cfflo )
{
   for(int i=0; i<dim; i++)
   {
      cffhi[i][0] = 1.0;
      cfflo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhi[i][j] = 0.0;
         cfflo[i][j] = 0.0;
      }
   }
}

void cmplx2_unit_series_vector
 ( int dim, int deg,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehi[i][0] = 1.0; cffrelo[i][0] = 0.0;
      cffimhi[i][0] = 0.0; cffimlo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffrehi[i][j] = 0.0; cffrelo[i][j] = 0.0;
         cffimhi[i][j] = 0.0; cffimlo[i][j] = 0.0;
      }
   }
}

void dbl2_unit_series_vectors
 ( int nbr, int dim, int deg, double ***cffhi, double ***cfflo )
{
   for(int i=0; i<nbr; i++)
      dbl2_unit_series_vector(dim,deg,cffhi[i],cfflo[i]);
}

void cmplx2_unit_series_vectors
 ( int nbr, int dim, int deg,
   double ***cffrehi, double ***cffrelo,
   double ***cffimhi, double ***cffimlo )
{
   for(int i=0; i<nbr; i++)
      cmplx2_unit_series_vector
         (dim,deg,cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i]);
}

void dbl2_update_series
 ( int dim, int degp1, int startidx, double **xhi, double **xlo,
   double **dxhi, double **dxlo, int vrblvl )
{
   if(vrblvl > 1) cout << scientific << setprecision(16);

   if(vrblvl > 1)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xhi[i][j] << "  " << xlo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         ddf_inc(&xhi[i][j],&xlo[i][j],dxhi[j][i],dxlo[j][i]);
      }

   if(vrblvl > 1)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xhi[i][j] << "  " << xlo[i][j] << endl;
      }
   }
}

void cmplx2_update_series
 ( int dim, int degp1, int startidx,
   double **xrehi, double **xrelo, double **ximhi, double **ximlo,
   double **dxrehi, double **dxrelo, double **dximhi, double **dximlo,
   int vrblvl )
{
   if(vrblvl > 1) cout << scientific << setprecision(16);

   if(vrblvl > 1)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehi[i][j] << "  " << xrelo[i][j] << endl << "  "
                 << ximhi[i][j] << "  " << ximlo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=startidx; j<degp1; j++) 
      for(int i=0; i<dim; i++)
      {
         // xre[i][j] = xre[i][j] + dxre[j][i];
         ddf_inc(&xrehi[i][j],&xrelo[i][j],dxrehi[j][i],dxrelo[j][i]);
         // xim[i][j] = xim[i][j] + dxim[j][i];
         ddf_inc(&ximhi[i][j],&ximlo[i][j],dximhi[j][i],dximlo[j][i]);
      }

   if(vrblvl > 1)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehi[i][j] << "  " << xrelo[i][j] << endl << "  "
                 << ximhi[i][j] << "  " << ximlo[i][j] << endl;
      }
   }
}

double dbl2_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datahi_h, double ***datalo_h,
   double ***datahi_d, double ***datalo_d,
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
                    << datahi_h[k][i][j] << "  "
                    << datalo_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datahi_d[k][i][j] << "  "
                    << datalo_d[k][i][j] << endl;
            }
            errsum += abs(datahi_h[k][i][j] - datahi_d[k][i][j])
                    + abs(datalo_h[k][i][j] - datalo_d[k][i][j]);
         }

   return errsum;
}

double cmplx2_error3sum
 ( int dim1, int dim2, int dim3,
   double ***datarehi_h, double ***datarelo_h,
   double ***dataimhi_h, double ***dataimlo_h,
   double ***datarehi_d, double ***datarelo_d,
   double ***dataimhi_d, double ***dataimlo_d,
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
                    << datarehi_h[k][i][j] << "  "
                    << datarelo_h[k][i][j] << endl << "  "
                    << dataimhi_h[k][i][j] << "  "
                    << dataimlo_h[k][i][j] << endl;
               cout << banner << "_d[" // "output_d["
                    << k << "][" << i << "][" << j << "] : "
                    << datarehi_d[k][i][j] << "  "
                    << datarelo_d[k][i][j] << endl << "  "
                    << dataimhi_d[k][i][j] << "  "
                    << dataimlo_d[k][i][j] << endl;
            }
            errsum += abs(datarehi_h[k][i][j] - datarehi_d[k][i][j])
                    + abs(datarelo_h[k][i][j] - datarelo_d[k][i][j])
                    + abs(dataimhi_h[k][i][j] - dataimhi_d[k][i][j])
                    + abs(dataimlo_h[k][i][j] - dataimlo_d[k][i][j]);
         }

   return errsum;
}

double dbl2_error2sum
 ( int nrows, int ncols,
   double **datahi_h, double **datalo_h,
   double **datahi_d, double **datalo_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datahi_h[i][j] << "  "
                 << datalo_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datahi_d[i][j] << "  "
                 << datalo_d[i][j] << endl;
         }
         errsum += abs(datahi_h[i][j] - datahi_d[i][j])
                 + abs(datalo_h[i][j] - datalo_d[i][j]);
      }

   return errsum;
}

double cmplx2_error2sum
 ( int nrows, int ncols,
   double **datarehi_h, double **datarelo_h,
   double **dataimhi_h, double **dataimlo_h,
   double **datarehi_d, double **datarelo_d,
   double **dataimhi_d, double **dataimlo_d,
   std::string banner, int vrblvl )
{
   double errsum = 0.0;

   for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
      {
         if(vrblvl > 1)
         {
            cout << banner << "_h[" << i << "][" << j << "] : "
                 << datarehi_h[i][j] << "  "
                 << datarelo_h[i][j] << endl << "  "
                 << dataimhi_h[i][j] << "  "
                 << dataimlo_h[i][j] << endl;
            cout << banner << "_d[" << i << "][" << j << "] : "
                 << datarehi_d[i][j] << "  "
                 << datarelo_d[i][j] << endl << "  "
                 << dataimhi_d[i][j] << "  "
                 << dataimlo_d[i][j] << endl;
         }
         errsum += abs(datarehi_h[i][j] - datarehi_d[i][j])
                 + abs(datarelo_h[i][j] - datarelo_d[i][j])
                 + abs(dataimhi_h[i][j] - dataimhi_d[i][j])
                 + abs(dataimlo_h[i][j] - dataimlo_d[i][j]);
      }

   return errsum;
}
