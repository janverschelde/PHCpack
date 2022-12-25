// The file dbl2_monomial_systems.cpp defines the functions specified in
// the file dbl2_monomial_systems.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include "octo_double_functions.h"
#include "random8_vectors.h"
#include "random8_series.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomial_systems.h"

using namespace std;

void make_real8_exponentials
 ( int dim, int  deg,
   double **shihihi, double **slohihi, double **shilohi, double **slolohi,
   double **shihilo, double **slohilo, double **shilolo, double **slololo )
{
   double rndhihihi,rndlohihi,rndhilohi,rndlolohi;
   double rndhihilo,rndlohilo,rndhilolo,rndlololo;
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   const double fac = 64.0;      // 1024.0;
   const double inc = 63.0/64.0; // 1023.0/1024.0;

   for(int i=0; i<dim; i++)
   {
      // rnd is in [-1, +1]
      random_octo_double(&rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                         &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo); 
      odf_div(rndhihihi,rndlohihi,rndhilohi,rndlolohi,
              rndhihilo,rndlohilo,rndhilolo,rndlololo,
              fac,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,  // acc is in
              &acchihilo,&acclohilo,&acchilolo,&acclololo); // [-1/fac, +1/fac]
      
      if(rndhihihi < 0) // if -1/fac <= rnd       < 0
      {                 // then   -1 <= rnd - inc < -inc
         odf_sub(acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo,
                 inc,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                 &rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                 &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
      }
      else              // if    0  <= rnd       <= 1/fac
      {                 // then inc <= rnd + inc <= 1
         odf_add(acchihihi,acclohihi,acchilohi,acclolohi,
                 acchihilo,acclohilo,acchilolo,acclololo,
                 inc,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                 &rndhihihi,&rndlohihi,&rndhilohi,&rndlolohi,
                 &rndhihilo,&rndlohilo,&rndhilolo,&rndlololo);
      }
      dbl8_exponential
         (deg,rndhihihi,rndlohihi,rndhilohi,rndlolohi,
              rndhihilo,rndlohilo,rndhilolo,rndlololo,
              shihihi[i],slohihi[i],shilohi[i],slolohi[i],
              shihilo[i],slohilo[i],shilolo[i],slololo[i]);
   }
}

void make_complex8_exponentials
 ( int dim, int deg,
   double **srehihihi, double **srelohihi,
   double **srehilohi, double **srelolohi,
   double **srehihilo, double **srelohilo,
   double **srehilolo, double **srelololo,
   double **simhihihi, double **simlohihi,
   double **simhilohi, double **simlolohi,
   double **simhihilo, double **simlohilo,
   double **simhilolo, double **simlololo )
{
   double xrehihihi,xrelohihi,xrehilohi,xrelolohi;
   double xrehihilo,xrelohilo,xrehilolo,xrelololo;
   double ximhihihi,ximlohihi,ximhilohi,ximlolohi;
   double ximhihilo,ximlohilo,ximhilolo,ximlololo;
   double yhihihi,ylohihi,yhilohi,ylolohi;
   double yhihilo,ylohilo,yhilolo,ylololo;

   for(int i=0; i<dim; i++)
   {                                            // cosine of some angle
      random_octo_double
         (&xrehihihi,&xrelohihi,&xrehilohi,&xrelolohi,
          &xrehihilo,&xrelohilo,&xrehilolo,&xrelololo);
 
      odf_sqr(xrehihihi,xrelohihi,xrehilohi,xrelolohi,
              xrehihilo,xrelohilo,xrehilolo,xrelololo,
              &yhihihi,&ylohihi,&yhilohi,&ylolohi,
              &yhihilo,&ylohilo,&yhilolo,&ylololo);        // y = cos^2
      odf_minus(&yhihihi,&ylohihi,&yhilohi,&ylolohi,
                &yhihilo,&ylohilo,&yhilolo,&ylololo);      // y = -cos^2
      odf_inc_d(&yhihihi,&ylohihi,&yhilohi,&ylolohi,
                &yhihilo,&ylohilo,&yhilolo,&ylololo,1.0);  // y = 1 - cos^2
      odf_sqrt(yhihihi,ylohihi,yhilohi,ylolohi,
               yhihilo,ylohilo,yhilolo,ylololo,
               &ximhihihi,&ximlohihi,&ximhilohi,&ximlolohi,
               &ximhihilo,&ximlohilo,&ximhilolo,&ximlololo);
     // sin is sqrt(1-cos^2)

      cmplx8_exponential
         (deg,xrehihihi,xrelohihi,xrehilohi,xrelolohi,
              xrehihilo,xrelohilo,xrehilolo,xrelololo,
              ximhihihi,ximlohihi,ximhilohi,ximlolohi,
              ximhihilo,ximlohilo,ximhilolo,ximlololo,
          srehihihi[i],srelohihi[i],srehilohi[i],srelolohi[i],
          srehihilo[i],srelohilo[i],srehilolo[i],srelololo[i],
          simhihihi[i],simlohihi[i],simhilohi[i],simlolohi[i],
          simhihilo[i],simlohilo[i],simhilolo[i],simlololo[i]);
   }
}

void evaluate_real8_monomials
 ( int dim, int deg, int **rowsA,
   double **shihihi, double **slohihi, double **shilohi, double **slolohi,
   double **shihilo, double **slohilo, double **shilolo, double **slololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo )
{
   const int degp1 = deg+1;

   double *acchihihi = new double[degp1]; // accumulates product
   double *acclohihi = new double[degp1];
   double *acchilohi = new double[degp1];
   double *acclolohi = new double[degp1];
   double *acchihilo = new double[degp1];
   double *acclohilo = new double[degp1];
   double *acchilolo = new double[degp1];
   double *acclololo = new double[degp1];

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_dbl8_product
                  (deg,shihihi[j],slohihi[j],shilohi[j],slolohi[j],
                       shihilo[j],slohilo[j],shilolo[j],slololo[j],
                   rhshihihi[i],rhslohihi[i],rhshilohi[i],rhslolohi[i],
                   rhshihilo[i],rhslohilo[i],rhshilolo[i],rhslololo[i],
                   acchihihi,acclohihi,acchilohi,acclolohi,
                   acchihilo,acclohilo,acchilolo,acclololo);

               for(int L=0; L<degp1; L++)
               {
                  rhshihihi[i][L] = acchihihi[L];
                  rhslohihi[i][L] = acclohihi[L];
                  rhshilohi[i][L] = acchilohi[L];
                  rhslolohi[i][L] = acclolohi[L];
                  rhshihilo[i][L] = acchihilo[L];
                  rhslohilo[i][L] = acclohilo[L];
                  rhshilolo[i][L] = acchilolo[L];
                  rhslololo[i][L] = acclololo[L];
               }
            }
         }
      }
   }
   free(acchihihi); free(acclohihi); free(acchilohi); free(acclolohi);
   free(acchihilo); free(acclohilo); free(acchilolo); free(acclololo);
}

void evaluate_complex8_monomials
 ( int dim, int deg, int **rowsA,
   double **srehihihi, double **srelohihi,
   double **srehilohi, double **srelolohi,
   double **srehihilo, double **srelohilo,
   double **srehilolo, double **srelololo,
   double **simhihihi, double **simlohihi,
   double **simhilohi, double **simlolohi,
   double **simhihilo, double **simlohilo,
   double **simhilolo, double **simlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo )
{
   const int degp1 = deg+1;

   double *accrehihihi = new double[degp1]; // accumulates product
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

   for(int i=0; i<dim; i++)    // run over all monomials
   {
      for(int j=0; j<dim; j++) // run over all variables
      {
         if(rowsA[i][j] > 0)   // only multiply if positive exponent
         {
            for(int k=0; k<rowsA[i][j]; k++)
            {
               CPU_cmplx8_product
                  (deg,srehihihi[j],srelohihi[j],srehilohi[j],srelolohi[j],
                       srehihilo[j],srelohilo[j],srehilolo[j],srelololo[j],
                       simhihihi[j],simlohihi[j],simhilohi[j],simlolohi[j],
                       simhihilo[j],simlohilo[j],simhilolo[j],simlololo[j],
                   rhsrehihihi[i],rhsrelohihi[i],rhsrehilohi[i],rhsrelolohi[i],
                   rhsrehihilo[i],rhsrelohilo[i],rhsrehilolo[i],rhsrelololo[i],
                   rhsimhihihi[i],rhsimlohihi[i],rhsimhilohi[i],rhsimlolohi[i],
                   rhsimhihilo[i],rhsimlohilo[i],rhsimhilolo[i],rhsimlololo[i],
                   accrehihihi,accrelohihi,accrehilohi,accrelolohi,
                   accrehihilo,accrelohilo,accrehilolo,accrelololo,
                   accimhihihi,accimlohihi,accimhilohi,accimlolohi,
                   accimhihilo,accimlohilo,accimhilolo,accimlololo);

               for(int L=0; L<degp1; L++)
               {
                  rhsrehihihi[i][L] = accrehihihi[L];
                  rhsrelohihi[i][L] = accrelohihi[L];
                  rhsrehilohi[i][L] = accrehilohi[L];
                  rhsrelolohi[i][L] = accrelolohi[L];
                  rhsrehihilo[i][L] = accrehihilo[L];
                  rhsrelohilo[i][L] = accrelohilo[L];
                  rhsrehilolo[i][L] = accrehilolo[L];
                  rhsrelololo[i][L] = accrelololo[L];
                  rhsimhihihi[i][L] = accimhihihi[L];
                  rhsimlohihi[i][L] = accimlohihi[L];
                  rhsimhilohi[i][L] = accimhilohi[L];
                  rhsimlolohi[i][L] = accimlolohi[L];
                  rhsimhihilo[i][L] = accimhihilo[L];
                  rhsimlohilo[i][L] = accimlohilo[L];
                  rhsimhilolo[i][L] = accimhilolo[L];
                  rhsimlololo[i][L] = accimlololo[L];
               }
            }
         }
      }
   }
   free(accrehihihi); free(accrelohihi); free(accrehilohi); free(accrelolohi);
   free(accrehihilo); free(accrelohilo); free(accrehilolo); free(accrelololo);
   free(accimhihihi); free(accimlohihi); free(accimhilohi); free(accimlolohi);
   free(accimhihilo); free(accimlohilo); free(accimhilolo); free(accimlololo);
}

void make_real8_coefficients
 ( int nbrcol, int dim,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo )
{
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         random_octo_double
            (&cffhihihi[i][j][0],&cfflohihi[i][j][0],
             &cffhilohi[i][j][0],&cfflolohi[i][j][0],
             &cffhihilo[i][j][0],&cfflohilo[i][j][0],
             &cffhilolo[i][j][0],&cfflololo[i][j][0]);
}

void make_complex8_coefficients
 ( int nbrcol, int dim,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo )
{
   double xrehihihi,xrelohihi,ximhihihi,ximlohihi,yhihihi,ylohihi;
   double xrehilohi,xrelolohi,ximhilohi,ximlolohi,yhilohi,ylolohi;
   double xrehihilo,xrelohilo,ximhihilo,ximlohilo,yhihilo,ylohilo;
   double xrehilolo,xrelololo,ximhilolo,ximlololo,yhilolo,ylololo;

   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
      {
         random_octo_double
            (&xrehihihi,&xrelohihi,&xrehilohi,&xrelolohi,
             &xrehihilo,&xrelohilo,&xrehilolo,&xrelololo);
 
         odf_sqr(xrehihihi,xrelohihi,xrehilohi,xrelolohi,
                 xrehihilo,xrelohilo,xrehilolo,xrelololo,
                 &yhihihi,&ylohihi,&yhilohi,&ylolohi,
                 &yhihilo,&ylohilo,&yhilolo,&ylololo);        // y = cos^2
         odf_minus(&yhihihi,&ylohihi,&yhilohi,&ylolohi,
                   &yhihilo,&ylohilo,&yhilolo,&ylololo);      // y = -cos^2
         odf_inc_d(&yhihihi,&ylohihi,&yhilohi,&ylolohi,
                   &yhihilo,&ylohilo,&yhilolo,&ylololo,1.0);  // y = 1 - cos^2
         odf_sqrt(yhihihi,ylohihi,yhilohi,ylolohi,
                  yhihilo,ylohilo,yhilolo,ylololo,
                  &ximhihihi,&ximlohihi,&ximhilohi,&ximlolohi,
                  &ximhihilo,&ximlohilo,&ximhilolo,&ximlololo);
         // sin is sqrt(1-cos^2)

         cffrehihihi[i][j][0] = xrehihihi;
         cffrelohihi[i][j][0] = xrelohihi;
         cffrehilohi[i][j][0] = xrehilohi;
         cffrelolohi[i][j][0] = xrelolohi;
         cffrehihilo[i][j][0] = xrehihilo;
         cffrelohilo[i][j][0] = xrelohilo;
         cffrehilolo[i][j][0] = xrehilolo;
         cffrelololo[i][j][0] = xrelololo;
         cffimhihihi[i][j][0] = ximhihihi;
         cffimlohihi[i][j][0] = ximlohihi;
         cffimhilohi[i][j][0] = ximhilohi;
         cffimlolohi[i][j][0] = ximlolohi;
         cffimhihilo[i][j][0] = ximhihilo;
         cffimlohilo[i][j][0] = ximlohilo;
         cffimhilolo[i][j][0] = ximhilolo;
         cffimlololo[i][j][0] = ximlololo;
      }
}

void evaluate_real8_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo, int vrblvl )
{
   const int degp1 = deg+1;

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "evaluating at the series x ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "x[" << i << "][" << j << "] : "
                 << xhihihi[i][j] << "  " << xlohihi[i][j] << endl << "  "
                 << xhilohi[i][j] << "  " << xlolohi[i][j] << endl << "  "
                 << xhihilo[i][j] << "  " << xlohilo[i][j] << endl << "  "
                 << xhilolo[i][j] << "  " << xlololo[i][j] << endl;
   }
   double **prdrhshihihi = new double*[dim];
   double **prdrhslohihi = new double*[dim];
   double **prdrhshilohi = new double*[dim];
   double **prdrhslolohi = new double*[dim];
   double **prdrhshihilo = new double*[dim];
   double **prdrhslohilo = new double*[dim];
   double **prdrhshilolo = new double*[dim];
   double **prdrhslololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhshihihi[i] = new double[degp1];
      prdrhslohihi[i] = new double[degp1];
      prdrhshilohi[i] = new double[degp1];
      prdrhslolohi[i] = new double[degp1];
      prdrhshihilo[i] = new double[degp1];
      prdrhslohilo[i] = new double[degp1];
      prdrhshilolo[i] = new double[degp1];
      prdrhslololo[i] = new double[degp1];

      for(int k=0; k<degp1; k++) // initialize sum to 0
      {
         rhshihihi[i][k] = 0.0; rhslohihi[i][k] = 0.0;
         rhshilohi[i][k] = 0.0; rhslolohi[i][k] = 0.0;
         rhshihilo[i][k] = 0.0; rhslohilo[i][k] = 0.0;
         rhshilolo[i][k] = 0.0; rhslololo[i][k] = 0.0;
      }
   }
   rhshihihi[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++)   // initialize product to coefficient
      {
         prdrhshihihi[i][0] = cffhihihi[col][i][0];
         prdrhslohihi[i][0] = cfflohihi[col][i][0];
         prdrhshilohi[i][0] = cffhilohi[col][i][0];
         prdrhslolohi[i][0] = cfflolohi[col][i][0];
         prdrhshihilo[i][0] = cffhihilo[col][i][0];
         prdrhslohilo[i][0] = cfflohilo[col][i][0];
         prdrhshilolo[i][0] = cffhilolo[col][i][0];
         prdrhslololo[i][0] = cfflololo[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhshihihi[i][k] = 0.0; prdrhslohihi[i][k] = 0.0;
            prdrhshilohi[i][k] = 0.0; prdrhslolohi[i][k] = 0.0;
            prdrhshihilo[i][k] = 0.0; prdrhslohilo[i][k] = 0.0;
            prdrhshilolo[i][k] = 0.0; prdrhslololo[i][k] = 0.0;
         }
      }
      if(vrblvl > 1)
         cout << "Evaluating at column " << col << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         for(int k=0; k<dim; k++) rowsA[i][k] = 0;
         for(int k=0; k<nvr[col][i]; k++) rowsA[i][idx[col][i][k]] = 1;
         if(vrblvl > 1)
         {
            for(int k=0; k<dim; k++) cout << " " << rowsA[i][k];
            cout << endl;
         }
      }
      evaluate_real8_monomials
         (dim,deg,rowsA,
          xhihihi,xlohihi,xhilohi,xlolohi,xhihilo,xlohilo,xhilolo,xlololo,
          prdrhshihihi,prdrhslohihi,prdrhshilohi,prdrhslolohi,
          prdrhshihilo,prdrhslohilo,prdrhshilolo,prdrhslololo);

      if(vrblvl > 1)
      {
         cout << scientific << setprecision(16);
         for(int i=0; i<dim; i++)
         {
            cout << "value at dimension " << i << " :" << endl;
            for(int j=0; j<degp1; j++)
               cout << "prdrhs[" << i << "][" << j << "] : "
                    << prdrhshihihi[i][j] << "  "
                    << prdrhslohihi[i][j] << endl << "  "
                    << prdrhshilohi[i][j] << "  "
                    << prdrhslolohi[i][j] << endl << "  "
                    << prdrhshihilo[i][j] << "  "
                    << prdrhslohilo[i][j] << endl << "  "
                    << prdrhshilolo[i][j] << "  "
                    << prdrhslololo[i][j] << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               // rhs[i][k] += prdrhs[i][k];
               odf_inc(&rhshihihi[i][k],&rhslohihi[i][k],
                       &rhshilohi[i][k],&rhslolohi[i][k],
                       &rhshihilo[i][k],&rhslohilo[i][k],
                       &rhshilolo[i][k],&rhslololo[i][k],
                       prdrhshihihi[i][k],prdrhslohihi[i][k],
                       prdrhshilohi[i][k],prdrhslolohi[i][k],
                       prdrhshihilo[i][k],prdrhslohilo[i][k],
                       prdrhshilolo[i][k],prdrhslololo[i][k]);
            }
      }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "the evaluated series ..." << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << rhshihihi[i][j] << "  " << rhslohihi[i][j] << endl << "  "
                 << rhshilohi[i][j] << "  " << rhslolohi[i][j] << endl << "  "
                 << rhshihilo[i][j] << "  " << rhslohilo[i][j] << endl << "  "
                 << rhshilolo[i][j] << "  " << rhslololo[i][j] << endl;
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhshihihi[i]); free(prdrhslohihi[i]);
      free(prdrhshilohi[i]); free(prdrhslolohi[i]);
      free(prdrhshihilo[i]); free(prdrhslohilo[i]);
      free(prdrhshilolo[i]); free(prdrhslololo[i]);
   }
   free(prdrhshihihi); free(prdrhslohihi);
   free(prdrhshilohi); free(prdrhslolohi);
   free(prdrhshihilo); free(prdrhslohilo);
   free(prdrhshilolo); free(prdrhslololo);
}

void evaluate_complex8_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo,
   double **xrehihihi, double **xrelohihi,
   double **xrehilohi, double **xrelolohi,
   double **xrehihilo, double **xrelohilo,
   double **xrehilolo, double **xrelololo,
   double **ximhihihi, double **ximlohihi,
   double **ximhilohi, double **ximlolohi,
   double **ximhihilo, double **ximlohilo,
   double **ximhilolo, double **ximlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo, int vrblvl )
{
   const int degp1 = deg+1;

   double **prdrhsrehihihi = new double*[dim];
   double **prdrhsrelohihi = new double*[dim];
   double **prdrhsrehilohi = new double*[dim];
   double **prdrhsrelolohi = new double*[dim];
   double **prdrhsrehihilo = new double*[dim];
   double **prdrhsrelohilo = new double*[dim];
   double **prdrhsrehilolo = new double*[dim];
   double **prdrhsrelololo = new double*[dim];
   double **prdrhsimhihihi = new double*[dim];
   double **prdrhsimlohihi = new double*[dim];
   double **prdrhsimhilohi = new double*[dim];
   double **prdrhsimlolohi = new double*[dim];
   double **prdrhsimhihilo = new double*[dim];
   double **prdrhsimlohilo = new double*[dim];
   double **prdrhsimhilolo = new double*[dim];
   double **prdrhsimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      prdrhsrehihihi[i] = new double[degp1];
      prdrhsrelohihi[i] = new double[degp1];
      prdrhsrehilohi[i] = new double[degp1];
      prdrhsrelolohi[i] = new double[degp1];
      prdrhsrehihilo[i] = new double[degp1];
      prdrhsrelohilo[i] = new double[degp1];
      prdrhsrehilolo[i] = new double[degp1];
      prdrhsrelololo[i] = new double[degp1];
      prdrhsimhihihi[i] = new double[degp1];
      prdrhsimlohihi[i] = new double[degp1];
      prdrhsimhilohi[i] = new double[degp1];
      prdrhsimlolohi[i] = new double[degp1];
      prdrhsimhihilo[i] = new double[degp1];
      prdrhsimlohilo[i] = new double[degp1];
      prdrhsimhilolo[i] = new double[degp1];
      prdrhsimlololo[i] = new double[degp1];

      for(int k=0; k<degp1; k++)  // initialize sum to zero
      {
         rhsrehihihi[i][k] = 0.0; rhsimhihihi[i][k] = 0.0;
         rhsrelohihi[i][k] = 0.0; rhsimlohihi[i][k] = 0.0;
         rhsrehilohi[i][k] = 0.0; rhsimhilohi[i][k] = 0.0;
         rhsrelolohi[i][k] = 0.0; rhsimlolohi[i][k] = 0.0;
         rhsrehihilo[i][k] = 0.0; rhsimhihilo[i][k] = 0.0;
         rhsrelohilo[i][k] = 0.0; rhsimlohilo[i][k] = 0.0;
         rhsrehilolo[i][k] = 0.0; rhsimhilolo[i][k] = 0.0;
         rhsrelololo[i][k] = 0.0; rhsimlololo[i][k] = 0.0;
      }
   }
   rhsrehihihi[dim-1][0] = -1.0; // last coefficient of cyclic n-roots is -1

   for(int col=0; col<nbrcol; col++)
   {
      for(int i=0; i<dim; i++) // initialize product to coefficient
      {
         prdrhsrehihihi[i][0] = cffrehihihi[col][i][0];
         prdrhsrelohihi[i][0] = cffrelohihi[col][i][0];
         prdrhsrehilohi[i][0] = cffrehilohi[col][i][0];
         prdrhsrelolohi[i][0] = cffrelolohi[col][i][0];
         prdrhsrehihilo[i][0] = cffrehihilo[col][i][0];
         prdrhsrelohilo[i][0] = cffrelohilo[col][i][0];
         prdrhsrehilolo[i][0] = cffrehilolo[col][i][0];
         prdrhsrelololo[i][0] = cffrelololo[col][i][0];
         prdrhsimhihihi[i][0] = cffimhihihi[col][i][0];
         prdrhsimlohihi[i][0] = cffimlohihi[col][i][0];
         prdrhsimhilohi[i][0] = cffimhilohi[col][i][0];
         prdrhsimlolohi[i][0] = cffimlolohi[col][i][0];
         prdrhsimhihilo[i][0] = cffimhihilo[col][i][0];
         prdrhsimlohilo[i][0] = cffimlohilo[col][i][0];
         prdrhsimhilolo[i][0] = cffimhilolo[col][i][0];
         prdrhsimlololo[i][0] = cffimlololo[col][i][0];

         for(int k=1; k<degp1; k++)
         {
            prdrhsrehihihi[i][k] = 0.0; prdrhsimhihihi[i][k] = 0.0;
            prdrhsrelohihi[i][k] = 0.0; prdrhsimlohihi[i][k] = 0.0;
            prdrhsrehilohi[i][k] = 0.0; prdrhsimhilohi[i][k] = 0.0;
            prdrhsrelolohi[i][k] = 0.0; prdrhsimlolohi[i][k] = 0.0;
            prdrhsrehihilo[i][k] = 0.0; prdrhsimhihilo[i][k] = 0.0;
            prdrhsrelohilo[i][k] = 0.0; prdrhsimlohilo[i][k] = 0.0;
            prdrhsrehilolo[i][k] = 0.0; prdrhsimhilolo[i][k] = 0.0;
            prdrhsrelololo[i][k] = 0.0; prdrhsimlololo[i][k] = 0.0;
         }
      }
      if(vrblvl > 1)
         cout << "Evaluating at column " << col << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         for(int k=0; k<dim; k++) rowsA[i][k] = 0;
         for(int k=0; k<nvr[col][i]; k++) rowsA[i][idx[col][i][k]] = 1;
         if(vrblvl > 1)
         {
            for(int k=0; k<dim; k++) cout << " " << rowsA[i][k];
            cout << endl;
         }
      }
      evaluate_complex8_monomials
         (dim,deg,rowsA,
          xrehihihi,xrelohihi,xrehilohi,xrelolohi,
          xrehihilo,xrelohilo,xrehilolo,xrelololo,
          ximhihihi,ximlohihi,ximhilohi,ximlolohi,
          ximhihilo,ximlohilo,ximhilolo,ximlololo,
          prdrhsrehihihi,prdrhsrelohihi,prdrhsrehilohi,prdrhsrelolohi,
          prdrhsrehihilo,prdrhsrelohilo,prdrhsrehilolo,prdrhsrelololo,
          prdrhsimhihihi,prdrhsimlohihi,prdrhsimhilohi,prdrhsimlolohi,
          prdrhsimhihilo,prdrhsimlohilo,prdrhsimhilolo,prdrhsimlololo);

      for(int i=0; i<dim; i++)
      {
         int rowsum = 0;  // check on row sum is a patch ...
         for(int j=0; j<dim; j++) rowsum += rowsA[i][j];
         if(rowsum != 0)
            for(int k=0; k<degp1; k++)
            {
               // rhsre[i][k] += prdrhsre[i][k];
               odf_inc(&rhsrehihihi[i][k],&rhsrelohihi[i][k],
                       &rhsrehilohi[i][k],&rhsrelolohi[i][k],
                       &rhsrehihilo[i][k],&rhsrelohilo[i][k],
                       &rhsrehilolo[i][k],&rhsrelololo[i][k],
                       prdrhsrehihihi[i][k],prdrhsrelohihi[i][k],
                       prdrhsrehilohi[i][k],prdrhsrelolohi[i][k],
                       prdrhsrehihilo[i][k],prdrhsrelohilo[i][k],
                       prdrhsrehilolo[i][k],prdrhsrelololo[i][k]);
               // rhsim[i][k] += prdrhsim[i][k];
               odf_inc(&rhsimhihihi[i][k],&rhsimlohihi[i][k],
                       &rhsimhilohi[i][k],&rhsimlolohi[i][k],
                       &rhsimhihilo[i][k],&rhsimlohilo[i][k],
                       &rhsimhilolo[i][k],&rhsimlololo[i][k],
                       prdrhsimhihihi[i][k],prdrhsimlohihi[i][k],
                       prdrhsimhilohi[i][k],prdrhsimlolohi[i][k],
                       prdrhsimhihilo[i][k],prdrhsimlohilo[i][k],
                       prdrhsimhilolo[i][k],prdrhsimlololo[i][k]);
            }
      }
   }
   for(int i=0; i<dim; i++)
   {
      free(prdrhsrehihihi[i]); free(prdrhsimhihihi[i]);
      free(prdrhsrelohihi[i]); free(prdrhsimlohihi[i]);
      free(prdrhsrehilohi[i]); free(prdrhsimhilohi[i]);
      free(prdrhsrelolohi[i]); free(prdrhsimlolohi[i]);
      free(prdrhsrehihilo[i]); free(prdrhsimhihilo[i]);
      free(prdrhsrelohilo[i]); free(prdrhsimlohilo[i]);
      free(prdrhsrehilolo[i]); free(prdrhsimhilolo[i]);
      free(prdrhsrelololo[i]); free(prdrhsimlololo[i]);
   }
   free(prdrhsrehihihi); free(prdrhsimhihihi);
   free(prdrhsrelohihi); free(prdrhsimlohihi);
   free(prdrhsrehilohi); free(prdrhsimhilohi);
   free(prdrhsrelolohi); free(prdrhsimlolohi);
   free(prdrhsrehihilo); free(prdrhsimhihilo);
   free(prdrhsrelohilo); free(prdrhsimlohilo);
   free(prdrhsrehilolo); free(prdrhsimhilolo);
   free(prdrhsrelololo); free(prdrhsimlololo);
}
