// The file dbl8_systems_host.cpp defines the functions with prototypes in
// the file dbl8_systems_host.h.

#include <iostream>
#include <iomanip>
#include "octo_double_functions.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"

using namespace std;

void CPU_dbl8_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo, int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_dbl8_product
                     (deg,inputhihihi[idxvar],inputlohihi[idxvar],
                          inputhilohi[idxvar],inputlolohi[idxvar],
                          inputhihilo[idxvar],inputlohilo[idxvar],
                          inputhilolo[idxvar],inputlololo[idxvar],
                      cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
                      cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
                      acchihihi,   acclohihi,   acchilohi,   acclolohi,
                      acchihilo,   acclohilo,   acchilolo,   acclololo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffhihihi[i][L] = acchihihi[L];
                     cfflohihi[i][L] = acclohihi[L];
                     cffhilohi[i][L] = acchilohi[L];
                     cfflolohi[i][L] = acclolohi[L];
                     cffhihilo[i][L] = acchihilo[L];
                     cfflohilo[i][L] = acclohilo[L];
                     cffhilolo[i][L] = acchilolo[L];
                     cfflololo[i][L] = acclololo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << cffhihihi[i][j] << "  " << cfflohihi[i][j] << endl;
            cout << cffhilohi[i][j] << "  " << cfflolohi[i][j] << endl;
            cout << cffhihilo[i][j] << "  " << cfflohilo[i][j] << endl;
            cout << cffhilolo[i][j] << "  " << cfflololo[i][j] << endl;
         }
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << inputhihihi[i][j] << "  " << inputlohihi[i][j] << endl;
            cout << inputhilohi[i][j] << "  " << inputlolohi[i][j] << endl;
            cout << inputhihilo[i][j] << "  " << inputlohilo[i][j] << endl;
            cout << inputhilolo[i][j] << "  " << inputlololo[i][j] << endl;
         }
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl8_evaldiff
         (dim,nvr[i],deg,idx[i],
          cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
          cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
           inputhihihi,    inputlohihi,    inputhilohi,    inputlolohi,
           inputhihilo,    inputlohilo,    inputhilolo,    inputlololo,
          outputhihihi[i],outputlohihi[i],outputhilohi[i],outputlolohi[i],
          outputhihilo[i],outputlohilo[i],outputhilolo[i],outputlololo[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << outputhihihi[i][dim][j] << "  "
                 << outputlohihi[i][dim][j] << endl;
            cout << outputhilohi[i][dim][j] << "  "
                 << outputlolohi[i][dim][j] << endl;
            cout << outputhihilo[i][dim][j] << "  "
                 << outputlohilo[i][dim][j] << endl;
            cout << outputhilolo[i][dim][j] << "  "
                 << outputlololo[i][dim][j] << endl;
         }
      }
   }
   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];
               double tmphihihi,tmplohihi,tmphilohi,tmplolohi;
               double tmphihilo,tmplohilo,tmphilolo,tmplololo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // output[i][idxvar][k] = factor*output[i][idxvar][k];
                  odf_mul_od_d
                     (outputhihihi[i][idxvar][k],outputlohihi[i][idxvar][k],
                      outputhilohi[i][idxvar][k],outputlolohi[i][idxvar][k],
                      outputhihilo[i][idxvar][k],outputlohilo[i][idxvar][k],
                      outputhilolo[i][idxvar][k],outputlololo[i][idxvar][k],
                      factor,&tmphihihi,&tmplohihi,&tmphilohi,&tmplolohi,
                             &tmphihilo,&tmplohilo,&tmphilolo,&tmplololo);
                  outputhihihi[i][idxvar][k] = tmphihihi;
                  outputlohihi[i][idxvar][k] = tmplohihi;
                  outputhilohi[i][idxvar][k] = tmphilohi;
                  outputlolohi[i][idxvar][k] = tmplolohi;
                  outputhihilo[i][idxvar][k] = tmphihilo;
                  outputlohilo[i][idxvar][k] = tmplohilo;
                  outputhilolo[i][idxvar][k] = tmphilolo;
                  outputlololo[i][idxvar][k] = tmplololo;
               }
            }
         }
      }
   }
}

void CPU_cmplx8_evaluate_monomials
 ( int dim, int deg,
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
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi, 
   double **inputimhilohi, double **inputimlolohi, 
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo, 
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi, 
   double ***outputrehihilo, double ***outputrelohilo, 
   double ***outputrehilolo, double ***outputrelololo, 
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo, int vrblvl )
{
   for(int i=0; i<dim; i++) // common factors in the coefficients
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];

               for(int k=0; k<expfac[i][j]; k++)
               {
                  CPU_cmplx8_product
                     (deg,inputrehihihi[idxvar],inputrelohihi[idxvar],
                          inputrehilohi[idxvar],inputrelolohi[idxvar],
                          inputrehihilo[idxvar],inputrelohilo[idxvar],
                          inputrehilolo[idxvar],inputrelololo[idxvar],
                          inputimhihihi[idxvar],inputimlohihi[idxvar],
                          inputimhilohi[idxvar],inputimlolohi[idxvar],
                          inputimhihilo[idxvar],inputimlohilo[idxvar],
                          inputimhilolo[idxvar],inputimlololo[idxvar],
                      cffrehihihi[i],cffrelohihi[i],
                      cffrehilohi[i],cffrelolohi[i],
                      cffrehihilo[i],cffrelohilo[i],
                      cffrehilolo[i],cffrelololo[i],
                      cffimhihihi[i],cffimlohihi[i],
                      cffimhilohi[i],cffimlolohi[i],
                      cffimhihilo[i],cffimlohilo[i],
                      cffimhilolo[i],cffimlololo[i],
                      accrehihihi,accrelohihi,accrehilohi,accrelolohi,
                      accrehihilo,accrelohilo,accrehilolo,accrelololo,
                      accimhihihi,accimlohihi,accimhilohi,accimlolohi,
                      accimhihilo,accimlohilo,accimhilolo,accimlololo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffrehihihi[i][L] = accrehihihi[L];
                     cffrelohihi[i][L] = accrelohihi[L];
                     cffrehilohi[i][L] = accrehilohi[L];
                     cffrelolohi[i][L] = accrelolohi[L];
                     cffrehihilo[i][L] = accrehihilo[L];
                     cffrelohilo[i][L] = accrelohilo[L];
                     cffrehilolo[i][L] = accrehilolo[L];
                     cffrelololo[i][L] = accrelololo[L];
                     cffimhihihi[i][L] = accimhihihi[L];
                     cffimlohihi[i][L] = accimlohihi[L];
                     cffimhilohi[i][L] = accimhilohi[L];
                     cffimlolohi[i][L] = accimlolohi[L];
                     cffimhihilo[i][L] = accimhihilo[L];
                     cffimlohilo[i][L] = accimlohilo[L];
                     cffimhilolo[i][L] = accimhilolo[L];
                     cffimlololo[i][L] = accimlololo[L];
                  }
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "coefficients for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << cffrehihihi[i][j] << "  " << cffrelohihi[i][j] << endl
                 << "  "
                 << cffrehilohi[i][j] << "  " << cffrelolohi[i][j] << endl
                 << "  "
                 << cffrehihilo[i][j] << "  " << cffrelohilo[i][j] << endl
                 << "  "
                 << cffrehilolo[i][j] << "  " << cffrelololo[i][j] << endl
                 << "  "
                 << cffimhihihi[i][j] << "  " << cffimlohihi[i][j] << endl
                 << "  "
                 << cffimhilohi[i][j] << "  " << cffimlolohi[i][j] << endl
                 << "  "
                 << cffimhihilo[i][j] << "  " << cffimlohilo[i][j] << endl
                 << "  "
                 << cffimhilolo[i][j] << "  " << cffimlololo[i][j] << endl;
      }
      cout << "dim : " << dim << "  nvr :";
      for(int i=0; i<dim; i++) cout << " " << nvr[i];
      cout << endl;
      cout << "deg : " << deg;
      for(int i=0; i<dim; i++) 
      {
         cout << "  idx[" << i << "] :";
         for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j];
      }
      cout << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehihihi[i][j] << "  " << inputrelohihi[i][j] << endl
                 << "  "
                 << inputrehilohi[i][j] << "  " << inputrelolohi[i][j] << endl
                 << "  "
                 << inputrehihilo[i][j] << "  " << inputrelohilo[i][j] << endl
                 << "  "
                 << inputrehilolo[i][j] << "  " << inputrelololo[i][j] << endl
                 << "  "
                 << inputimhihihi[i][j] << "  " << inputimlohihi[i][j] << endl
                 << "  "
                 << inputimhilohi[i][j] << "  " << inputimlolohi[i][j] << endl
                 << "  "
                 << inputimhihilo[i][j] << "  " << inputimlohilo[i][j] << endl
                 << "  "
                 << inputimhilolo[i][j] << "  " << inputimlololo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_cmplx8_evaldiff
         (dim,nvr[i],deg,idx[i],
          cffrehihihi[i],cffrelohihi[i],cffrehilohi[i],cffrelolohi[i],
          cffrehihilo[i],cffrelohilo[i],cffrehilolo[i],cffrelololo[i],
          cffimhihihi[i],cffimlohihi[i],cffimhilohi[i],cffimlolohi[i],
          cffimhihilo[i],cffimlohilo[i],cffimhilolo[i],cffimlololo[i],
          inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
          inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
          inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
          inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
          outputrehihihi[i],outputrelohihi[i],
          outputrehilohi[i],outputrelolohi[i],
          outputrehihilo[i],outputrelohilo[i],
          outputrehilolo[i],outputrelololo[i],
          outputimhihihi[i],outputimlohihi[i],
          outputimhilohi[i],outputimlolohi[i],
          outputimhihilo[i],outputimlohilo[i],
          outputimhilolo[i],outputimlololo[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehihihi[i][dim][j] << "  "
                 << outputrelohihi[i][dim][j] << endl << "  "
                 << outputrehilohi[i][dim][j] << "  "
                 << outputrelolohi[i][dim][j] << endl << "  "
                 << outputrehihilo[i][dim][j] << "  "
                 << outputrelohilo[i][dim][j] << endl << "  "
                 << outputrehilolo[i][dim][j] << "  "
                 << outputrelololo[i][dim][j] << endl << "  "
                 << outputimhihihi[i][dim][j] << "  "
                 << outputimlohihi[i][dim][j] << endl << "  "
                 << outputimhilohi[i][dim][j] << "  "
                 << outputimlolohi[i][dim][j] << endl << "  "
                 << outputimhihilo[i][dim][j] << "  "
                 << outputimlohilo[i][dim][j] << endl << "  "
                 << outputimhilolo[i][dim][j] << "  "
                 << outputimlololo[i][dim][j] << endl;
      }
   }
   for(int i=0; i<dim; i++) // multiply derivatives with the powers
   {
      if(nbrfac[i] > 0) // there are common factors in monomial i
      {
         for(int j=0; j<nvr[i]; j++) // run over all exponents
         {
            if(expfac[i][j] > 0) // the j-th exponent with variable idx[i][j]
            {
               int idxvar = idx[i][j];
               double factor = (double) exp[i][j];
               double acchihihi,acclohihi,acchilohi,acclolohi;
               double acchihilo,acclohilo,acchilolo,acclololo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // outputre[i][idxvar][k] = factor*outputre[i][idxvar][k];
                  odf_mul_od_d
                     (outputrehihihi[i][idxvar][k],
                      outputrelohihi[i][idxvar][k],
                      outputrehilohi[i][idxvar][k],
                      outputrelolohi[i][idxvar][k],
                      outputrehihilo[i][idxvar][k],
                      outputrelohilo[i][idxvar][k],
                      outputrehilolo[i][idxvar][k],
                      outputrelololo[i][idxvar][k],factor,
                      &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                      &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  outputrehihihi[i][idxvar][k] = acchihihi;
                  outputrelohihi[i][idxvar][k] = acclohihi;
                  outputrehilohi[i][idxvar][k] = acchilohi;
                  outputrelolohi[i][idxvar][k] = acclolohi;
                  outputrehihilo[i][idxvar][k] = acchihilo;
                  outputrelohilo[i][idxvar][k] = acclohilo;
                  outputrehilolo[i][idxvar][k] = acchilolo;
                  outputrelololo[i][idxvar][k] = acclololo;
                  // outputim[i][idxvar][k] = factor*outputim[i][idxvar][k];
                  odf_mul_od_d
                     (outputimhihihi[i][idxvar][k],
                      outputimlohihi[i][idxvar][k],
                      outputimhilohi[i][idxvar][k],
                      outputimlolohi[i][idxvar][k],
                      outputimhihilo[i][idxvar][k],
                      outputimlohilo[i][idxvar][k],
                      outputimhilolo[i][idxvar][k],
                      outputimlololo[i][idxvar][k],factor,
                      &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                      &acchihilo,&acclohilo,&acchilolo,&acclololo);
                  outputimhihihi[i][idxvar][k] = acchihihi;
                  outputimlohihi[i][idxvar][k] = acclohihi;
                  outputimhilohi[i][idxvar][k] = acchilohi;
                  outputimlolohi[i][idxvar][k] = acclolohi;
                  outputimhihilo[i][idxvar][k] = acchihilo;
                  outputimlohilo[i][idxvar][k] = acclohilo;
                  outputimhilolo[i][idxvar][k] = acchilolo;
                  outputimlololo[i][idxvar][k] = acclololo;
               }
            }
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "after multiplication with the factors ..." << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehihihi[i][dim][j] << "  "
                 << outputrelohihi[i][dim][j] << endl << "  "
                 << outputrehilohi[i][dim][j] << "  "
                 << outputrelolohi[i][dim][j] << endl << "  "
                 << outputrehihilo[i][dim][j] << "  "
                 << outputrelohilo[i][dim][j] << endl << "  "
                 << outputrehilolo[i][dim][j] << "  "
                 << outputrelololo[i][dim][j] << endl << "  "
                 << outputimhihihi[i][dim][j] << "  "
                 << outputimlohihi[i][dim][j] << endl << "  "
                 << outputimhilohi[i][dim][j] << "  "
                 << outputimlolohi[i][dim][j] << endl << "  "
                 << outputimhihilo[i][dim][j] << "  "
                 << outputimlohilo[i][dim][j] << endl << "  "
                 << outputimhilolo[i][dim][j] << "  "
                 << outputimlololo[i][dim][j] << endl;
      }
   }
}

void CPU_dbl8_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **acchihihi, double **acclohihi,
   double **acchilohi, double **acclolohi,
   double **acchihilo, double **acclohilo,
   double **acchilolo, double **acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhihihi[i][j] = 0.0; funvallohihi[i][j] = 0.0;
         funvalhilohi[i][j] = 0.0; funvallolohi[i][j] = 0.0;
         funvalhihilo[i][j] = 0.0; funvallohilo[i][j] = 0.0;
         funvalhilolo[i][j] = 0.0; funvallololo[i][j] = 0.0;
      }
   funvalhihihi[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalhihihi[k][i][j] = 0.0; jacvallohihi[k][i][j] = 0.0;
            jacvalhilohi[k][i][j] = 0.0; jacvallolohi[k][i][j] = 0.0;
            jacvalhihilo[k][i][j] = 0.0; jacvallohilo[k][i][j] = 0.0;
            jacvalhilolo[k][i][j] = 0.0; jacvallololo[k][i][j] = 0.0;
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffhihihi[i][j][k] << "  "
                    << cfflohihi[i][j][k] << endl
                    << "  "
                    << cffhilohi[i][j][k] << "  "
                    << cfflolohi[i][j][k] << endl
                    << "  "
                    << cffhihilo[i][j][k] << "  "
                    << cfflohilo[i][j][k] << endl
                    << "  "
                    << cffhilolo[i][j][k] << "  "
                    << cfflololo[i][j][k] << endl;
         }
      }
      cout << "dim : " << dim << endl;
      for(int i=0; i<nbrcol; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";
            cout << "  idx[" << i << "][" << j << "] :";
            for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
            cout << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputhihihi[i][j] << "  " << inputlohihi[i][j] << endl
                 << "  "
                 << inputhilohi[i][j] << "  " << inputlolohi[i][j] << endl
                 << "  "
                 << inputhihilo[i][j] << "  " << inputlohilo[i][j] << endl
                 << "  "
                 << inputhilolo[i][j] << "  " << inputlololo[i][j] << endl;
      }
   }
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // evaluate j-th monomial in column i
         {
            for(int k=0; k<=dim; k++)
               for(int L=0; L<degp1; L++)
               {
                  acchihihi[k][L] = 0.0; acclohihi[k][L] = 0.0;
                  acchilohi[k][L] = 0.0; acclolohi[k][L] = 0.0;
                  acchihilo[k][L] = 0.0; acclohilo[k][L] = 0.0;
                  acchilolo[k][L] = 0.0; acclololo[k][L] = 0.0;
               }

            CPU_dbl8_evaldiff
               (dim,nvr[i][j],deg,idx[i][j],
                cffhihihi[i][j],cfflohihi[i][j],
                cffhilohi[i][j],cfflolohi[i][j],
                cffhihilo[i][j],cfflohilo[i][j],
                cffhilolo[i][j],cfflololo[i][j],
                inputhihihi,inputlohihi,inputhilohi,inputlolohi,
                inputhihilo,inputlohilo,inputhilolo,inputlololo,
                acchihihi,acclohihi,acchilohi,acclolohi,
                acchihilo,acclohilo,acchilolo,acclololo);

            for(int L=0; L<degp1; L++) // funval[j][L] += acc[dim][L];
            {
               odf_inc(&funvalhihihi[j][L],&funvallohihi[j][L],
                       &funvalhilohi[j][L],&funvallolohi[j][L],
                       &funvalhihilo[j][L],&funvallohilo[j][L],
                       &funvalhilolo[j][L],&funvallololo[j][L],
                       acchihihi[dim][L],acclohihi[dim][L],
                       acchilohi[dim][L],acclolohi[dim][L],
                       acchihilo[dim][L],acclohilo[dim][L],
                       acchilolo[dim][L],acclololo[dim][L]);
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  // jacval[L][j][idxval] += acc[idxval][L];
                  odf_inc(&jacvalhihihi[L][j][idxval],
                          &jacvallohihi[L][j][idxval],
                          &jacvalhilohi[L][j][idxval],
                          &jacvallolohi[L][j][idxval],
                          &jacvalhihilo[L][j][idxval],
                          &jacvallohilo[L][j][idxval],
                          &jacvalhilolo[L][j][idxval],
                          &jacvallololo[L][j][idxval],
                          acchihihi[idxval][L],acclohihi[idxval][L],
                          acchilohi[idxval][L],acclolohi[idxval][L],
                          acchihilo[idxval][L],acclohilo[idxval][L],
                          acchilolo[idxval][L],acclololo[idxval][L]);
               }
            }
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<degp1; j++)
            cout << funvalhihihi[i][j] << "  " << funvallohihi[i][j] << endl
                 << "  "
                 << funvalhilohi[i][j] << "  " << funvallolohi[i][j] << endl
                 << "  "
                 << funvalhihilo[i][j] << "  " << funvallohilo[i][j] << endl
                 << "  "
                 << funvalhilolo[i][j] << "  " << funvallololo[i][j] << endl;
      }
   }
}

void CPU_cmplx8_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffrehihihi, double ***cffrelohihi,
   double ***cffrehilohi, double ***cffrelolohi,
   double ***cffrehihilo, double ***cffrelohilo,
   double ***cffrehilolo, double ***cffrelololo,
   double ***cffimhihihi, double ***cffimlohihi,
   double ***cffimhilohi, double ***cffimlolohi,
   double ***cffimhihilo, double ***cffimlohilo,
   double ***cffimhilolo, double ***cffimlololo,
   double **accrehihihi, double **accrelohihi,
   double **accrehilohi, double **accrelolohi,
   double **accrehihilo, double **accrelohilo,
   double **accrehilolo, double **accrelololo,
   double **accimhihihi, double **accimlohihi,
   double **accimhilohi, double **accimlolohi,
   double **accimhihilo, double **accimlohilo,
   double **accimhilolo, double **accimlololo,
   double **inputrehihihi, double **inputrelohihi,
   double **inputrehilohi, double **inputrelolohi,
   double **inputrehihilo, double **inputrelohilo,
   double **inputrehilolo, double **inputrelololo,
   double **inputimhihihi, double **inputimlohihi,
   double **inputimhilohi, double **inputimlolohi,
   double **inputimhihilo, double **inputimlohilo,
   double **inputimhilolo, double **inputimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehihihi[i][j] = 0.0; funvalrelohihi[i][j] = 0.0;
         funvalrehilohi[i][j] = 0.0; funvalrelolohi[i][j] = 0.0;
         funvalrehihilo[i][j] = 0.0; funvalrelohilo[i][j] = 0.0;
         funvalrehilolo[i][j] = 0.0; funvalrelololo[i][j] = 0.0;
         funvalimhihihi[i][j] = 0.0; funvalimlohihi[i][j] = 0.0;
         funvalimhilohi[i][j] = 0.0; funvalimlolohi[i][j] = 0.0;
         funvalimhihilo[i][j] = 0.0; funvalimlohilo[i][j] = 0.0;
         funvalimhilolo[i][j] = 0.0; funvalimlololo[i][j] = 0.0;
      }
   funvalrehihihi[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalrehihihi[k][i][j] = 0.0; jacvalrelohihi[k][i][j] = 0.0;
            jacvalrehilohi[k][i][j] = 0.0; jacvalrelolohi[k][i][j] = 0.0;
            jacvalrehihilo[k][i][j] = 0.0; jacvalrelohilo[k][i][j] = 0.0;
            jacvalrehilolo[k][i][j] = 0.0; jacvalrelololo[k][i][j] = 0.0;
            jacvalimhihihi[k][i][j] = 0.0; jacvalimlohihi[k][i][j] = 0.0;
            jacvalimhilohi[k][i][j] = 0.0; jacvalimlolohi[k][i][j] = 0.0;
            jacvalimhihilo[k][i][j] = 0.0; jacvalimlohilo[k][i][j] = 0.0;
            jacvalimhilolo[k][i][j] = 0.0; jacvalimlololo[k][i][j] = 0.0;
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++)
               cout << cffrehihihi[i][j][k] << "  "
                    << cffrelohihi[i][j][k] << endl
                    << "  "
                    << cffrehilohi[i][j][k] << "  "
                    << cffrelolohi[i][j][k] << endl
                    << "  "
                    << cffrehihilo[i][j][k] << "  "
                    << cffrelohilo[i][j][k] << endl
                    << "  "
                    << cffrehilolo[i][j][k] << "  "
                    << cffrelololo[i][j][k] << endl
                    << "  "
                    << cffimhihihi[i][j][k] << "  "
                    << cffimlohihi[i][j][k] << endl
                    << "  "
                    << cffimhilohi[i][j][k] << "  "
                    << cffimlolohi[i][j][k] << endl
                    << "  "
                    << cffimhihilo[i][j][k] << "  "
                    << cffimlohilo[i][j][k] << endl
                    << "  "
                    << cffimhilolo[i][j][k] << "  "
                    << cffimlololo[i][j][k] << endl;
         }
      }
      cout << "dim : " << dim << endl;
      for(int i=0; i<nbrcol; i++)
      {
         for(int j=0; j<dim; j++)
         {
            cout << "nvr[" << i << "][" << j << "] : " << nvr[i][j] << " :";
            cout << "  idx[" << i << "][" << j << "] :";
            for(int k=0; k<nvr[i][j]; k++) cout << " " << idx[i][j][k];
            cout << endl;
         }
      }
      for(int i=0; i<dim; i++)
      {
         cout << "input series for variable " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << inputrehihihi[i][j] << "  " << inputrelohihi[i][j] << endl
                 << "  "
                 << inputrehilohi[i][j] << "  " << inputrelolohi[i][j] << endl
                 << "  "
                 << inputrehihilo[i][j] << "  " << inputrelohilo[i][j] << endl
                 << "  "
                 << inputrehilolo[i][j] << "  " << inputrelololo[i][j] << endl
                 << "  "
                 << inputimhihihi[i][j] << "  " << inputimlohihi[i][j] << endl
                 << "  "
                 << inputimhilohi[i][j] << "  " << inputimlolohi[i][j] << endl
                 << "  "
                 << inputimhihilo[i][j] << "  " << inputimlohilo[i][j] << endl
                 << "  "
                 << inputimhilolo[i][j] << "  " << inputimlololo[i][j] << endl;
      }
   }
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // evaluate j-th monomial in column i
         {
            CPU_cmplx8_evaldiff
               (dim,nvr[i][j],deg,idx[i][j],
                cffrehihihi[i][j],cffrelohihi[i][j],
                cffrehilohi[i][j],cffrelolohi[i][j],
                cffrehihilo[i][j],cffrelohilo[i][j],
                cffrehilolo[i][j],cffrelololo[i][j],
                cffimhihihi[i][j],cffimlohihi[i][j],
                cffimhilohi[i][j],cffimlolohi[i][j],
                cffimhihilo[i][j],cffimlohilo[i][j],
                cffimhilolo[i][j],cffimlololo[i][j],
                inputrehihihi,inputrelohihi,inputrehilohi,inputrelolohi,
                inputrehihilo,inputrelohilo,inputrehilolo,inputrelololo,
                inputimhihihi,inputimlohihi,inputimhilohi,inputimlolohi,
                inputimhihilo,inputimlohilo,inputimhilolo,inputimlololo,
                accrehihihi,accrelohihi,accrehilohi,accrelolohi,
                accrehihilo,accrelohilo,accrehilolo,accrelololo,
                accimhihihi,accimlohihi,accimhilohi,accimlolohi,
                accimhihilo,accimlohilo,accimhilolo,accimlololo);

            for(int L=0; L<degp1; L++)
            {
               // funvalre[j][L] += accre[dim][L];
               odf_inc(&funvalrehihihi[j][L],&funvalrelohihi[j][L],
                       &funvalrehilohi[j][L],&funvalrelolohi[j][L],
                       &funvalrehihilo[j][L],&funvalrelohilo[j][L],
                       &funvalrehilolo[j][L],&funvalrelololo[j][L],
                       accrehihihi[dim][L],accrelohihi[dim][L],
                       accrehilohi[dim][L],accrelolohi[dim][L],
                       accrehihilo[dim][L],accrelohilo[dim][L],
                       accrehilolo[dim][L],accrelololo[dim][L]);
               // funvalim[j][L] += accim[dim][L];
               odf_inc(&funvalimhihihi[j][L],&funvalimlohihi[j][L],
                       &funvalimhilohi[j][L],&funvalimlolohi[j][L],
                       &funvalimhihilo[j][L],&funvalimlohilo[j][L],
                       &funvalimhilolo[j][L],&funvalimlololo[j][L],
                       accimhihihi[dim][L],accimlohihi[dim][L],
                       accimhilohi[dim][L],accimlolohi[dim][L],
                       accimhihilo[dim][L],accimlohilo[dim][L],
                       accimhilolo[dim][L],accimlololo[dim][L]);
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  // jacvalre[L][j][idxval] += accre[idxval][L];
                  odf_inc(&jacvalrehihihi[L][j][idxval],
                          &jacvalrelohihi[L][j][idxval],
                          &jacvalrehilohi[L][j][idxval],
                          &jacvalrelolohi[L][j][idxval],
                          &jacvalrehihilo[L][j][idxval],
                          &jacvalrelohilo[L][j][idxval],
                          &jacvalrehilolo[L][j][idxval],
                          &jacvalrelololo[L][j][idxval],
                          accrehihihi[idxval][L],accrelohihi[idxval][L],
                          accrehilohi[idxval][L],accrelolohi[idxval][L],
                          accrehihilo[idxval][L],accrelohilo[idxval][L],
                          accrehilolo[idxval][L],accrelololo[idxval][L]);
                  // jacvalim[L][j][idxval] += accim[idxval][L];
                  odf_inc(&jacvalimhihihi[L][j][idxval],
                          &jacvalimlohihi[L][j][idxval],
                          &jacvalimhilohi[L][j][idxval],
                          &jacvalimlolohi[L][j][idxval],
                          &jacvalimhihilo[L][j][idxval],
                          &jacvalimlohilo[L][j][idxval],
                          &jacvalimhilolo[L][j][idxval],
                          &jacvalimlololo[L][j][idxval],
                          accimhihihi[idxval][L],accimlohihi[idxval][L],
                          accimhilohi[idxval][L],accimlolohi[idxval][L],
                          accimhihilo[idxval][L],accimlohilo[idxval][L],
                          accimhilolo[idxval][L],accimlololo[idxval][L]);
               }
            }
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << funvalrehihihi[i][j] << "  "
                 << funvalrelohihi[i][j] << endl << "  "
                 << funvalrehilohi[i][j] << "  "
                 << funvalrelolohi[i][j] << endl << "  "
                 << funvalrehihilo[i][j] << "  "
                 << funvalrelohilo[i][j] << endl << "  "
                 << funvalrehilolo[i][j] << "  "
                 << funvalrelololo[i][j] << endl << "  "
                 << funvalimhihihi[i][j] << "  "
                 << funvalimlohihi[i][j] << endl << "  "
                 << funvalimhilohi[i][j] << "  "
                 << funvalimlolohi[i][j] << endl << "  "
                 << funvalimhihilo[i][j] << "  "
                 << funvalimlohilo[i][j] << endl << "  "
                 << funvalimhilolo[i][j] << "  "
                 << funvalimlololo[i][j] << endl;
      }
   }
}

void dbl8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbhihihi, double **mblohihi, double **mbhilohi, double **mblolohi,
   double **mbhihilo, double **mblohilo, double **mbhilolo, double **mblololo,
   double damper,
   double ***outputhihihi, double ***outputlohihi,
   double ***outputhilohi, double ***outputlolohi,
   double ***outputhihilo, double ***outputlohilo,
   double ***outputhilolo, double ***outputlololo,
   double **funvalhihihi, double **funvallohihi, 
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo, 
   double **funvalhilolo, double **funvallololo, 
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double ***jacvalhihihi, double ***jacvallohihi,
   double ***jacvalhilohi, double ***jacvallolohi,
   double ***jacvalhihilo, double ***jacvallohilo,
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhihihi[i][j] = outputhihihi[i][dim][j];
         funvallohihi[i][j] = outputlohihi[i][dim][j];
         funvalhilohi[i][j] = outputhilohi[i][dim][j];
         funvallolohi[i][j] = outputlolohi[i][dim][j];
         funvalhihilo[i][j] = outputhihilo[i][dim][j];
         funvallohilo[i][j] = outputlohilo[i][dim][j];
         funvalhilolo[i][j] = outputhilolo[i][dim][j];
         funvallololo[i][j] = outputlololo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << funvalhihihi[i][0] << "  "
                            << funvallohihi[i][0] << endl;
         cout << "     " << funvalhilohi[i][0] << "  "
                         << funvallolohi[i][0] << endl;
         cout << "     " << funvalhihilo[i][0] << "  "
                         << funvallohilo[i][0] << endl;
         cout << "     " << funvalhilolo[i][0] << "  "
                         << funvallololo[i][0] << endl;
      }
   }
   /*
    * Linearize the function values in the rhs and swap sign,
    * but keep in mind that the right hand side is mbre + I*mbim,
    * as the monomial system is x^U = rhs, or x^U - rhs = 0,
    * so from the funval we substract rhs and then flip sign
    * in the computation of the rhs of the linear system.
    */
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         // rhs[i][j] = -(funval[j][i] - mb[j][i]);
         odf_sub(funvalhihihi[j][i],funvallohihi[j][i],
                 funvalhilohi[j][i],funvallolohi[j][i],
                 funvalhihilo[j][i],funvallohilo[j][i],
                 funvalhilolo[j][i],funvallololo[j][i],
                 mbhihihi[j][i],mblohihi[j][i],mbhilohi[j][i],mblolohi[j][i],
                 mbhihilo[j][i],mblohilo[j][i],mbhilolo[j][i],mblololo[j][i],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhshihihi[i][j] = -acchihihi;
         rhslohihi[i][j] = -acclohihi;
         rhshilohi[i][j] = -acchilohi;
         rhslolohi[i][j] = -acclolohi;
         rhshihilo[i][j] = -acchihilo;
         rhslohilo[i][j] = -acclohilo;
         rhshilolo[i][j] = -acchilolo;
         rhslololo[i][j] = -acclololo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << rhshihihi[i][j] << "  " << rhslohihi[i][j] << endl;
            cout << rhshilohi[i][j] << "  " << rhslolohi[i][j] << endl;
            cout << rhshihilo[i][j] << "  " << rhslohilo[i][j] << endl;
            cout << rhshilolo[i][j] << "  " << rhslololo[i][j] << endl;
         }
      }
   }
   // output[i][idx[i][k]] for k from 0 to nvr[i] has the coefficients
   // of the series of the derivative of the i-th monomial with respect
   // to the variable idx[i][k].
   for(int i=0; i<dim; i++)          // the i-th monomial
   {
      int *indexes = idx[i];         // indices of the variables
      for(int k=0; k<nvr[i]; k++)    // derivative w.r.t. idx[i][k]
      {                              // has j-th coefficient
         int idxval = indexes[k];

         for(int j=0; j<degp1; j++) 
         {
            jacvalhihihi[j][i][idxval] = outputhihihi[i][idxval][j];
            jacvallohihi[j][i][idxval] = outputlohihi[i][idxval][j];
            jacvalhilohi[j][i][idxval] = outputhilohi[i][idxval][j];
            jacvallolohi[j][i][idxval] = outputlolohi[i][idxval][j];
            jacvalhihilo[j][i][idxval] = outputhihilo[i][idxval][j];
            jacvallohilo[j][i][idxval] = outputlohilo[i][idxval][j];
            jacvalhilolo[j][i][idxval] = outputhilolo[i][idxval][j];
            jacvallololo[j][i][idxval] = outputlololo[i][idxval][j];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
         {
            cout << jacvalhihihi[0][i][j]
                 << "  " << jacvallohihi[0][i][j] << endl;
            cout << jacvalhilohi[0][i][j]
                 << "  " << jacvallolohi[0][i][j] << endl;
            cout << jacvalhihilo[0][i][j]
                 << "  " << jacvallohilo[0][i][j] << endl;
            cout << jacvalhilolo[0][i][j]
                 << "  " << jacvallololo[0][i][j] << endl;
         }
      }
   }
}

void cmplx8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbrehihihi, double **mbrelohihi,
   double **mbrehilohi, double **mbrelolohi,
   double **mbrehihilo, double **mbrelohilo,
   double **mbrehilolo, double **mbrelololo,
   double **mbimhihihi, double **mbimlohihi,
   double **mbimhilohi, double **mbimlolohi,
   double **mbimhihilo, double **mbimlohilo,
   double **mbimhilolo, double **mbimlololo, double damper,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo, int vrblvl )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehihihi[i][j] = outputrehihihi[i][dim][j];
         funvalrelohihi[i][j] = outputrelohihi[i][dim][j];
         funvalrehilohi[i][j] = outputrehilohi[i][dim][j];
         funvalrelolohi[i][j] = outputrelolohi[i][dim][j];
         funvalrehihilo[i][j] = outputrehihilo[i][dim][j];
         funvalrelohilo[i][j] = outputrelohilo[i][dim][j];
         funvalrehilolo[i][j] = outputrehilolo[i][dim][j];
         funvalrelololo[i][j] = outputrelololo[i][dim][j];
         funvalimhihihi[i][j] = outputimhihihi[i][dim][j];
         funvalimlohihi[i][j] = outputimlohihi[i][dim][j];
         funvalimhilohi[i][j] = outputimhilohi[i][dim][j];
         funvalimlolohi[i][j] = outputimlolohi[i][dim][j];
         funvalimhihilo[i][j] = outputimhihilo[i][dim][j];
         funvalimlohilo[i][j] = outputimlohilo[i][dim][j];
         funvalimhilolo[i][j] = outputimhilolo[i][dim][j];
         funvalimlololo[i][j] = outputimlololo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalrehihihi[i][0] << "  " << funvalrelohihi[i][0] << endl
              << "  "
              << funvalrehilohi[i][0] << "  " << funvalrelolohi[i][0] << endl
              << "  "
              << funvalrehihilo[i][0] << "  " << funvalrelohilo[i][0] << endl
              << "  "
              << funvalrehilolo[i][0] << "  " << funvalrelololo[i][0] << endl
              << "  "
              << funvalimhihihi[i][0] << "  " << funvalimlohihi[i][0] << endl
              << "  "
              << funvalimhilohi[i][0] << "  " << funvalimlolohi[i][0] << endl
              << "  "
              << funvalimhihilo[i][0] << "  " << funvalimlohilo[i][0] << endl
              << "  "
              << funvalimhilolo[i][0] << "  " << funvalimlololo[i][0] << endl;
   }
   /*
    * Linearize the function values in the rhs and swap sign,
    * but keep in mind that the right hand side is mbre + I*mbim,
    * as the monomial system is x^U = rhs, or x^U - rhs = 0,
    * so from the funval we substract rhs and then flip sign
    * in the computation of the rhs of the linear system.
    */
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         // rhsre[i][j] = -(funvalre[j][i] - mbre[j][i]);
         odf_sub(funvalrehihihi[j][i],funvalrelohihi[j][i],
                 funvalrehilohi[j][i],funvalrelolohi[j][i],
                 funvalrehihilo[j][i],funvalrelohilo[j][i],
                 funvalrehilolo[j][i],funvalrelololo[j][i],
                 mbrehihihi[j][i],mbrelohihi[j][i],
                 mbrehilohi[j][i],mbrelolohi[j][i],
                 mbrehihilo[j][i],mbrelohilo[j][i],
                 mbrehilolo[j][i],mbrelololo[j][i],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhsrehihihi[i][j] = -acchihihi;
         rhsrelohihi[i][j] = -acclohihi;
         rhsrehilohi[i][j] = -acchilohi;
         rhsrelolohi[i][j] = -acclolohi;
         rhsrehihilo[i][j] = -acchihilo;
         rhsrelohilo[i][j] = -acclohilo;
         rhsrehilolo[i][j] = -acchilolo;
         rhsrelololo[i][j] = -acclololo;
         // rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
         odf_sub(funvalimhihihi[j][i],funvalimlohihi[j][i],
                 funvalimhilohi[j][i],funvalimlolohi[j][i],
                 funvalimhihilo[j][i],funvalimlohilo[j][i],
                 funvalimhilolo[j][i],funvalimlololo[j][i],
                 mbimhihihi[j][i],mbimlohihi[j][i],
                 mbimhilohi[j][i],mbimlolohi[j][i],
                 mbimhihilo[j][i],mbimlohilo[j][i],
                 mbimhilolo[j][i],mbimlololo[j][i],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhsimhihihi[i][j] = -acchihihi;
         rhsimlohihi[i][j] = -acclohihi;
         rhsimhilohi[i][j] = -acchilohi;
         rhsimlolohi[i][j] = -acclolohi;
         rhsimhihilo[i][j] = -acchihilo;
         rhsimlohilo[i][j] = -acclohilo;
         rhsimhilolo[i][j] = -acchilolo;
         rhsimlololo[i][j] = -acclololo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsrehihihi[i][j] << "  " << rhsrelohihi[i][j] << endl
                 << "  "
                 << rhsrehilohi[i][j] << "  " << rhsrelolohi[i][j] << endl
                 << "  "
                 << rhsrehihilo[i][j] << "  " << rhsrelohilo[i][j] << endl
                 << "  "
                 << rhsrehilolo[i][j] << "  " << rhsrelololo[i][j] << endl
                 << "  "
                 << rhsimhihihi[i][j] << "  " << rhsimlohihi[i][j] << endl
                 << "  "
                 << rhsimhilohi[i][j] << "  " << rhsimlolohi[i][j] << endl
                 << "  "
                 << rhsimhihilo[i][j] << "  " << rhsimlohilo[i][j] << endl
                 << "  "
                 << rhsimhilolo[i][j] << "  " << rhsimlololo[i][j] << endl;
      }
   }
   // output[i][idx[i][k]] for k from 0 to nvr[i] has the coefficients
   // of the series of the derivative of the i-th monomial with respect
   // to the variable idx[i][k].
   for(int i=0; i<dim; i++)          // the i-th monomial
   {
      int *indexes = idx[i];         // indices of the variables
      for(int k=0; k<nvr[i]; k++)    // derivative w.r.t. idx[i][k]
      {                              // has j-th coefficient
         int idxval = indexes[k];

         for(int j=0; j<degp1; j++) 
         {
            jacvalrehihihi[j][i][idxval] = outputrehihihi[i][idxval][j];
            jacvalrelohihi[j][i][idxval] = outputrelohihi[i][idxval][j];
            jacvalrehilohi[j][i][idxval] = outputrehilohi[i][idxval][j];
            jacvalrelolohi[j][i][idxval] = outputrelolohi[i][idxval][j];
            jacvalrehihilo[j][i][idxval] = outputrehihilo[i][idxval][j];
            jacvalrelohilo[j][i][idxval] = outputrelohilo[i][idxval][j];
            jacvalrehilolo[j][i][idxval] = outputrehilolo[i][idxval][j];
            jacvalrelololo[j][i][idxval] = outputrelololo[i][idxval][j];
            jacvalimhihihi[j][i][idxval] = outputimhihihi[i][idxval][j];
            jacvalimlohihi[j][i][idxval] = outputimlohihi[i][idxval][j];
            jacvalimhilohi[j][i][idxval] = outputimhilohi[i][idxval][j];
            jacvalimlolohi[j][i][idxval] = outputimlolohi[i][idxval][j];
            jacvalimhihilo[j][i][idxval] = outputimhihilo[i][idxval][j];
            jacvalimlohilo[j][i][idxval] = outputimlohilo[i][idxval][j];
            jacvalimhilolo[j][i][idxval] = outputimhilolo[i][idxval][j];
            jacvalimlololo[j][i][idxval] = outputimlololo[i][idxval][j];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalrehihihi[0][i][j] << "  "
                 << jacvalrelohihi[0][i][j] << endl << "  "
                 << jacvalrehilohi[0][i][j] << "  "
                 << jacvalrelolohi[0][i][j] << endl << "  "
                 << jacvalrehihilo[0][i][j] << "  "
                 << jacvalrelohilo[0][i][j] << endl << "  "
                 << jacvalrehilolo[0][i][j] << "  "
                 << jacvalrelololo[0][i][j] << endl << "  "
                 << jacvalimhihihi[0][i][j] << "  "
                 << jacvalimlohihi[0][i][j] << endl << "  "
                 << jacvalimhilohi[0][i][j] << "  "
                 << jacvalimlolohi[0][i][j] << endl << "  "
                 << jacvalimhihilo[0][i][j] << "  "
                 << jacvalimlohilo[0][i][j] << endl << "  "
                 << jacvalimhilolo[0][i][j] << "  "
                 << jacvalimlololo[0][i][j] << endl;
      }
   }
}

void dbl8_define_rhs
 ( int dim, int degp1,
   double **mbhihihi, double **mblohihi, double **mbhilohi, double **mblolohi,
   double **mbhihilo, double **mblohilo, double **mbhilolo, double **mblololo,
   double **funvalhihihi, double **funvallohihi,
   double **funvalhilohi, double **funvallolohi,
   double **funvalhihilo, double **funvallohilo,
   double **funvalhilolo, double **funvallololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo, int vrblvl )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalhihihi[i][0] << "  " << funvallohihi[i][0] << endl
              << "  "
              << funvalhilohi[i][0] << "  " << funvallolohi[i][0] << endl
              << "  "
              << funvalhihilo[i][0] << "  " << funvallohilo[i][0] << endl
              << "  "
              << funvalhilolo[i][0] << "  " << funvallololo[i][0] << endl;
   }
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++) // rhs[i][j] = -(funval[j][i] - mb[j][i]);
      {
         odf_sub(funvalhihihi[j][i],funvallohihi[j][i],
                 funvalhilohi[j][i],funvallolohi[j][i],
                 funvalhihilo[j][i],funvallohilo[j][i],
                 funvalhilolo[j][i],funvallololo[j][i],
                 mbhihihi[j][i],mblohihi[j][i],mbhilohi[j][i],mblolohi[j][i],
                 mbhihilo[j][i],mblohilo[j][i],mbhilolo[j][i],mblololo[j][i],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhshihihi[i][j] = -acchihihi;
         rhslohihi[i][j] = -acclohihi;
         rhshilohi[i][j] = -acchilohi;
         rhslolohi[i][j] = -acclolohi;
         rhshihilo[i][j] = -acchihilo;
         rhslohilo[i][j] = -acclohilo;
         rhshilolo[i][j] = -acchilolo;
         rhslololo[i][j] = -acclololo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhshihihi[i][j] << "  " << rhslohihi[i][j] << endl
                 << "  "
                 << rhshilohi[i][j] << "  " << rhslolohi[i][j] << endl
                 << "  "
                 << rhshihilo[i][j] << "  " << rhslohilo[i][j] << endl
                 << "  "
                 << rhshilolo[i][j] << "  " << rhslololo[i][j] << endl;
      }
   }
}

void cmplx8_define_rhs
 ( int dim, int degp1,
   double **mbrehihihi, double **mbrelohihi,
   double **mbrehilohi, double **mbrelolohi,
   double **mbrehihilo, double **mbrelohilo,
   double **mbrehilolo, double **mbrelololo,
   double **mbimhihihi, double **mbimlohihi,
   double **mbimhilohi, double **mbimlolohi,
   double **mbimhihilo, double **mbimlohilo,
   double **mbimhilolo, double **mbimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double **rhsrehihihi, double **rhsrelohihi,
   double **rhsrehilohi, double **rhsrelolohi,
   double **rhsrehihilo, double **rhsrelohilo,
   double **rhsrehilolo, double **rhsrelololo,
   double **rhsimhihihi, double **rhsimlohihi,
   double **rhsimhilohi, double **rhsimlolohi,
   double **rhsimhihilo, double **rhsimlohilo,
   double **rhsimhilolo, double **rhsimlololo, int vrblvl )
{
   double acchihihi,acclohihi,acchilohi,acclolohi;
   double acchihilo,acclohilo,acchilolo,acclololo;

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " 
              << funvalrehihihi[i][0] << "  " << funvalrelohihi[i][0] << endl
              << "  " 
              << funvalrehilohi[i][0] << "  " << funvalrelolohi[i][0] << endl
              << "  " 
              << funvalrehihilo[i][0] << "  " << funvalrelohilo[i][0] << endl
              << "  " 
              << funvalrehilolo[i][0] << "  " << funvalrelololo[i][0] << endl
              << "  " 
              << funvalimhihihi[i][0] << "  " << funvalimlohihi[i][0] << endl
              << "  " 
              << funvalimhilohi[i][0] << "  " << funvalimlolohi[i][0] << endl
              << "  " 
              << funvalimhihilo[i][0] << "  " << funvalimlohilo[i][0] << endl
              << "  " 
              << funvalimhilolo[i][0] << "  " << funvalimlololo[i][0] << endl;
   }
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         // rhsre[i][j] = -(funvalre[j][i] - mbre[j][i]);
         odf_sub(funvalrehihihi[j][i],funvalrelohihi[j][i],
                 funvalrehilohi[j][i],funvalrelolohi[j][i],
                 funvalrehihilo[j][i],funvalrelohilo[j][i],
                 funvalrehilolo[j][i],funvalrelololo[j][i],
                 mbrehihihi[j][i],mbrelohihi[j][i],
                 mbrehilohi[j][i],mbrelolohi[j][i],
                 mbrehihilo[j][i],mbrelohilo[j][i],
                 mbrehilolo[j][i],mbrelololo[j][i],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhsrehihihi[i][j] = -acchihihi;
         rhsrelohihi[i][j] = -acclohihi;
         rhsrehilohi[i][j] = -acchilohi;
         rhsrelolohi[i][j] = -acclolohi;
         rhsrehihilo[i][j] = -acchihilo;
         rhsrelohilo[i][j] = -acclohilo;
         rhsrehilolo[i][j] = -acchilolo;
         rhsrelololo[i][j] = -acclololo;
         // rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
         odf_sub(funvalimhihihi[j][i],funvalimlohihi[j][i],
                 funvalimhilohi[j][i],funvalimlolohi[j][i],
                 funvalimhihilo[j][i],funvalimlohilo[j][i],
                 funvalimhilolo[j][i],funvalimlololo[j][i],
                 mbimhihihi[j][i],mbimlohihi[j][i],
                 mbimhilohi[j][i],mbimlolohi[j][i],
                 mbimhihilo[j][i],mbimlohilo[j][i],
                 mbimhilolo[j][i],mbimlololo[j][i],
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhsimhihihi[i][j] = -acchihihi;
         rhsimlohihi[i][j] = -acclohihi;
         rhsimhilolo[i][j] = -acchilolo;
         rhsimlololo[i][j] = -acclololo;
         rhsimhihihi[i][j] = -acchihihi;
         rhsimlohihi[i][j] = -acclohihi;
         rhsimhilolo[i][j] = -acchilolo;
         rhsimlololo[i][j] = -acclololo;
      }
 
   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsrehihihi[i][j] << "  " << rhsrelohihi[i][j] << endl
                 << "  "
                 << rhsrehilohi[i][j] << "  " << rhsrelolohi[i][j] << endl
                 << "  "
                 << rhsrehihilo[i][j] << "  " << rhsrelohilo[i][j] << endl
                 << "  "
                 << rhsrehilolo[i][j] << "  " << rhsrelololo[i][j] << endl
                 << "  "
                 << rhsimhihihi[i][j] << "  " << rhsimlohihi[i][j] << endl
                 << "  "
                 << rhsimhilohi[i][j] << "  " << rhsimlolohi[i][j] << endl
                 << "  "
                 << rhsimhihilo[i][j] << "  " << rhsimlohilo[i][j] << endl
                 << "  "
                 << rhsimhilolo[i][j] << "  " << rhsimlololo[i][j] << endl;
      }
   }
}

void dbl8_map_evaldiff_output
 ( int dim, int deg,
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
   double ***jacvalhilolo, double ***jacvallololo, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         funvalhihihi[i][j] = outputhihihi[i][dim][j];
         funvallohihi[i][j] = outputlohihi[i][dim][j];
         funvalhilohi[i][j] = outputhilohi[i][dim][j];
         funvallolohi[i][j] = outputlolohi[i][dim][j];
         funvalhihilo[i][j] = outputhihilo[i][dim][j];
         funvallohilo[i][j] = outputlohilo[i][dim][j];
         funvalhilolo[i][j] = outputhilolo[i][dim][j];
         funvallololo[i][j] = outputlololo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the evaluated series :" << endl;

      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalhihihi[i][0] << "  "
                            << funvallohihi[i][0] << endl
                            << funvalhilohi[i][0] << "  "
                            << funvallolohi[i][0] << endl
                            << funvalhihilo[i][0] << "  "
                            << funvallohilo[i][0] << endl
                            << funvalhilolo[i][0] << "  "
                            << funvallololo[i][0] << endl;
   }
   // output[i][j][k] is the k-th coefficient in the series
   // the derivative of the i-th polynomial with respect to variable j.

   for(int i=0; i<dim; i++)          // the i-th polynomial
   {
      for(int j=0; j<dim; j++)       // derivative w.r.t. j-th variable
      {
         for(int k=0; k<=deg; k++) 
         {
            jacvalhihihi[k][i][j] = outputhihihi[i][j][k];
            jacvallohihi[k][i][j] = outputlohihi[i][j][k];
            jacvalhilohi[k][i][j] = outputhilohi[i][j][k];
            jacvallolohi[k][i][j] = outputlolohi[i][j][k];
            jacvalhihilo[k][i][j] = outputhihilo[i][j][k];
            jacvallohilo[k][i][j] = outputlohilo[i][j][k];
            jacvalhilolo[k][i][j] = outputhilolo[i][j][k];
            jacvallololo[k][i][j] = outputlololo[i][j][k];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalhihihi[0][i][j] << "  "
                 << jacvallohihi[0][i][j] << endl
                 << jacvalhilohi[0][i][j] << "  "
                 << jacvallolohi[0][i][j] << endl
                 << jacvalhihilo[0][i][j] << "  "
                 << jacvallohilo[0][i][j] << endl
                 << jacvalhilolo[0][i][j] << "  "
                 << jacvallololo[0][i][j] << endl;
      }
   }
}

void cmplx8_map_evaldiff_output
 ( int dim, int deg,
   double ***outputrehihihi, double ***outputrelohihi,
   double ***outputrehilohi, double ***outputrelolohi,
   double ***outputrehihilo, double ***outputrelohilo,
   double ***outputrehilolo, double ***outputrelololo,
   double ***outputimhihihi, double ***outputimlohihi,
   double ***outputimhilohi, double ***outputimlolohi,
   double ***outputimhihilo, double ***outputimlohilo,
   double ***outputimhilolo, double ***outputimlololo,
   double **funvalrehihihi, double **funvalrelohihi,
   double **funvalrehilohi, double **funvalrelolohi,
   double **funvalrehihilo, double **funvalrelohilo,
   double **funvalrehilolo, double **funvalrelololo,
   double **funvalimhihihi, double **funvalimlohihi,
   double **funvalimhilohi, double **funvalimlolohi,
   double **funvalimhihilo, double **funvalimlohilo,
   double **funvalimhilolo, double **funvalimlololo,
   double ***jacvalrehihihi, double ***jacvalrelohihi,
   double ***jacvalrehilohi, double ***jacvalrelolohi,
   double ***jacvalrehihilo, double ***jacvalrelohilo,
   double ***jacvalrehilolo, double ***jacvalrelololo,
   double ***jacvalimhihihi, double ***jacvalimlohihi,
   double ***jacvalimhilohi, double ***jacvalimlolohi,
   double ***jacvalimhihilo, double ***jacvalimlohilo,
   double ***jacvalimhilolo, double ***jacvalimlololo, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         funvalrehihihi[i][j] = outputrehihihi[i][dim][j];
         funvalrelohihi[i][j] = outputrelohihi[i][dim][j];
         funvalrehilohi[i][j] = outputrehilohi[i][dim][j];
         funvalrelolohi[i][j] = outputrelolohi[i][dim][j];
         funvalrehihilo[i][j] = outputrehihilo[i][dim][j];
         funvalrelohilo[i][j] = outputrelohilo[i][dim][j];
         funvalrehilolo[i][j] = outputrehilolo[i][dim][j];
         funvalrelololo[i][j] = outputrelololo[i][dim][j];
         funvalimhihihi[i][j] = outputimhihihi[i][dim][j];
         funvalimlohihi[i][j] = outputimlohihi[i][dim][j];
         funvalimhilohi[i][j] = outputimhilohi[i][dim][j];
         funvalimlolohi[i][j] = outputimlolohi[i][dim][j];
         funvalimhihilo[i][j] = outputimhihilo[i][dim][j];
         funvalimlohilo[i][j] = outputimlohilo[i][dim][j];
         funvalimhilolo[i][j] = outputimhilolo[i][dim][j];
         funvalimlololo[i][j] = outputimlololo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the evaluated series :" << endl;

      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalrehihihi[i][0] << "  "
                            << funvalrelohihi[i][0] << endl
                            << funvalrehilohi[i][0] << "  "
                            << funvalrelolohi[i][0] << endl
                            << funvalrehihilo[i][0] << "  "
                            << funvalrelohilo[i][0] << endl
                            << funvalrehilolo[i][0] << "  "
                            << funvalrelololo[i][0] << endl
                            << funvalimhihihi[i][0] << "  "
                            << funvalimlohihi[i][0] << endl
                            << funvalimhilohi[i][0] << "  "
                            << funvalimlolohi[i][0] << endl
                            << funvalimhihilo[i][0] << "  "
                            << funvalimlohilo[i][0] << endl
                            << funvalimhilolo[i][0] << "  "
                            << funvalimlololo[i][0] << endl;
   }
   // output[i][j][k] is the k-th coefficient in the series
   // the derivative of the i-th polynomial with respect to variable j.

   for(int i=0; i<dim; i++)          // the i-th polynomial
   {
      for(int j=0; j<dim; j++)       // derivative w.r.t. j-th variable
      {
         for(int k=0; k<=deg; k++) 
         {
            jacvalrehihihi[k][i][j] = outputrehihihi[i][j][k];
            jacvalrelohihi[k][i][j] = outputrelohihi[i][j][k];
            jacvalrehilohi[k][i][j] = outputrehilohi[i][j][k];
            jacvalrelolohi[k][i][j] = outputrelolohi[i][j][k];
            jacvalrehihilo[k][i][j] = outputrehihilo[i][j][k];
            jacvalrelohilo[k][i][j] = outputrelohilo[i][j][k];
            jacvalrehilolo[k][i][j] = outputrehilolo[i][j][k];
            jacvalrelololo[k][i][j] = outputrelololo[i][j][k];
            jacvalimhihihi[k][i][j] = outputimhihihi[i][j][k];
            jacvalimlohihi[k][i][j] = outputimlohihi[i][j][k];
            jacvalimhilohi[k][i][j] = outputimhilohi[i][j][k];
            jacvalimlolohi[k][i][j] = outputimlolohi[i][j][k];
            jacvalimhihilo[k][i][j] = outputimhihilo[i][j][k];
            jacvalimlohilo[k][i][j] = outputimlohilo[i][j][k];
            jacvalimhilolo[k][i][j] = outputimhilolo[i][j][k];
            jacvalimlololo[k][i][j] = outputimlololo[i][j][k];
         }
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacvalrehihihi[0][i][j] << "  "
                 << jacvalrelohihi[0][i][j] << endl
                 << jacvalrehilohi[0][i][j] << "  "
                 << jacvalrelolohi[0][i][j] << endl
                 << jacvalrehihilo[0][i][j] << "  "
                 << jacvalrelohilo[0][i][j] << endl
                 << jacvalrehilolo[0][i][j] << "  "
                 << jacvalrelololo[0][i][j] << endl
                 << jacvalimhihihi[0][i][j] << "  "
                 << jacvalimlohihi[0][i][j] << endl
                 << jacvalimhilohi[0][i][j] << "  "
                 << jacvalimlolohi[0][i][j] << endl
                 << jacvalimhihilo[0][i][j] << "  "
                 << jacvalimlohilo[0][i][j] << endl
                 << jacvalimhilolo[0][i][j] << "  "
                 << jacvalimlololo[0][i][j] << endl;
      }
   }
}
