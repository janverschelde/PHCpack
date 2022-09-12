// The file dbl8_systems_host.cpp defines the functions with prototypes in
// the file dbl8_systems_host.h.

#include <iostream>
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
   if(vrblvl > 0)
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

   if(vrblvl > 0)
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
   if(vrblvl > 0)
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

   if(vrblvl > 0)
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
}

void dbl8_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
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

   if(vrblvl > 0)
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
   // Linearize the function values in the rhs and swap sign,
   // but keep in mind that the right hand side is 1 - t,
   // so we subtract 1 and add t to the rhs.
   for(int j=0; j<dim; j++) 
   {                                 // rhs[0][j] = -(funval[j][0] - 1.0);
      rhshihihi[0][j] = -funvalhihihi[j][0];
      rhslohihi[0][j] = -funvallohihi[j][0];
      rhshilohi[0][j] = -funvalhilohi[j][0];
      rhslolohi[0][j] = -funvallolohi[j][0];
      rhshihilo[0][j] = -funvalhihilo[j][0];
      rhslohilo[0][j] = -funvallohilo[j][0];
      rhshilolo[0][j] = -funvalhilolo[j][0];
      rhslololo[0][j] = -funvallololo[j][0];
      odf_inc_d(&rhshihihi[0][j],&rhslohihi[0][j],
                &rhshilohi[0][j],&rhslolohi[0][j],
                &rhshihilo[0][j],&rhslohilo[0][j],
                &rhshilolo[0][j],&rhslololo[0][j],1.0);
   }
   if(degp1 > 1)
   {
      for(int j=0; j<dim; j++)       // rhs[1][j] = -(funval[j][1] + 1.0);
      {
         rhshihihi[1][j] = -funvalhihihi[j][1];
         rhslohihi[1][j] = -funvallohihi[j][1];
         rhshilohi[1][j] = -funvalhilohi[j][1];
         rhslolohi[1][j] = -funvallolohi[j][1];
         rhshihilo[1][j] = -funvalhihilo[j][1];
         rhslohilo[1][j] = -funvallohilo[j][1];
         rhshilolo[1][j] = -funvalhilolo[j][1];
         rhslololo[1][j] = -funvallololo[j][1];
         odf_dec(&rhshihihi[1][j],&rhslohihi[1][j],
                 &rhshilohi[1][j],&rhslolohi[1][j],
                 &rhshihilo[1][j],&rhslohilo[1][j],
                 &rhshilolo[1][j],&rhslololo[1][j],
                 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
      }
      for(int i=2; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshihihi[i][j] = -funvalhihihi[j][i];
            rhslohihi[i][j] = -funvallohihi[j][i];
            rhshilohi[i][j] = -funvalhilohi[j][i];
            rhslolohi[i][j] = -funvallolohi[j][i];
            rhshihilo[i][j] = -funvalhihilo[j][i];
            rhslohilo[i][j] = -funvallohilo[j][i];
            rhshilolo[i][j] = -funvalhilolo[j][i];
            rhslololo[i][j] = -funvallololo[j][i];
         }
   }
   if(vrblvl > 0)
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
   if(vrblvl > 0)
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

   if(vrblvl > 0)
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
   // Linearize the function values in the rhs and swap sign,
   // but keep in mind that the right hand side is 1 - t,
   // so we subtract 1 and add t to the rhs.
   for(int j=0; j<dim; j++)
   {
      // rhsre[0][j] = -(funvalre[j][0] - 1.0);
      odf_sub(funvalrehihihi[j][0],funvalrelohihi[j][0],
              funvalrehilohi[j][0],funvalrelolohi[j][0],
              funvalrehihilo[j][0],funvalrelohilo[j][0],
              funvalrehilolo[j][0],funvalrelololo[j][0],
              1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
              &acchihihi,&acclohihi,&acchilohi,&acclolohi,
              &acchihilo,&acclohilo,&acchilolo,&acclololo);
      rhsrehihihi[0][j] = -acchihihi;
      rhsrelohihi[0][j] = -acclohihi;
      rhsrehilohi[0][j] = -acchilohi;
      rhsrelolohi[0][j] = -acclolohi;
      rhsrehihilo[0][j] = -acchihilo;
      rhsrelohilo[0][j] = -acclohilo;
      rhsrehilolo[0][j] = -acchilolo;
      rhsrelololo[0][j] = -acclololo;
      // rhsim[0][j] = -(funvalim[j][0] - 0.0);
      rhsimhihihi[0][j] = -funvalimhihihi[j][0];
      rhsimlohihi[0][j] = -funvalimlohihi[j][0];
      rhsimhilohi[0][j] = -funvalimhilohi[j][0];
      rhsimlolohi[0][j] = -funvalimlolohi[j][0];
      rhsimhihilo[0][j] = -funvalimhihilo[j][0];
      rhsimlohilo[0][j] = -funvalimlohilo[j][0];
      rhsimhilolo[0][j] = -funvalimhilolo[j][0];
      rhsimlololo[0][j] = -funvalimlololo[j][0];
   }
   if(degp1 > 1)
   {
      for(int j=0; j<dim; j++)
      {
         // rhsre[1][j] = -(funvalre[j][1] + 1.0);
         odf_add(funvalrehihihi[j][1],funvalrelohihi[j][1],
                 funvalrehilohi[j][1],funvalrelolohi[j][1],
                 funvalrehihilo[j][1],funvalrelohilo[j][1],
                 funvalrehilolo[j][1],funvalrelololo[j][1],
                 1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                 &acchihihi,&acclohihi,&acchilohi,&acclolohi,
                 &acchihilo,&acclohilo,&acchilolo,&acclololo);
         rhsrehihihi[1][j] = -acchihihi;
         rhsrelohihi[1][j] = -acclohihi;
         rhsrehilohi[1][j] = -acchilohi;
         rhsrelolohi[1][j] = -acclolohi;
         rhsrehihilo[1][j] = -acchihilo;
         rhsrelohilo[1][j] = -acclohilo;
         rhsrehilolo[1][j] = -acchilolo;
         rhsrelololo[1][j] = -acclololo;
         // rhsim[1][j] = -(funvalim[j][1] + 0.0);
         rhsimhihihi[1][j] = -funvalimhihihi[j][1];
         rhsimlohihi[1][j] = -funvalimlohihi[j][1];
         rhsimhilohi[1][j] = -funvalimhilohi[j][1];
         rhsimlolohi[1][j] = -funvalimlolohi[j][1];
         rhsimhihilo[1][j] = -funvalimhihilo[j][1];
         rhsimlohilo[1][j] = -funvalimlohilo[j][1];
         rhsimhilolo[1][j] = -funvalimhilolo[j][1];
         rhsimlololo[1][j] = -funvalimlololo[j][1];
      }
      for(int i=2; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsrehihihi[i][j] = -funvalrehihihi[j][i];
            rhsrelohihi[i][j] = -funvalrelohihi[j][i];
            rhsrehilohi[i][j] = -funvalrehilohi[j][i];
            rhsrelolohi[i][j] = -funvalrelolohi[j][i];
            rhsrehihilo[i][j] = -funvalrehihilo[j][i];
            rhsrelohilo[i][j] = -funvalrelohilo[j][i];
            rhsrehilolo[i][j] = -funvalrehilolo[j][i];
            rhsrelololo[i][j] = -funvalrelololo[j][i];
            rhsimhihihi[i][j] = -funvalimhihihi[j][i];
            rhsimlohihi[i][j] = -funvalimlohihi[j][i];
            rhsimhilohi[i][j] = -funvalimhilohi[j][i];
            rhsimlolohi[i][j] = -funvalimlolohi[j][i];
            rhsimhihilo[i][j] = -funvalimhihilo[j][i];
            rhsimlohilo[i][j] = -funvalimlohilo[j][i];
            rhsimhilolo[i][j] = -funvalimhilolo[j][i];
            rhsimlololo[i][j] = -funvalimlololo[j][i];
         }
   }
   if(vrblvl > 0)
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
   if(vrblvl > 0)
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
