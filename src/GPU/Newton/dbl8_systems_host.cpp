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
