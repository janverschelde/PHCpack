// The file dbl4_systems_host.cpp defines the functions with prototypes in
// the file dbl4_systems_host.h.

#include <iostream>
#include "quad_double_functions.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"

using namespace std;

void CPU_dbl4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo, int vrblvl )
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
                  CPU_dbl4_product
                     (deg,inputhihi[idxvar],inputlohi[idxvar],
                          inputhilo[idxvar],inputlolo[idxvar],
                      cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
                      acchihi,   acclohi,   acchilo,   acclolo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffhihi[i][L] = acchihi[L];
                     cfflohi[i][L] = acclohi[L];
                     cffhilo[i][L] = acchilo[L];
                     cfflolo[i][L] = acclolo[L];
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
            cout << cffhihi[i][j] << "  " << cfflohi[i][j] << endl;
            cout << cffhilo[i][j] << "  " << cfflolo[i][j] << endl;
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
            cout << inputhihi[i][j] << "  " << inputlohi[i][j] << endl;
            cout << inputhilo[i][j] << "  " << inputlolo[i][j] << endl;
         }
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl4_evaldiff
         (dim,nvr[i],deg,idx[i],cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
           inputhihi,    inputlohi,    inputhilo,    inputlolo,
          outputhihi[i],outputlohi[i],outputhilo[i],outputlolo[i]);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
         {
            cout << outputhihi[i][dim][j] << "  "
                 << outputlohi[i][dim][j] << endl;
            cout << outputhilo[i][dim][j] << "  "
                 << outputlolo[i][dim][j] << endl;
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
               double tmphihi,tmplohi,tmphilo,tmplolo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // output[i][idxvar][k] = factor*output[i][idxvar][k];
                  qdf_mul_qd_d
                     (outputhihi[i][idxvar][k],outputlohi[i][idxvar][k],
                      outputhilo[i][idxvar][k],outputlolo[i][idxvar][k],
                      factor,&tmphihi,&tmplohi,&tmphilo,&tmplolo);
                  outputhihi[i][idxvar][k] = tmphihi;
                  outputlohi[i][idxvar][k] = tmplohi;
                  outputhilo[i][idxvar][k] = tmphilo;
                  outputlolo[i][idxvar][k] = tmplolo;
               }
            }
         }
      }
   }
}

void CPU_cmplx4_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi, double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi, double *accimhilo, double *accimlolo,
   double **inputrehihi, double **inputrelohi,
   double **inputrehilo, double **inputrelolo,
   double **inputimhihi, double **inputimlohi, 
   double **inputimhilo, double **inputimlolo, 
   double ***outputrehihi, double ***outputrelohi, 
   double ***outputrehilo, double ***outputrelolo, 
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo, int vrblvl )
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
                  CPU_cmplx4_product
                     (deg,inputrehihi[idxvar],inputrelohi[idxvar],
                          inputrehilo[idxvar],inputrelolo[idxvar],
                          inputimhihi[idxvar],inputimlohi[idxvar],
                          inputimhilo[idxvar],inputimlolo[idxvar],
                      cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
                      cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
                      accrehihi,accrelohi,accrehilo,accrelolo,
                      accimhihi,accimlohi,accimhilo,accimlolo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffrehihi[i][L] = accrehihi[L];
                     cffrelohi[i][L] = accrelohi[L];
                     cffrehilo[i][L] = accrehilo[L];
                     cffrelolo[i][L] = accrelolo[L];
                     cffimhihi[i][L] = accimhihi[L];
                     cffimlohi[i][L] = accimlohi[L];
                     cffimhilo[i][L] = accimhilo[L];
                     cffimlolo[i][L] = accimlolo[L];
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
            cout << cffrehihi[i][j] << "  " << cffrelohi[i][j] << endl << "  "
                 << cffrehilo[i][j] << "  " << cffrelolo[i][j] << endl << "  "
                 << cffimhihi[i][j] << "  " << cffimlohi[i][j] << endl << "  "
                 << cffimhilo[i][j] << "  " << cffimlolo[i][j] << endl;
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
            cout << inputrehihi[i][j] << "  " << inputrelohi[i][j] << endl
                 << "  "
                 << inputrehilo[i][j] << "  " << inputrelolo[i][j] << endl
                 << "  "
                 << inputimhihi[i][j] << "  " << inputimlohi[i][j] << endl
                 << "  "
                 << inputimhilo[i][j] << "  " << inputimlolo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_cmplx4_evaldiff
         (dim,nvr[i],deg,idx[i],
          cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
          cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
          inputrehihi,inputrelohi,inputrehilo,inputrelolo,
          inputimhihi,inputimlohi,inputimhilo,inputimlolo,
          outputrehihi[i],outputrelohi[i],outputrehilo[i],outputrelolo[i],
          outputimhihi[i],outputimlohi[i],outputimhilo[i],outputimlolo[i]);

   if(vrblvl > 0)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehihi[i][dim][j] << "  "
                 << outputrelohi[i][dim][j] << endl << "  "
                 << outputrehilo[i][dim][j] << "  "
                 << outputrelolo[i][dim][j] << endl << "  "
                 << outputimhihi[i][dim][j] << "  "
                 << outputimlohi[i][dim][j] << endl << "  "
                 << outputimhilo[i][dim][j] << "  "
                 << outputimlolo[i][dim][j] << endl;
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
               double acchihi,acclohi,acchilo,acclolo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // outputre[i][idxvar][k] = factor*outputre[i][idxvar][k];
                  qdf_mul_qd_d
                     (outputrehihi[i][idxvar][k],outputrelohi[i][idxvar][k],
                      outputrehilo[i][idxvar][k],outputrelolo[i][idxvar][k],
                      factor,&acchihi,&acclohi,&acchilo,&acclolo);
                  outputrehihi[i][idxvar][k] = acchihi;
                  outputrelohi[i][idxvar][k] = acclohi;
                  outputrehilo[i][idxvar][k] = acchilo;
                  outputrelolo[i][idxvar][k] = acclolo;
                  // outputim[i][idxvar][k] = factor*outputim[i][idxvar][k];
                  qdf_mul_qd_d
                     (outputimhihi[i][idxvar][k],outputimlohi[i][idxvar][k],
                      outputimhilo[i][idxvar][k],outputimlolo[i][idxvar][k],
                      factor,&acchihi,&acclohi,&acchilo,&acclolo);
                  outputimhihi[i][idxvar][k] = acchihi;
                  outputimlohi[i][idxvar][k] = acclohi;
                  outputimhilo[i][idxvar][k] = acchilo;
                  outputimlolo[i][idxvar][k] = acclolo;
               }
            }
         }
      }
   }
   if(vrblvl > 0)
   {
      cout << "after multiplication with the factors ..." << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehihi[i][dim][j] << "  "
                 << outputrelohi[i][dim][j] << endl << "  "
                 << outputrehilo[i][dim][j] << "  "
                 << outputrelolo[i][dim][j] << endl << "  "
                 << outputimhihi[i][dim][j] << "  "
                 << outputimlohi[i][dim][j] << endl << "  "
                 << outputimhilo[i][dim][j] << "  "
                 << outputimlolo[i][dim][j] << endl;
      }
   }
}

void dbl4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi, 
   double **funvalhilo, double **funvallolo, 
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhihi[i][j] = outputhihi[i][dim][j];
         funvallohi[i][j] = outputlohi[i][dim][j];
         funvalhilo[i][j] = outputhilo[i][dim][j];
         funvallolo[i][j] = outputlolo[i][dim][j];
      }

   if(vrblvl > 0)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << funvalhihi[i][0] << "  "
                            << funvallohi[i][0] << endl;
         cout << "     " << funvalhilo[i][0] << "  "
                         << funvallolo[i][0] << endl;
      }
   }
   // Linearize the function values in the rhs and swap sign,
   // but keep in mind that the right hand side is 1 - t,
   // so we subtract 1 and add t to the rhs.
   for(int j=0; j<dim; j++) 
   {                                 // rhs[0][j] = -(funval[j][0] - 1.0);
      rhshihi[0][j] = -funvalhihi[j][0];
      rhslohi[0][j] = -funvallohi[j][0];
      rhshilo[0][j] = -funvalhilo[j][0];
      rhslolo[0][j] = -funvallolo[j][0];
      qdf_inc_d(&rhshihi[0][j],&rhslohi[0][j],
                &rhshilo[0][j],&rhslolo[0][j],1.0);
   }
   if(degp1 > 1)
   {
      for(int j=0; j<dim; j++)       // rhs[1][j] = -(funval[j][1] + 1.0);
      {
         rhshihi[1][j] = -funvalhihi[j][1];
         rhslohi[1][j] = -funvallohi[j][1];
         rhshilo[1][j] = -funvalhilo[j][1];
         rhslolo[1][j] = -funvallolo[j][1];
         qdf_dec(&rhshihi[1][j],&rhslohi[1][j],&rhshilo[1][j],&rhslolo[1][j],
                 1.0,0.0,0.0,0.0);
      }
      for(int i=2; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshihi[i][j] = -funvalhihi[j][i];
            rhslohi[i][j] = -funvallohi[j][i];
            rhshilo[i][j] = -funvalhilo[j][i];
            rhslolo[i][j] = -funvallolo[j][i];
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
            cout << rhshihi[i][j] << "  " << rhslohi[i][j] << endl;
            cout << rhshilo[i][j] << "  " << rhslolo[i][j] << endl;
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
            jacvalhihi[j][i][idxval] = outputhihi[i][idxval][j];
            jacvallohi[j][i][idxval] = outputlohi[i][idxval][j];
            jacvalhilo[j][i][idxval] = outputhilo[i][idxval][j];
            jacvallolo[j][i][idxval] = outputlolo[i][idxval][j];
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
            cout << jacvalhihi[0][i][j]
                 << "  " << jacvallohi[0][i][j] << endl;
            cout << jacvalhilo[0][i][j]
                 << "  " << jacvallolo[0][i][j] << endl;
         }
      }
   }
}

void cmplx4_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double ***outputrehihi, double ***outputrelohi,
   double ***outputrehilo, double ***outputrelolo,
   double ***outputimhihi, double ***outputimlohi,
   double ***outputimhilo, double ***outputimlolo,
   double **funvalrehihi, double **funvalrelohi,
   double **funvalrehilo, double **funvalrelolo,
   double **funvalimhihi, double **funvalimlohi,
   double **funvalimhilo, double **funvalimlolo,
   double **rhsrehihi, double **rhsrelohi,
   double **rhsrehilo, double **rhsrelolo,
   double **rhsimhihi, double **rhsimlohi,
   double **rhsimhilo, double **rhsimlolo,
   double ***jacvalrehihi, double ***jacvalrelohi,
   double ***jacvalrehilo, double ***jacvalrelolo,
   double ***jacvalimhihi, double ***jacvalimlohi,
   double ***jacvalimhilo, double ***jacvalimlolo, int vrblvl )
{
   double acchihi,acclohi,acchilo,acclolo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehihi[i][j] = outputrehihi[i][dim][j];
         funvalrelohi[i][j] = outputrelohi[i][dim][j];
         funvalrehilo[i][j] = outputrehilo[i][dim][j];
         funvalrelolo[i][j] = outputrelolo[i][dim][j];
         funvalimhihi[i][j] = outputimhihi[i][dim][j];
         funvalimlohi[i][j] = outputimlohi[i][dim][j];
         funvalimhilo[i][j] = outputimhilo[i][dim][j];
         funvalimlolo[i][j] = outputimlolo[i][dim][j];
      }

   if(vrblvl > 0)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalrehihi[i][0] << "  " << funvalrelohi[i][0] << endl
              << "  "
              << funvalrehilo[i][0] << "  " << funvalrelolo[i][0] << endl
              << "  "
              << funvalimhihi[i][0] << "  " << funvalimlohi[i][0] << endl
              << "  "
              << funvalimhilo[i][0] << "  " << funvalimlolo[i][0] << endl;
   }
   // Linearize the function values in the rhs and swap sign,
   // but keep in mind that the right hand side is 1 - t,
   // so we subtract 1 and add t to the rhs.
   for(int j=0; j<dim; j++)
   {
      // rhsre[0][j] = -(funvalre[j][0] - 1.0);
      qdf_sub(funvalrehihi[j][0],funvalrelohi[j][0],
              funvalrehilo[j][0],funvalrelolo[j][0],1.0,0.0,0.0,0.0,
              &acchihi,&acclohi,&acchilo,&acclolo);
      rhsrehihi[0][j] = -acchihi;
      rhsrelohi[0][j] = -acclohi;
      rhsrehilo[0][j] = -acchilo;
      rhsrelolo[0][j] = -acclolo;
      // rhsim[0][j] = -(funvalim[j][0] - 0.0);
      rhsimhihi[0][j] = -funvalimhihi[j][0];
      rhsimlohi[0][j] = -funvalimlohi[j][0];
      rhsimhilo[0][j] = -funvalimhilo[j][0];
      rhsimlolo[0][j] = -funvalimlolo[j][0];
   }
   if(degp1 > 1)
   {
      for(int j=0; j<dim; j++)
      {
         // rhsre[1][j] = -(funvalre[j][1] + 1.0);
         qdf_add(funvalrehihi[j][1],funvalrelohi[j][1],
                 funvalrehilo[j][1],funvalrelolo[j][1],1.0,0.0,0.0,0.0,
                 &acchihi,&acclohi,&acchilo,&acclolo);
         rhsrehihi[1][j] = -acchihi;
         rhsrelohi[1][j] = -acclohi;
         rhsrehilo[1][j] = -acchilo;
         rhsrelolo[1][j] = -acclolo;
         // rhsim[1][j] = -(funvalim[j][1] + 0.0);
         rhsimhihi[1][j] = -funvalimhihi[j][1];
         rhsimlohi[1][j] = -funvalimlohi[j][1];
         rhsimhilo[1][j] = -funvalimhilo[j][1];
         rhsimlolo[1][j] = -funvalimlolo[j][1];
      }
      for(int i=2; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsrehihi[i][j] = -funvalrehihi[j][i];
            rhsrelohi[i][j] = -funvalrelohi[j][i];
            rhsrehilo[i][j] = -funvalrehilo[j][i];
            rhsrelolo[i][j] = -funvalrelolo[j][i];
            rhsimhihi[i][j] = -funvalimhihi[j][i];
            rhsimlohi[i][j] = -funvalimlohi[j][i];
            rhsimhilo[i][j] = -funvalimhilo[j][i];
            rhsimlolo[i][j] = -funvalimlolo[j][i];
         }
   }
   if(vrblvl > 0)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsrehihi[i][j] << "  " << rhsrelohi[i][j] << endl << "  "
                 << rhsrehilo[i][j] << "  " << rhsrelolo[i][j] << endl << "  "
                 << rhsimhihi[i][j] << "  " << rhsimlohi[i][j] << endl << "  "
                 << rhsimhilo[i][j] << "  " << rhsimlolo[i][j] << endl;
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
            jacvalrehihi[j][i][idxval] = outputrehihi[i][idxval][j];
            jacvalrelohi[j][i][idxval] = outputrelohi[i][idxval][j];
            jacvalrehilo[j][i][idxval] = outputrehilo[i][idxval][j];
            jacvalrelolo[j][i][idxval] = outputrelolo[i][idxval][j];
            jacvalimhihi[j][i][idxval] = outputimhihi[i][idxval][j];
            jacvalimlohi[j][i][idxval] = outputimlohi[i][idxval][j];
            jacvalimhilo[j][i][idxval] = outputimhilo[i][idxval][j];
            jacvalimlolo[j][i][idxval] = outputimlolo[i][idxval][j];
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
            cout << jacvalrehihi[0][i][j] << "  "
                 << jacvalrelohi[0][i][j] << endl << "  "
                 << jacvalrehilo[0][i][j] << "  "
                 << jacvalrelolo[0][i][j] << endl << "  "
                 << jacvalimhihi[0][i][j] << "  "
                 << jacvalimlohi[0][i][j] << endl << "  "
                 << jacvalimhilo[0][i][j] << "  "
                 << jacvalimlolo[0][i][j] << endl;
      }
   }
}
