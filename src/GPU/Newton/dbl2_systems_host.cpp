// The file dbl2_systems_host.cpp defines the functions with prototypes in
// the file dbl2_systems_host.h.

#include <iostream>
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"

using namespace std;

void CPU_dbl2_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhi, double **cfflo, double *acchi, double *acclo,
   double **inputhi, double **inputlo,
   double ***outputhi, double ***outputlo, int vrblvl )
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
                  CPU_dbl2_product
                     (deg,inputhi[idxvar],inputlo[idxvar],
                      cffhi[i],cfflo[i],acchi,acclo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffhi[i][L] = acchi[L];
                     cfflo[i][L] = acclo[L];
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
            cout << cffhi[i][j] << "  " << cfflo[i][j] << endl;
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
            cout << inputhi[i][j] << "  " << inputlo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl2_evaldiff
         (dim,nvr[i],deg,idx[i],cffhi[i],cfflo[i],
          inputhi,inputlo,outputhi[i],outputlo[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputhi[i][dim][j] << "  "
                 << outputlo[i][dim][j] << endl;
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
               double tmphi,tmplo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // output[i][idxvar][k] = factor*output[i][idxvar][k];
                  ddf_mul_d_dd(factor,outputhi[i][idxvar][k],
                                      outputlo[i][idxvar][k],&tmphi,&tmplo);
                  outputhi[i][idxvar][k] = tmphi;
                  outputlo[i][idxvar][k] = tmplo;
               }
            }
         }
      }
   }
}

void CPU_cmplx2_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehi, double **cffrelo, double **cffimhi, double **cffimlo,
   double *accrehi, double *accrelo, double *accimhi, double *accimlo,
   double **inputrehi, double **inputrelo,
   double **inputimhi, double **inputimlo, 
   double ***outputrehi, double ***outputrelo, 
   double ***outputimhi, double ***outputimlo, int vrblvl )
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
                  CPU_cmplx2_product
                     (deg,inputrehi[idxvar],inputrelo[idxvar],
                          inputimhi[idxvar],inputimlo[idxvar],
                      cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
                      accrehi,accrelo,accimhi,accimlo);

                  for(int L=0; L<=deg; L++)
                  {
                     cffrehi[i][L] = accrehi[L]; cffrelo[i][L] = accrelo[L];
                     cffimhi[i][L] = accimhi[L]; cffimlo[i][L] = accimlo[L];
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
            cout << cffrehi[i][j] << "  " << cffrelo[i][j] << endl << "  "
                 << cffimhi[i][j] << "  " << cffimlo[i][j] << endl;
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
            cout << inputrehi[i][j] << "  " << inputrelo[i][j] << "  "
                 << inputimhi[i][j] << "  " << inputimlo[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_cmplx2_evaldiff
         (dim,nvr[i],deg,idx[i],
          cffrehi[i],cffrelo[i],cffimhi[i],cffimlo[i],
          inputrehi,inputrelo,inputimhi,inputimlo,
          outputrehi[i],outputrelo[i],outputimhi[i],outputimlo[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputrehi[i][dim][j] << "  "
                 << outputrelo[i][dim][j] << endl << "  "
                 << outputimhi[i][dim][j] << "  "
                 << outputimlo[i][dim][j] << endl;
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
               double acchi,acclo;

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  // outputre[i][idxvar][k] = factor*outputre[i][idxvar][k];
                  ddf_mul_d_dd
                     (factor,outputrehi[i][idxvar][k],
                             outputrelo[i][idxvar][k],&acchi,&acclo);
                  outputrehi[i][idxvar][k] = acchi;
                  outputrelo[i][idxvar][k] = acclo;
                  // outputim[i][idxvar][k] = factor*outputim[i][idxvar][k];
                  ddf_mul_d_dd
                     (factor,outputimhi[i][idxvar][k],
                             outputimlo[i][idxvar][k],&acchi,&acclo);
                  outputimhi[i][idxvar][k] = acchi;
                  outputimlo[i][idxvar][k] = acclo;
               }
            }
         }
      }
   }
}

void dbl2_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbhi, double **mblo, double damper,
   double ***outputhi, double ***outputlo,
   double **funvalhi, double **funvallo, 
   double **rhshi, double **rhslo, double ***jacvalhi, double ***jacvallo,
   int vrblvl )
{
   double acchi,acclo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalhi[i][j] = outputhi[i][dim][j];
         funvallo[i][j] = outputlo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalhi[i][0] << "  "
                            << funvallo[i][0] << endl;
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
         ddf_sub(funvalhi[j][i],funvallo[j][i],
                 mbhi[j][i],mblo[j][i],&acchi,&acclo);
         rhshi[i][j] = -acchi;
         rhslo[i][j] = -acclo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhshi[i][j] << "  " << rhslo[i][j] << endl;
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
            jacvalhi[j][i][idxval] = outputhi[i][idxval][j];
            jacvallo[j][i][idxval] = outputlo[i][idxval][j];
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
            cout << jacvalhi[0][i][j] << "  " << jacvallo[0][i][j] << endl;
      }
   }
}

void cmplx2_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbrehi, double **mbrelo, double **mbimhi, double **mbimlo,
   double damper,
   double ***outputrehi, double ***outputrelo,
   double ***outputimhi, double ***outputimlo,
   double **funvalrehi, double **funvalrelo,
   double **funvalimhi, double **funvalimlo,
   double **rhsrehi, double **rhsrelo, double **rhsimhi, double **rhsimlo,
   double ***jacvalrehi, double ***jacvalrelo,
   double ***jacvalimhi, double ***jacvalimlo, int vrblvl )
{
   double acchi,acclo;

   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalrehi[i][j] = outputrehi[i][dim][j];
         funvalrelo[i][j] = outputrelo[i][dim][j];
         funvalimhi[i][j] = outputimhi[i][dim][j];
         funvalimlo[i][j] = outputimlo[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalrehi[i][0] << "  " << funvalrelo[i][0] << endl << "  "
              << funvalimhi[i][0] << "  " << funvalimlo[i][0] << endl;
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
         ddf_sub(funvalrehi[j][i],funvalrelo[j][i],
                 mbrehi[j][i],mbrelo[j][i],&acchi,&acclo);
         rhsrehi[i][j] = -acchi;
         rhsrelo[i][j] = -acclo;
         // rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
         ddf_sub(funvalimhi[j][i],funvalimlo[j][i],
                 mbimhi[j][i],mbimlo[j][i],&acchi,&acclo);
         rhsimhi[i][j] = -acchi;
         rhsimlo[i][j] = -acclo;
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsrehi[i][j] << "  " << rhsrelo[i][j] << endl << "  "
                 << rhsimhi[i][j] << "  " << rhsimlo[i][j] << endl;
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
            jacvalrehi[j][i][idxval] = outputrehi[i][idxval][j];
            jacvalrelo[j][i][idxval] = outputrelo[i][idxval][j];
            jacvalimhi[j][i][idxval] = outputimhi[i][idxval][j];
            jacvalimlo[j][i][idxval] = outputimlo[i][idxval][j];
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
            cout << jacvalrehi[0][i][j] << "  "
                 << jacvalrelo[0][i][j] << endl << "  "
                 << jacvalimhi[0][i][j] << "  "
                 << jacvalimlo[0][i][j] << endl;
      }
   }
}
