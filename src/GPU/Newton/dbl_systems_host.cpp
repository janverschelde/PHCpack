// The file dbl_systems_host.cpp defines the functions with prototypes in
// the file dbl_systems_host.h.

#include <iostream>
#include <iomanip>
#include "dbl_convolutions_host.h"
#include "dbl_monomials_host.h"

using namespace std;

void CPU_dbl_evaluate_monomials
 ( int dim, int deg, int *nvr, int **idx, int **exp, int *nbrfac,
   int **expfac, double **cff, double *acc, double **input,
   double ***output, int vrblvl )
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
                  CPU_dbl_product(deg,input[idxvar],cff[i],acc);
                  for(int L=0; L<=deg; L++) cff[i][L] = acc[L];
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
         for(int j=0; j<=deg; j++) cout << cff[i][j] << endl;
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
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_dbl_evaldiff(dim,nvr[i],deg,idx[i],cff[i],input,output[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++) cout << output[i][dim][j] << endl;
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

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
                  output[i][idxvar][k] = factor*output[i][idxvar][k];
            }
         }
      }
   }
}

void CPU_cmplx_evaluate_monomials
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffre, double **cffim, double *accre, double *accim,
   double **inputre, double **inputim, double ***outputre, double ***outputim,
   int vrblvl )
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
                  CPU_cmplx_product
                     (deg,inputre[idxvar],inputim[idxvar],
                      cffre[i],cffim[i],accre,accim);

                  for(int L=0; L<=deg; L++)
                  {
                     cffre[i][L] = accre[L];
                     cffim[i][L] = accim[L];
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
            cout << cffre[i][j] << "  " << cffim[i][j] << endl;
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
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
      }
   }
   for(int i=0; i<dim; i++)
      CPU_cmplx_evaldiff
         (dim,nvr[i],deg,idx[i],cffre[i],cffim[i],
          inputre,inputim,outputre[i],outputim[i]);

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for monomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << outputre[i][dim][j] << "  "
                 << outputim[i][dim][j] << endl;
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

               // multiply derivative w.r.t. idxvar with factor
               for(int k=0; k<=deg; k++)
               {
                  outputre[i][idxvar][k] = factor*outputre[i][idxvar][k];
                  outputim[i][idxvar][k] = factor*outputim[i][idxvar][k];
               }
            }
         }
      }
   }
}

void CPU_dbl_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cff, double **acc, double **input,
   double **funval, double ***jacval, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++) funval[i][j] = 0.0;
   funval[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++) jacval[k][i][j] = 0.0;

   if(vrblvl > 1)
   {
      for(int i=0; i<nbrcol; i++)
      {
         cout << "coefficients for column " << i << " :" << endl;
         for(int j=0; j<dim; j++)
         {
            cout << "coefficients for monomial " << j << " :" << endl;
            for(int k=0; k<=deg; k++) cout << cff[i][j][k] << endl;
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
         for(int j=0; j<=deg; j++) cout << input[i][j] << endl;
      }
   }
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // evaluate j-th monomial in column i
         {
            for(int k=0; k<=dim; k++)
               for(int L=0; L<degp1; L++) acc[k][L] = 0.0;

            CPU_dbl_evaldiff(dim,nvr[i][j],deg,idx[i][j],cff[i][j],input,acc);

            for(int L=0; L<degp1; L++) funval[j][L] += acc[dim][L];

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
                  jacval[L][j][idxval] += acc[idxval][L];
            }
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<degp1; j++) cout << funval[i][j] << endl;
      }
   }
}

void CPU_cmplx_evaluate_columns
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx,
   double ***cffre, double ***cffim, double **accre, double **accim,
   double **inputre, double **inputim, double **funvalre, double **funvalim,
   double ***jacvalre, double ***jacvalim, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalre[i][j] = 0.0;
         funvalim[i][j] = 0.0;
      }
   funvalre[dim-1][0] = -1.0; // constant of last eq in cyclic dim-roots

   for(int k=0; k<degp1; k++)  // the Jacobian is linearized
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            jacvalre[k][i][j] = 0.0;
            jacvalim[k][i][j] = 0.0;
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
               cout << cffre[i][j][k] << "  " << cffim[i][j][k] << endl;
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
            cout << inputre[i][j] << "  " << inputim[i][j] << endl;
      }
   }
   for(int i=0; i<nbrcol; i++)
      for(int j=0; j<dim; j++)
         if(nvr[i][j] > 0)       // evaluate j-th monomial in column i
         {
            CPU_cmplx_evaldiff
               (dim,nvr[i][j],deg,idx[i][j],cffre[i][j],cffim[i][j],
                inputre,inputim,accre,accim);

            for(int L=0; L<degp1; L++)
            {
               funvalre[j][L] += accre[dim][L];
               funvalim[j][L] += accim[dim][L];
            }

            int *indexes = idx[i][j];      // indices of the variables
            for(int k=0; k<nvr[i][j]; k++) // derivative w.r.t. idx[i][j][k]
            {                              // has j-th coefficient
               int idxval = indexes[k];
               for(int L=0; L<degp1; L++) 
               {
                  jacvalre[L][j][idxval] += accre[idxval][L];
                  jacvalim[L][j][idxval] += accim[idxval][L];
               }
            }
         }

   if(vrblvl > 1)
   {
      for(int i=0; i<dim; i++)
      {
         cout << "output series for polynomial " << i << " :" << endl;
         for(int j=0; j<=deg; j++)
            cout << funvalre[i][j] << "  "
                 << funvalim[i][j] << endl;
      }
   }
}

void dbl_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx, double **mb, double damper,
   double ***output, double **funval, double **rhs, double ***jacval,
   int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++) funval[i][j] = output[i][dim][j];

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++) cout << i << " : " << funval[i][0] << endl;
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
         rhs[i][j] = -(funval[j][i] - mb[j][i]);
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << rhs[i][j] << endl;
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
            jacval[j][i][idxval] = output[i][idxval][j];
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacval[0][i][j] << endl;
      }
   }
}

void cmplx_linearize_evaldiff_output
 ( int dim, int degp1, int *nvr, int **idx,
   double **mbre, double **mbim, double damper,
   double ***outputre, double ***outputim,
   double **funvalre, double **funvalim,
   double **rhsre, double **rhsim, double ***jacvalre, double ***jacvalim,
   int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         funvalre[i][j] = outputre[i][dim][j];
         funvalim[i][j] = outputim[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << funvalre[i][0] << "  " << funvalim[i][0] << endl;
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
         rhsre[i][j] = -(funvalre[j][i] - mbre[j][i]);
         rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
      }

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsre[i][j] << "  " << rhsim[i][j] << endl;
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
            jacvalre[j][i][idxval] = outputre[i][idxval][j];
            jacvalim[j][i][idxval] = outputim[i][idxval][j];
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
            cout << jacvalre[0][i][j] << "  " << jacvalim[0][i][j] << endl;
      }
   }
}

void dbl_define_rhs
 ( int dim, int degp1, double **mb, double **funval, double **rhs,
   int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++) cout << i << " : " << funval[i][0] << endl;
   }
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++) rhs[i][j] = -(funval[j][i] - mb[j][i]);

   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++) cout << rhs[i][j] << endl;
      }
   }
}

void cmplx_define_rhs
 ( int dim, int degp1, double **mbre, double **mbim,
   double **funvalre, double **funvalim,
   double **rhsre, double **rhsim, int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the evaluated series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " 
              << funvalre[i][0] << "  " << funvalim[i][0] << endl;
   }
   for(int i=0; i<degp1; i++)
      for(int j=0; j<dim; j++)
      {
         rhsre[i][j] = -(funvalre[j][i] - mbre[j][i]);
         rhsim[i][j] = -(funvalim[j][i] - mbim[j][i]);
      }
 
   if(vrblvl > 1)
   {
      cout << "The right hand side series :" << endl;
      for(int i=0; i<degp1; i++)
      {
         cout << "coefficient vector " << i << " :" << endl;
         for(int j=0; j<dim; j++)
            cout << rhsre[i][j] << "  " << rhsim[i][j] << endl;
      }
   }
}

void dbl_map_evaldiff_output
 ( int dim, int deg, double ***output, double **funval, double ***jacval,
   int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
         funval[i][j] = output[i][dim][j];

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the evaluated series :" << endl;

      for(int i=0; i<dim; i++)
         cout << i << " : " << funval[i][0] << endl;
   }
   // output[i][j][k] is the k-th coefficient in the series
   // the derivative of the i-th polynomial with respect to variable j.

   for(int i=0; i<dim; i++)          // the i-th polynomial
   {
      for(int j=0; j<dim; j++)       // derivative w.r.t. j-th variable
      {
         for(int k=0; k<=deg; k++) 
            jacval[k][i][j] = output[i][j][k];
      }
   }
   if(vrblvl > 1)
   {
      cout << "The leading coefficients of the Jacobian matrix : " << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "row " << i << " : " << endl;
         for(int j=0; j<dim; j++)
            cout << jacval[0][i][j] << endl;
      }
   }
}

void cmplx_map_evaldiff_output
 ( int dim, int deg, double ***outputre, double ***outputim,
   double **funvalre, double **funvalim,
   double ***jacvalre, double ***jacvalim, int vrblvl )
{
   // The coefficients of the series for the function values
   // at the input are in output[i][dim].
   for(int i=0; i<dim; i++)
      for(int j=0; j<=deg; j++)
      {
         funvalre[i][j] = outputre[i][dim][j];
         funvalim[i][j] = outputim[i][dim][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the evaluated series :" << endl;

      for(int i=0; i<dim; i++)
         cout << i << " : " << funvalre[i][0] << "  "
                            << funvalim[i][0] << endl;
   }
   // output[i][j][k] is the k-th coefficient in the series
   // the derivative of the i-th polynomial with respect to variable j.

   for(int i=0; i<dim; i++)          // the i-th polynomial
   {
      for(int j=0; j<dim; j++)       // derivative w.r.t. j-th variable
      {
         for(int k=0; k<=deg; k++) 
         {
            jacvalre[k][i][j] = outputre[i][j][k];
            jacvalim[k][i][j] = outputim[i][j][k];
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
            cout << jacvalre[0][i][j] << "  "
                 << jacvalim[0][i][j] << endl;
      }
   }
}
