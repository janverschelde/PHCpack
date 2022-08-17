// The file dbl4_newton_testers.cpp defines the functions with prototypes in
// the file dbl4_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "quad_double_functions.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"
#include "dbl4_factorizations.h"
#include "dbl4_bals_host.h"

using namespace std;

void dbl4_unit_series_vector
 ( int dim, int deg,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihi[i][0] = 1.0;
      cfflohi[i][0] = 0.0;
      cffhilo[i][0] = 0.0;
      cfflolo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhihi[i][j] = 0.0;
         cfflohi[i][j] = 0.0;
         cffhilo[i][j] = 0.0;
         cfflolo[i][j] = 0.0;
      }
   }
}

void dbl4_evaluate_monomials
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

void dbl4_update_series
 ( int dim, int degp1,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **dxhihi, double **dxlohi, double **dxhilo, double **dxlolo,
   int vrblvl )
{
   if(vrblvl > 0)
   {
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << xhihi[i][j] << "  " << xlohi[i][j] << endl;
            cout << xhilo[i][j] << "  " << xlolo[i][j] << endl;
         }
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         qdf_inc(&xhihi[i][j],&xlohi[i][j],&xhilo[i][j],&xlolo[i][j],
                 dxhihi[j][i],dxlohi[j][i],dxhilo[j][i],dxlolo[j][i]);
      }

   if(vrblvl > 0)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << xhihi[i][j] << "  " << xlohi[i][j] << endl;
            cout << xhilo[i][j] << "  " << xlolo[i][j] << endl;
         }
      }
   }
}

void dbl4_newton_step
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi, double **inputlohi,
   double **inputhilo, double **inputlolo,
   double ***outputhihi, double ***outputlohi,
   double ***outputhilo, double ***outputlolo,
   double **funvalhihi, double **funvallohi,
   double **funvalhilo, double **funvallolo,
   double ***jacvalhihi, double ***jacvallohi,
   double ***jacvalhilo, double ***jacvallolo,
   double **rhshihi, double **rhslohi, double **rhshilo, double **rhslolo,
   double **solhihi, double **sollohi, double **solhilo, double **sollolo,
   double **workmathihi, double **workmatlohi,
   double **workmathilo, double **workmatlolo,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **workrhshihi, double **workrhslohi,
   double **workrhshilo, double **workrhslolo,
   double **resvechihi, double **resveclohi,
   double **resvechilo, double **resveclolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int *ipvt, int vrblvl )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   for(int i=0; i<dim; i++)
   {
      cffhihi[i][0] = 1.0;
      cfflohi[i][0] = 0.0;
      cffhilo[i][0] = 0.0;
      cfflolo[i][0] = 0.0;

      for(int j=1; j<degp1; j++)
      {
         cffhihi[i][j] = 0.0;
         cfflohi[i][j] = 0.0;
         cffhilo[i][j] = 0.0;
         cfflolo[i][j] = 0.0;
      }
   }
   dbl4_evaluate_monomials
      (dim,deg,nvr,idx,exp,nbrfac,expfac,
       cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
        inputhihi, inputlohi, inputhilo, inputlolo,
       outputhihi,outputlohi,outputhilo,outputlolo,vrblvl);

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhihi[i][j][k] = 0.0;
            jacvallohi[i][j][k] = 0.0;
            jacvalhilo[i][j][k] = 0.0;
            jacvallolo[i][j][k] = 0.0;
         }

   dbl4_linearize_evaldiff_output
      (dim,degp1,nvr,idx,outputhihi,outputlohi,outputhilo,outputlolo,
       funvalhihi,funvallohi,funvalhilo,funvallolo,
       rhshihi,rhslohi,rhshilo,rhslolo,
       jacvalhihi,jacvallohi,jacvalhilo,jacvallolo,vrblvl);

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         workrhshihi[i][j] = rhshihi[i][j];
         workrhslohi[i][j] = rhslohi[i][j];
         workrhshilo[i][j] = rhshilo[i][j];
         workrhslolo[i][j] = rhslolo[i][j];
      }
   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhihi[i][j] = 0.0;
         sollohi[i][j] = 0.0;
         solhilo[i][j] = 0.0;
         sollolo[i][j] = 0.0;
      }
 
   CPU_dbl4_linear_solve
      (dim,degp1,jacvalhihi,jacvallohi,jacvalhilo,jacvallolo,
       workrhshihi,workrhslohi,workrhshilo,workrhslolo,
       solhihi,sollohi,solhilo,sollolo,
       workmathihi,workmatlohi,workmathilo,workmatlolo,
       workvechihi,workveclohi,workvechilo,workveclolo,ipvt,0); // vrblvl);

   CPU_dbl4_linear_residue
      (dim,degp1,jacvalhihi,jacvallohi,jacvalhilo,jacvallolo,
       rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
       resvechihi,resveclohi,resvechilo,resveclolo,
       resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);

   if(vrblvl > 0)
      cout << "maximum residual : "
           << *resmaxhihi << "  " << *resmaxlohi << endl
           << *resmaxhilo << "  " << *resmaxlolo << endl;

   dbl4_update_series
      (dim,degp1,inputhihi,inputlohi,inputhilo,inputlolo,
       solhihi,sollohi,solhilo,sollolo,vrblvl);
}
