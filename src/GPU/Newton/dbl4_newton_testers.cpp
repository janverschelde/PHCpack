// The file dbl4_newton_testers.cpp defines the functions with prototypes in
// the file dbl4_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "quad_double_functions.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"
#include "dbl4_factorizations.h"
#include "dbl4_bals_host.h"
#include "dbl4_bals_kernels.h"
#include "dbl4_systems_host.h"
#include "dbl4_systems_kernels.h"

using namespace std;

void dbl4_unit_series_vector
 ( int dim, int deg,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihi[i][0] = 1.000005;
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

void cmplx4_unit_series_vector
 ( int dim, int deg,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo )
{
   for(int i=0; i<dim; i++)
   {
      cffrehihi[i][0] = 1.000005; cffrelohi[i][0] = 0.0;
      cffrehilo[i][0] = 0.0; cffrelolo[i][0] = 0.0;
      cffimhihi[i][0] = 0.000005; cffimlohi[i][0] = 0.0;
      cffimhilo[i][0] = 0.0; cffimlolo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffrehihi[i][j] = 0.0; cffrelohi[i][j] = 0.0;
         cffrehilo[i][j] = 0.0; cffrelolo[i][j] = 0.0;
         cffimhihi[i][j] = 0.0; cffimlohi[i][j] = 0.0;
         cffimhilo[i][j] = 0.0; cffimlolo[i][j] = 0.0;
      }
   }
}

void dbl4_update_series
 ( int dim, int degp1,
   double **xhihi, double **xlohi, double **xhilo, double **xlolo,
   double **dxhihi, double **dxlohi, double **dxhilo, double **dxlolo,
   int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
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
      cout << scientific << setprecision(16);
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

void cmplx4_update_series
 ( int dim, int degp1,
   double **xrehihi, double **xrelohi, double **xrehilo, double **xrelolo,
   double **ximhihi, double **ximlohi, double **ximhilo, double **ximlolo,
   double **dxrehihi, double **dxrelohi, double **dxrehilo, double **dxrelolo,
   double **dximhihi, double **dximlohi, double **dximhilo, double **dximlolo,
   int vrblvl )
{
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The series before the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehihi[i][j] << "  " << xrelohi[i][j] << endl << "  "
                 << xrehilo[i][j] << "  " << xrelolo[i][j] << endl << "  "
                 << ximhihi[i][j] << "  " << ximlohi[i][j] << endl << "  "
                 << ximhilo[i][j] << "  " << ximlolo[i][j] << endl;
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++)
      {
         // xre[i][j] = xre[i][j] + dxre[j][i];
         qdf_inc(&xrehihi[i][j],&xrelohi[i][j],&xrehilo[i][j],&xrelolo[i][j],
                 dxrehihi[j][i],dxrelohi[j][i],dxrehilo[j][i],dxrelolo[j][i]);
         // xim[i][j] = xim[i][j] + dxim[j][i];
         qdf_inc(&ximhihi[i][j],&ximlohi[i][j],&ximhilo[i][j],&ximlolo[i][j],
                 dximhihi[j][i],dximlohi[j][i],dximhilo[j][i],dximlolo[j][i]);
      }
 
   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
            cout << xrehihi[i][j] << "  " << xrelohi[i][j] << endl << "  "
                 << xrehilo[i][j] << "  " << xrelolo[i][j] << endl << "  "
                 << ximhihi[i][j] << "  " << ximlohi[i][j] << endl << "  "
                 << ximhilo[i][j] << "  " << ximlolo[i][j] << endl;
      }
   }
}

void dbl4_newton_lustep
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
   CPU_dbl4_evaluate_monomials
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
 
   CPU_dbl4_lusb_solve
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

   if(vrblvl > 0) cout << "maximum residual : " << *resmaxhihi << endl;

   dbl4_update_series
      (dim,degp1,inputhihi,inputlohi,inputhilo,inputlolo,
       solhihi,sollohi,solhilo,sollolo,vrblvl);
}

void dbl4_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihi, double **cfflohi, double **cffhilo, double **cfflolo,
   double *acchihi, double *acclohi, double *acchilo, double *acclolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h,
   double ***outputhilo_h, double ***outputlolo_h,
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d,
   double **funvalhihi_h, double **funvallohi_h,
   double **funvalhilo_h, double **funvallolo_h,
   double **funvalhihi_d, double **funvallohi_d,
   double **funvalhilo_d, double **funvallolo_d,
   double ***jacvalhihi_h, double ***jacvallohi_h,
   double ***jacvalhilo_h, double ***jacvallolo_h,
   double ***jacvalhihi_d, double ***jacvallohi_d,
   double ***jacvalhilo_d, double ***jacvallolo_d,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double **workmathihi, double **workmatlohi,
   double **workmathilo, double **workmatlolo,
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl, int mode )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   if((mode == 1) || (mode == 2))
   {
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
      CPU_dbl4_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
          inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
          outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
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
      GPU_dbl4_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;
      for(int k=0; k<dim; k++) // monomial k
         for(int i=0; i<=dim; i++)
            for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "output_h[" << k << "][" << i << "][" << j << "] : "
                    << outputhihi_h[k][i][j] << "  "
                    << outputlohi_h[k][i][j] << endl << "  "
                    << outputhilo_h[k][i][j] << "  "
                    << outputlolo_h[k][i][j] << endl;
               cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                    << outputhihi_d[k][i][j] << "  "
                    << outputlohi_d[k][i][j] << endl << "  "
                    << outputhilo_d[k][i][j] << "  "
                    << outputlolo_d[k][i][j] << endl;
            }
            errsum += abs(outputhihi_h[k][i][j] - outputhihi_d[k][i][j])
                    + abs(outputlohi_h[k][i][j] - outputlohi_d[k][i][j])
                    + abs(outputhilo_h[k][i][j] - outputhilo_d[k][i][j])
                    + abs(outputlolo_h[k][i][j] - outputlolo_d[k][i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }

   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhihi_h[i][j][k] = 0.0;
            jacvallohi_h[i][j][k] = 0.0;
            jacvalhilo_h[i][j][k] = 0.0;
            jacvallolo_h[i][j][k] = 0.0;
            jacvalhihi_d[i][j][k] = 0.0;
            jacvallohi_d[i][j][k] = 0.0;
            jacvalhilo_d[i][j][k] = 0.0;
            jacvallolo_d[i][j][k] = 0.0;
         }

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
      dbl4_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
          funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
          rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
          jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,vrblvl);
   if((mode == 1) || (mode == 2))
      dbl4_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
          funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
          rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
          jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,vrblvl);

   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "funval_h[" << i << "][" << j << "] : "
                    << funvalhihi_h[i][j] << "  "
                    << funvallohi_h[i][j] << endl << "  "
                    << funvalhilo_h[i][j] << "  "
                    << funvallolo_h[i][j] << endl;
               cout << "funval_d[" << i << "][" << j << "] : "
                    << funvalhihi_d[i][j] << "  "
                    << funvallohi_d[i][j] << endl << "  "
                    << funvalhilo_d[i][j] << "  "
                    << funvallolo_d[i][j] << endl;
            }
            errsum += abs(funvalhihi_h[i][j] - funvalhihi_d[i][j])
                    + abs(funvallohi_h[i][j] - funvallohi_d[i][j])
                    + abs(funvalhilo_h[i][j] - funvalhilo_d[i][j])
                    + abs(funvallolo_h[i][j] - funvallolo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            for(int k=0; k<dim; k++)
            {
               if(vrblvl > 1)
               {
                  cout << "jacval_h[" << i << "][" << j << "][" << k << "] : "
                       << jacvalhihi_h[i][j][k] << "  "
                       << jacvallohi_h[i][j][k] << endl << "  "
                       << jacvalhilo_h[i][j][k] << "  "
                       << jacvallolo_h[i][j][k] << endl;
                  cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                       << jacvalhihi_d[i][j][k] << "  "
                       << jacvallohi_d[i][j][k] << endl << "  "
                       << jacvalhilo_d[i][j][k] << "  "
                       << jacvallolo_d[i][j][k] << endl;
               }
               errsum += abs(jacvalhihi_h[i][j][k] - jacvalhihi_d[i][j][k])
                       + abs(jacvallohi_h[i][j][k] - jacvallohi_d[i][j][k])
                       + abs(jacvalhilo_h[i][j][k] - jacvalhilo_d[i][j][k])
                       + abs(jacvallolo_h[i][j][k] - jacvallolo_d[i][j][k]);
            }
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "rhs_h[" << i << "][" << j << "] : "
                    << rhshihi_h[i][j] << "  " << rhslohi_h[i][j] << endl
                    << "  "
                    << rhshilo_h[i][j] << "  " << rhslolo_h[i][j] << endl;
               cout << "rhs_d[" << i << "][" << j << "] : "
                    << rhshihi_d[i][j] << "  " << rhslohi_d[i][j] << endl
                    << "  "
                    << rhshilo_d[i][j] << "  " << rhslolo_d[i][j] << endl;
            }
            errsum += abs(rhshihi_h[i][j] - rhshihi_d[i][j])
                    + abs(rhslohi_h[i][j] - rhslohi_d[i][j])
                    + abs(rhshilo_h[i][j] - rhshilo_d[i][j])
                    + abs(rhslolo_h[i][j] - rhslolo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhshihi_h[i][j] = rhshihi_h[i][j];
         urhslohi_h[i][j] = rhslohi_h[i][j];
         urhshilo_h[i][j] = rhshilo_h[i][j];
         urhslolo_h[i][j] = rhslolo_h[i][j];
         urhshihi_d[i][j] = rhshihi_d[i][j];
         urhslohi_d[i][j] = rhslohi_d[i][j];
         urhshilo_d[i][j] = rhshilo_d[i][j];
         urhslolo_d[i][j] = rhslolo_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhihi_h[i][j] = 0.0; sollohi_h[i][j] = 0.0;
         solhilo_h[i][j] = 0.0; sollolo_h[i][j] = 0.0;
         solhihi_d[i][j] = 0.0; sollohi_d[i][j] = 0.0;
         solhilo_d[i][j] = 0.0; sollolo_d[i][j] = 0.0;
      }
 
   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling CPU_dbl4_qrbs_solve ..." << endl;
      CPU_dbl4_qrbs_solve
         (dim,degp1,jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
          urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
          solhihi_h,sollohi_h,solhilo_h,sollolo_h,
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
          workvechihi,workveclohi,workvechilo,workveclolo,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl4_linear_residue ..." << endl;

         CPU_dbl4_linear_residue
            (dim,degp1,jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
             rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
             solhihi_h,sollohi_h,solhilo_h,sollolo_h,
             resvechihi,resveclohi,resvechilo,resveclolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);

         cout << "maximum residual : " << *resmaxhihi << endl;
      }
      dbl4_update_series
         (dim,degp1,inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
          solhihi_h,sollohi_h,solhilo_h,sollolo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling GPU_dbl4_bals_solve ..." << endl;

      GPU_dbl4_bals_solve
         (dim,degp1,szt,nbt,
          jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
          Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
          urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
          solhihi_d,sollohi_d,solhilo_d,sollolo_d,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl4_linear_residue ..." << endl;

         CPU_dbl4_linear_residue
            (dim,degp1,jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
             rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
             solhihi_d,sollohi_d,solhilo_d,sollolo_d,
             resvechihi,resveclohi,resvechilo,resveclolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);

         cout << "maximum residual : " << *resmaxhihi << "  " << endl;
      }
      dbl4_update_series
         (dim,degp1,inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          solhihi_d,sollohi_d,solhilo_d,sollolo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "Q_h[" << i << "][" << j << "] : "
                    << Qhihi_h[i][j] << "  " << Qlohi_h[i][j] << endl << "  "
                    << Qhilo_h[i][j] << "  " << Qlolo_h[i][j] << endl;
               cout << "Q_d[" << i << "][" << j << "] : "
                    << Qhihi_d[i][j] << "  " << Qlohi_d[i][j] << endl << "  "
                    << Qhilo_d[i][j] << "  " << Qlolo_d[i][j] << endl;
            }
            errsum += abs(Qhihi_h[i][j] - Qhihi_d[i][j])
                    + abs(Qlohi_h[i][j] - Qlohi_d[i][j])
                    + abs(Qhilo_h[i][j] - Qhilo_d[i][j])
                    + abs(Qlolo_h[i][j] - Qlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "R_h[" << i << "][" << j << "] : "
                    << Rhihi_h[i][j] << "  " << Rlohi_h[i][j] << endl << "  "
                    << Rhilo_h[i][j] << "  " << Rlolo_h[i][j] << endl;
               cout << "R_d[" << i << "][" << j << "] : "
                    << Rhihi_d[i][j] << "  " << Rlohi_d[i][j] << endl << "  "
                    << Rhilo_d[i][j] << "  " << Rlolo_d[i][j] << endl;
            }
            errsum += abs(Rhihi_h[i][j] - Rhihi_d[i][j])
                    + abs(Rlohi_h[i][j] - Rlohi_d[i][j])
                    + abs(Rhilo_h[i][j] - Rhilo_d[i][j])
                    + abs(Rlolo_h[i][j] - Rlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "urhs_h[" << i << "][" << j << "] : "
                    << urhshihi_h[i][j] << "  " << urhslohi_h[i][j]
                    << endl << "  "
                    << urhshilo_h[i][j] << "  " << urhslolo_h[i][j] << endl;
               cout << "urhs_d[" << i << "][" << j << "] : "
                    << urhshihi_d[i][j] << "  " << urhslohi_d[i][j]
                    << endl << "  "
                    << urhshilo_d[i][j] << "  " << urhslolo_d[i][j] << endl;
            }
            errsum += abs(urhshihi_h[i][j] - urhshihi_d[i][j])
                    + abs(urhslohi_h[i][j] - urhslohi_d[i][j])
                    + abs(urhshilo_h[i][j] - urhshilo_d[i][j])
                    + abs(urhslolo_h[i][j] - urhslolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "sol_h[" << i << "][" << j << "] : "
                    << solhihi_h[i][j] << "  " << sollohi_h[i][j] << endl
                    << "  "
                    << solhilo_h[i][j] << "  " << sollolo_h[i][j] << endl;
               cout << "sol_d[" << i << "][" << j << "] : "
                    << solhihi_d[i][j] << "  " << sollohi_d[i][j] << endl
                    << "  "
                    << solhilo_d[i][j] << "  " << sollolo_d[i][j] << endl;
            }
            errsum += abs(solhihi_h[i][j] - solhihi_d[i][j])
                    + abs(sollohi_h[i][j] - sollohi_d[i][j])
                    + abs(solhilo_h[i][j] - solhilo_d[i][j])
                    + abs(sollolo_h[i][j] - sollolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "input_h[" << i << "][" << j << "] : "
                    << inputhihi_h[i][j] << "  "
                    << inputlohi_h[i][j] << endl << "  "
                    << inputhilo_h[i][j] << "  "
                    << inputlolo_h[i][j] << endl;
               cout << "input_d[" << i << "][" << j << "] : "
                    << inputhihi_d[i][j] << "  "
                    << inputlohi_d[i][j] << endl << "  "
                    << inputhilo_d[i][j] << "  "
                    << inputlolo_d[i][j] << endl;
            }
            errsum += abs(inputhihi_h[i][j] - inputhihi_d[i][j])
                    + abs(inputlohi_h[i][j] - inputlohi_d[i][j])
                    + abs(inputhilo_h[i][j] - inputhilo_d[i][j])
                    + abs(inputlolo_h[i][j] - inputlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}

void cmplx4_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffrehihi, double **cffrelohi,
   double **cffrehilo, double **cffrelolo,
   double **cffimhihi, double **cffimlohi,
   double **cffimhilo, double **cffimlolo,
   double *accrehihi, double *accrelohi,
   double *accrehilo, double *accrelolo,
   double *accimhihi, double *accimlohi,
   double *accimhilo, double *accimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d,
   double ***outputrehihi_h, double ***outputrelohi_h,
   double ***outputrehilo_h, double ***outputrelolo_h,
   double ***outputimhihi_h, double ***outputimlohi_h,
   double ***outputimhilo_h, double ***outputimlolo_h,
   double ***outputrehihi_d, double ***outputrelohi_d,
   double ***outputrehilo_d, double ***outputrelolo_d,
   double ***outputimhihi_d, double ***outputimlohi_d,
   double ***outputimhilo_d, double ***outputimlolo_d,
   double **funvalrehihi_h, double **funvalrelohi_h,
   double **funvalrehilo_h, double **funvalrelolo_h,
   double **funvalimhihi_h, double **funvalimlohi_h,
   double **funvalimhilo_h, double **funvalimlolo_h,
   double **funvalrehihi_d, double **funvalrelohi_d,
   double **funvalrehilo_d, double **funvalrelolo_d,
   double **funvalimhihi_d, double **funvalimlohi_d,
   double **funvalimhilo_d, double **funvalimlolo_d,
   double ***jacvalrehihi_h, double ***jacvalrelohi_h,
   double ***jacvalrehilo_h, double ***jacvalrelolo_h,
   double ***jacvalimhihi_h, double ***jacvalimlohi_h,
   double ***jacvalimhilo_h, double ***jacvalimlolo_h,
   double ***jacvalrehihi_d, double ***jacvalrelohi_d,
   double ***jacvalrehilo_d, double ***jacvalrelolo_d,
   double ***jacvalimhihi_d, double ***jacvalimlohi_d,
   double ***jacvalimhilo_d, double ***jacvalimlolo_d,
   double **rhsrehihi_h, double **rhsrelohi_h,
   double **rhsrehilo_h, double **rhsrelolo_h,
   double **rhsimhihi_h, double **rhsimlohi_h,
   double **rhsimhilo_h, double **rhsimlolo_h,
   double **rhsrehihi_d, double **rhsrelohi_d, 
   double **rhsrehilo_d, double **rhsrelolo_d, 
   double **rhsimhihi_d, double **rhsimlohi_d,
   double **rhsimhilo_d, double **rhsimlolo_d,
   double **urhsrehihi_h, double **urhsrelohi_h,
   double **urhsrehilo_h, double **urhsrelolo_h,
   double **urhsimhihi_h, double **urhsimlohi_h,
   double **urhsimhilo_h, double **urhsimlolo_h,
   double **urhsrehihi_d, double **urhsrelohi_d,
   double **urhsrehilo_d, double **urhsrelolo_d,
   double **urhsimhihi_d, double **urhsimlohi_d,
   double **urhsimhilo_d, double **urhsimlolo_d,
   double **solrehihi_h, double **solrelohi_h,
   double **solrehilo_h, double **solrelolo_h,
   double **solimhihi_h, double **solimlohi_h, 
   double **solimhilo_h, double **solimlolo_h, 
   double **solrehihi_d, double **solrelohi_d, 
   double **solrehilo_d, double **solrelolo_d, 
   double **solimhihi_d, double **solimlohi_d,
   double **solimhilo_d, double **solimlolo_d,
   double **Qrehihi_h, double **Qrelohi_h,
   double **Qrehilo_h, double **Qrelolo_h,
   double **Qimhihi_h, double **Qimlohi_h,
   double **Qimhilo_h, double **Qimlolo_h,
   double **Qrehihi_d, double **Qrelohi_d,
   double **Qrehilo_d, double **Qrelolo_d,
   double **Qimhihi_d, double **Qimlohi_d, 
   double **Qimhilo_d, double **Qimlolo_d, 
   double **Rrehihi_h, double **Rrelohi_h,
   double **Rrehilo_h, double **Rrelolo_h,
   double **Rimhihi_h, double **Rimlohi_h, 
   double **Rimhilo_h, double **Rimlolo_h, 
   double **Rrehihi_d, double **Rrelohi_d,
   double **Rrehilo_d, double **Rrelolo_d,
   double **Rimhihi_d, double **Rimlohi_d,
   double **Rimhilo_d, double **Rimlolo_d,
   double **workmatrehihi, double **workmatrelohi,
   double **workmatrehilo, double **workmatrelolo,
   double **workmatimhihi, double **workmatimlohi,
   double **workmatimhilo, double **workmatimlolo,
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo, int vrblvl, int mode )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         cffrehihi[i][0] = 1.0; cffrelohi[i][0] = 0.0;
         cffrehilo[i][0] = 0.0; cffrelolo[i][0] = 0.0;
         cffimhihi[i][0] = 0.0; cffimlohi[i][0] = 0.0;
         cffimhilo[i][0] = 0.0; cffimlolo[i][0] = 0.0;

         for(int j=1; j<degp1; j++)
         {
            cffrehihi[i][j] = 0.0; cffrelohi[i][j] = 0.0;
            cffrehilo[i][j] = 0.0; cffrelolo[i][j] = 0.0;
            cffimhihi[i][j] = 0.0; cffimlohi[i][j] = 0.0;
            cffimhilo[i][j] = 0.0; cffimlolo[i][j] = 0.0;
         }
      }
      CPU_cmplx4_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,
          accrehihi,accrelohi,accrehilo,accrelolo,
          accimhihi,accimlohi,accimhilo,accimlolo,
          inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
          inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
          outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
          outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)  // reset the coefficients
      {
         cffrehihi[i][0] = 1.0; cffrelohi[i][0] = 0.0;
         cffrehilo[i][0] = 0.0; cffrelolo[i][0] = 0.0;
         cffimhihi[i][0] = 0.0; cffimlohi[i][0] = 0.0;
         cffimhilo[i][0] = 0.0; cffimlolo[i][0] = 0.0;

         for(int j=1; j<degp1; j++)
         {
            cffrehihi[i][j] = 0.0; cffrelohi[i][j] = 0.0;
            cffrehilo[i][j] = 0.0; cffrelolo[i][j] = 0.0;
            cffimhihi[i][j] = 0.0; cffimlohi[i][j] = 0.0;
            cffimhilo[i][j] = 0.0; cffimlolo[i][j] = 0.0;
         }
      }
      GPU_cmplx4_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,
          accrehihi,accrelohi,accrehilo,accrelolo,
          accimhihi,accimlohi,accimhilo,accimlolo,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
          outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
          outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
          vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;
      for(int k=0; k<dim; k++) // monomial k
         for(int i=0; i<=dim; i++)
            for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "output_h[" << k << "][" << i << "][" << j << "] : "
                    << outputrehihi_h[k][i][j] << "  "
                    << outputrelohi_h[k][i][j] << endl << "  "
                    << outputrehilo_h[k][i][j] << "  "
                    << outputrelolo_h[k][i][j] << endl << "  "
                    << outputimhihi_h[k][i][j] << "  "
                    << outputimlohi_h[k][i][j] << endl << "  "
                    << outputimhilo_h[k][i][j] << "  "
                    << outputimlolo_h[k][i][j] << endl;
               cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                    << outputrehihi_d[k][i][j] << "  "
                    << outputrelohi_d[k][i][j] << endl << "  "
                    << outputrehilo_d[k][i][j] << "  "
                    << outputrelolo_d[k][i][j] << endl << "  "
                    << outputimhihi_d[k][i][j] << "  "
                    << outputimlohi_d[k][i][j] << endl << "  "
                    << outputimhilo_d[k][i][j] << "  "
                    << outputimlolo_d[k][i][j] << endl;
            }
            errsum += abs(outputrehihi_h[k][i][j] - outputrehihi_d[k][i][j])
                    + abs(outputrelohi_h[k][i][j] - outputrelohi_d[k][i][j])
                    + abs(outputrehilo_h[k][i][j] - outputrehilo_d[k][i][j])
                    + abs(outputrelolo_h[k][i][j] - outputrelolo_d[k][i][j])
                    + abs(outputimhihi_h[k][i][j] - outputimhihi_d[k][i][j])
                    + abs(outputimlohi_h[k][i][j] - outputimlohi_d[k][i][j])
                    + abs(outputimhilo_h[k][i][j] - outputimhilo_d[k][i][j])
                    + abs(outputimlolo_h[k][i][j] - outputimlolo_d[k][i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalrehihi_h[i][j][k] = 0.0; jacvalimhihi_h[i][j][k] = 0.0;
            jacvalrelohi_h[i][j][k] = 0.0; jacvalimlohi_h[i][j][k] = 0.0;
            jacvalrehilo_h[i][j][k] = 0.0; jacvalimhilo_h[i][j][k] = 0.0;
            jacvalrelolo_h[i][j][k] = 0.0; jacvalimlolo_h[i][j][k] = 0.0;
            jacvalrehihi_d[i][j][k] = 0.0; jacvalimhihi_d[i][j][k] = 0.0;
            jacvalrelohi_d[i][j][k] = 0.0; jacvalimlohi_d[i][j][k] = 0.0;
            jacvalrehilo_d[i][j][k] = 0.0; jacvalimhilo_d[i][j][k] = 0.0;
            jacvalrelolo_d[i][j][k] = 0.0; jacvalimlolo_d[i][j][k] = 0.0;
         }

   if((mode == 1) || (mode == 2))
      cmplx4_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
          outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
          funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
          funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
          rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
          rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
          jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
          jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
          vrblvl);

   if((mode == 0) || (mode == 2))
      cmplx4_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
          outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
          funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
          funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
          rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
          rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
          jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
          jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
          vrblvl);

   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "funval_h[" << i << "][" << j << "] : "
                    << funvalrehihi_h[i][j] << "  "
                    << funvalrelohi_h[i][j] << endl << "  "
                    << funvalrehilo_h[i][j] << "  "
                    << funvalrelolo_h[i][j] << endl << "  "
                    << funvalimhihi_h[i][j] << "  "
                    << funvalimlohi_h[i][j] << endl << "  "
                    << funvalimhilo_h[i][j] << "  "
                    << funvalimlolo_h[i][j] << endl;
               cout << "funval_d[" << i << "][" << j << "] : "
                    << funvalrehihi_d[i][j] << "  "
                    << funvalrelohi_d[i][j] << endl << "  "
                    << funvalrehilo_d[i][j] << "  "
                    << funvalrelolo_d[i][j] << endl << "  "
                    << funvalimhihi_d[i][j] << "  "
                    << funvalimlohi_d[i][j] << endl << "  "
                    << funvalimhilo_d[i][j] << "  "
                    << funvalimlolo_d[i][j] << endl;
            }
            errsum += abs(funvalrehihi_h[i][j] - funvalrehihi_d[i][j])
                    + abs(funvalrelohi_h[i][j] - funvalrelohi_d[i][j])
                    + abs(funvalrehilo_h[i][j] - funvalrehilo_d[i][j])
                    + abs(funvalrelolo_h[i][j] - funvalrelolo_d[i][j])
                    + abs(funvalimhihi_h[i][j] - funvalimhihi_d[i][j])
                    + abs(funvalimlohi_h[i][j] - funvalimlohi_d[i][j])
                    + abs(funvalimhilo_h[i][j] - funvalimhilo_d[i][j])
                    + abs(funvalimlolo_h[i][j] - funvalimlolo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            for(int k=0; k<dim; k++)
            {
               if(vrblvl > 1)
               {
                  cout << "jacval_h[" << i << "][" << j << "][" << k << "] : "
                       << jacvalrehihi_h[i][j][k] << "  "
                       << jacvalrelohi_h[i][j][k] << endl << "  "
                       << jacvalrehilo_h[i][j][k] << "  "
                       << jacvalrelolo_h[i][j][k] << endl << "  "
                       << jacvalimhihi_h[i][j][k] << "  "
                       << jacvalimlohi_h[i][j][k] << endl << "  "
                       << jacvalimhilo_h[i][j][k] << "  "
                       << jacvalimlolo_h[i][j][k] << endl;
                  cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                       << jacvalrehihi_d[i][j][k] << "  "
                       << jacvalrelohi_d[i][j][k] << endl << "  "
                       << jacvalrehilo_d[i][j][k] << "  "
                       << jacvalrelolo_d[i][j][k] << endl << "  "
                       << jacvalimhihi_d[i][j][k] << "  "
                       << jacvalimlohi_d[i][j][k] << endl << "  "
                       << jacvalimhilo_d[i][j][k] << "  "
                       << jacvalimlolo_d[i][j][k] << endl;
               }
             errsum += abs(jacvalrehihi_h[i][j][k] - jacvalrehihi_d[i][j][k])
                     + abs(jacvalrelohi_h[i][j][k] - jacvalrelohi_d[i][j][k])
                     + abs(jacvalrehilo_h[i][j][k] - jacvalrehilo_d[i][j][k])
                     + abs(jacvalrelolo_h[i][j][k] - jacvalrelolo_d[i][j][k])
                     + abs(jacvalimhihi_h[i][j][k] - jacvalimhihi_d[i][j][k])
                     + abs(jacvalimlohi_h[i][j][k] - jacvalimlohi_d[i][j][k])
                     + abs(jacvalimhilo_h[i][j][k] - jacvalimhilo_d[i][j][k])
                     + abs(jacvalimlolo_h[i][j][k] - jacvalimlolo_d[i][j][k]);
            }
         }
      }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = 0.0;
      for(int i=0; i<degp1; i++)
      {
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "rhs_h[" << i << "][" << j << "] : "
                    << rhsrehihi_h[i][j] << "  "
                    << rhsrelohi_h[i][j] << endl << "  "
                    << rhsrehilo_h[i][j] << "  "
                    << rhsrelolo_h[i][j] << endl << "  "
                    << rhsimhihi_h[i][j] << "  "
                    << rhsimlohi_h[i][j] << endl << "  "
                    << rhsimhilo_h[i][j] << "  "
                    << rhsimlolo_h[i][j] << endl;
               cout << "rhs_d[" << i << "][" << j << "] : "
                    << rhsrehihi_d[i][j] << "  "
                    << rhsrelohi_d[i][j] << endl << "  "
                    << rhsrehilo_d[i][j] << "  "
                    << rhsrelolo_d[i][j] << endl << "  "
                    << rhsimhihi_d[i][j] << "  "
                    << rhsimlohi_d[i][j] << endl << "  "
                    << rhsimhilo_d[i][j] << "  "
                    << rhsimlolo_d[i][j] << endl;
            }
            errsum += abs(rhsrehihi_h[i][j] - rhsrehihi_d[i][j])
                    + abs(rhsrelohi_h[i][j] - rhsrelohi_d[i][j])
                    + abs(rhsrehilo_h[i][j] - rhsrehilo_d[i][j])
                    + abs(rhsrelolo_h[i][j] - rhsrelolo_d[i][j])
                    + abs(rhsimhihi_h[i][j] - rhsimhihi_d[i][j])
                    + abs(rhsimlohi_h[i][j] - rhsimlohi_d[i][j])
                    + abs(rhsimhilo_h[i][j] - rhsimhilo_d[i][j])
                    + abs(rhsimlolo_h[i][j] - rhsimlolo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }
   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhsrehihi_h[i][j] = rhsrehihi_h[i][j];
         urhsrehilo_h[i][j] = rhsrehilo_h[i][j];
         urhsimhihi_h[i][j] = rhsimhihi_h[i][j];
         urhsimhilo_h[i][j] = rhsimhilo_h[i][j];
         urhsrelohi_h[i][j] = rhsrelohi_h[i][j];
         urhsrelolo_h[i][j] = rhsrelolo_h[i][j];
         urhsimlohi_h[i][j] = rhsimlohi_h[i][j];
         urhsimlolo_h[i][j] = rhsimlolo_h[i][j];
         urhsrehihi_d[i][j] = rhsrehihi_d[i][j];
         urhsrehilo_d[i][j] = rhsrehilo_d[i][j];
         urhsimhihi_d[i][j] = rhsimhihi_d[i][j];
         urhsimhilo_d[i][j] = rhsimhilo_d[i][j];
         urhsrelohi_d[i][j] = rhsrelohi_d[i][j];
         urhsrelolo_d[i][j] = rhsrelolo_d[i][j];
         urhsimlohi_d[i][j] = rhsimlohi_d[i][j];
         urhsimlolo_d[i][j] = rhsimlolo_d[i][j];
      }

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solrehihi_h[i][j] = 0.0; solimhihi_h[i][j] = 0.0;
         solrelohi_h[i][j] = 0.0; solimlohi_h[i][j] = 0.0;
         solrehilo_h[i][j] = 0.0; solimhilo_h[i][j] = 0.0;
         solrelolo_h[i][j] = 0.0; solimlolo_h[i][j] = 0.0;
         solrehihi_d[i][j] = 0.0; solimhihi_d[i][j] = 0.0;
         solrelohi_d[i][j] = 0.0; solimlohi_d[i][j] = 0.0;
         solrehilo_d[i][j] = 0.0; solimhilo_d[i][j] = 0.0;
         solrelolo_d[i][j] = 0.0; solimlolo_d[i][j] = 0.0;
      }

   if((mode == 1) || (mode == 2))
   {
      CPU_cmplx4_qrbs_solve
         (dim,degp1,
          jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
          jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
          urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
          urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
          solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
          solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
          workmatrehihi,workmatrelohi,workmatrehilo,workmatrelolo,
          workmatimhihi,workmatimlohi,workmatimhilo,workmatimlolo,
          Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
          Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
          Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
          Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
          workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
          workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,vrblvl);
 
      if(vrblvl > 0)
      {
         CPU_cmplx4_linear_residue
            (dim,degp1,
             jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
             jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
             rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
             rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
             solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
             solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
             resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
             resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);
         cout << "maximum residual : " << *resmaxhihi << endl;
      }
      cmplx4_update_series
         (dim,degp1,
          inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
          inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
          solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
          solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      GPU_cmplx4_bals_solve
         (dim,degp1,szt,nbt,
          jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
          jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
          Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
          Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
          Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
          Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
          urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
          urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
          solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
          solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,vrblvl);

      if(vrblvl > 0)
      {
         CPU_cmplx4_linear_residue
            (dim,degp1,
             jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
             jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
             rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
             rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
             solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
             solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
             resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
             resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);
         cout << "maximum residual : " << *resmaxhihi << endl;
      }
      cmplx4_update_series
         (dim,degp1,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
          solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
          solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "Q_h[" << i << "][" << j << "] : "
                    << Qrehihi_h[i][j] << "  "
                    << Qrelohi_h[i][j] << endl << "  "
                    << Qrehilo_h[i][j] << "  "
                    << Qrelolo_h[i][j] << endl << "  "
                    << Qimhihi_h[i][j] << "  "
                    << Qimlohi_h[i][j] << endl << "  "
                    << Qimhilo_h[i][j] << "  "
                    << Qimlolo_h[i][j] << endl;
               cout << "Q_d[" << i << "][" << j << "] : "
                    << Qrehihi_d[i][j] << "  "
                    << Qrelohi_d[i][j] << endl << "  "
                    << Qrehilo_d[i][j] << "  "
                    << Qrelolo_d[i][j] << endl << "  "
                    << Qimhihi_d[i][j] << "  "
                    << Qimlohi_d[i][j] << endl << "  "
                    << Qimhilo_d[i][j] << "  "
                    << Qimlolo_d[i][j] << endl;
            }
            errsum += abs(Qrehihi_h[i][j] - Qrehihi_d[i][j])
                    + abs(Qrelohi_h[i][j] - Qrelohi_d[i][j])
                    + abs(Qrehilo_h[i][j] - Qrehilo_d[i][j])
                    + abs(Qrelolo_h[i][j] - Qrelolo_d[i][j])
                    + abs(Qimhihi_h[i][j] - Qimhihi_d[i][j])
                    + abs(Qimlohi_h[i][j] - Qimlohi_d[i][j])
                    + abs(Qimhilo_h[i][j] - Qimhilo_d[i][j])
                    + abs(Qimlolo_h[i][j] - Qimlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "R_h[" << i << "][" << j << "] : "
                    << Rrehihi_h[i][j] << "  "
                    << Rrelohi_h[i][j] << endl << "  "
                    << Rrehilo_h[i][j] << "  "
                    << Rrelolo_h[i][j] << endl << "  "
                    << Rimhihi_h[i][j] << "  "
                    << Rimlohi_h[i][j] << endl << "  "
                    << Rimhilo_h[i][j] << "  "
                    << Rimlolo_h[i][j] << endl;
               cout << "R_d[" << i << "][" << j << "] : "
                    << Rrehihi_d[i][j] << "  "
                    << Rrelohi_d[i][j] << endl << "  "
                    << Rrehilo_d[i][j] << "  "
                    << Rrelolo_d[i][j] << endl << "  "
                    << Rimhihi_d[i][j] << "  "
                    << Rimlohi_d[i][j] << endl << "  "
                    << Rimhilo_d[i][j] << "  "
                    << Rimlolo_d[i][j] << endl;
            }
            errsum += abs(Rrehihi_h[i][j] - Rrehihi_d[i][j])
                    + abs(Rrelohi_h[i][j] - Rrelohi_d[i][j])
                    + abs(Rrehilo_h[i][j] - Rrehilo_d[i][j])
                    + abs(Rrelolo_h[i][j] - Rrelolo_d[i][j])
                    + abs(Rimhihi_h[i][j] - Rimhihi_d[i][j])
                    + abs(Rimlohi_h[i][j] - Rimlohi_d[i][j])
                    + abs(Rimhilo_h[i][j] - Rimhilo_d[i][j])
                    + abs(Rimlolo_h[i][j] - Rimlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "urhs_h[" << i << "][" << j << "] : "
                    << urhsrehihi_h[i][j] << "  "
                    << urhsrelohi_h[i][j] << endl << "  "
                    << urhsrehilo_h[i][j] << "  "
                    << urhsrelolo_h[i][j] << endl << "  "
                    << urhsimhihi_h[i][j] << "  "
                    << urhsimlohi_h[i][j] << endl << "  "
                    << urhsimhilo_h[i][j] << "  "
                    << urhsimlolo_h[i][j] << endl;
               cout << "urhs_d[" << i << "][" << j << "] : "
                    << urhsrehihi_d[i][j] << "  "
                    << urhsrelohi_d[i][j] << endl << "  "
                    << urhsrehilo_d[i][j] << "  "
                    << urhsrelolo_d[i][j] << endl << "  "
                    << urhsimhihi_d[i][j] << "  "
                    << urhsimlohi_d[i][j] << endl << "  "
                    << urhsimhilo_d[i][j] << "  "
                    << urhsimlolo_d[i][j] << endl;
            }
            errsum += abs(urhsrehihi_h[i][j] - urhsrehihi_d[i][j])
                    + abs(urhsrelohi_h[i][j] - urhsrelohi_d[i][j])
                    + abs(urhsrehilo_h[i][j] - urhsrehilo_d[i][j])
                    + abs(urhsrelolo_h[i][j] - urhsrelolo_d[i][j])
                    + abs(urhsimhihi_h[i][j] - urhsimhihi_d[i][j])
                    + abs(urhsimlohi_h[i][j] - urhsimlohi_d[i][j])
                    + abs(urhsimhilo_h[i][j] - urhsimhilo_d[i][j])
                    + abs(urhsimlolo_h[i][j] - urhsimlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            if(vrblvl > 1)
            {
               cout << "sol_h[" << i << "][" << j << "] : "
                    << solrehihi_h[i][j] << "  "
                    << solrelohi_h[i][j] << endl << "  "
                    << solrehilo_h[i][j] << "  "
                    << solrelolo_h[i][j] << endl << "  "
                    << solimhihi_h[i][j] << "  "
                    << solimlohi_h[i][j] << endl << "  "
                    << solimhilo_h[i][j] << "  "
                    << solimlolo_h[i][j] << endl;
               cout << "sol_d[" << i << "][" << j << "] : "
                    << solrehihi_d[i][j] << "  "
                    << solrelohi_d[i][j] << endl << "  "
                    << solrehilo_d[i][j] << "  "
                    << solrelolo_d[i][j] << endl << "  "
                    << solimhihi_d[i][j] << "  "
                    << solimlohi_d[i][j] << endl << "  "
                    << solimhilo_d[i][j] << "  "
                    << solimlolo_d[i][j] << endl;
            }
            errsum += abs(solrehihi_h[i][j] - solrehihi_d[i][j])
                    + abs(solrelohi_h[i][j] - solrelohi_d[i][j])
                    + abs(solrehilo_h[i][j] - solrehilo_d[i][j])
                    + abs(solrelolo_h[i][j] - solrelolo_d[i][j])
                    + abs(solimhihi_h[i][j] - solimhihi_d[i][j])
                    + abs(solimlohi_h[i][j] - solimlohi_d[i][j])
                    + abs(solimhilo_h[i][j] - solimhilo_d[i][j])
                    + abs(solimlolo_h[i][j] - solimlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
            if(vrblvl > 1)
            {
               cout << "input_h[" << i << "][" << j << "] : "
                    << inputrehihi_h[i][j] << "  "
                    << inputrelohi_h[i][j] << endl << "  "
                    << inputrehilo_h[i][j] << "  "
                    << inputrelolo_h[i][j] << endl << "  "
                    << inputimhihi_h[i][j] << "  "
                    << inputimlohi_h[i][j] << endl << "  "
                    << inputimhilo_h[i][j] << "  "
                    << inputimlolo_h[i][j] << endl;
               cout << "input_d[" << i << "][" << j << "] : "
                    << inputrehihi_d[i][j] << "  "
                    << inputrelohi_d[i][j] << endl << "  "
                    << inputrehilo_d[i][j] << "  "
                    << inputrelolo_d[i][j] << endl << "  "
                    << inputimhihi_d[i][j] << "  "
                    << inputimlohi_d[i][j] << endl << "  "
                    << inputimhilo_d[i][j] << "  "
                    << inputimlolo_d[i][j] << endl;
            }
            errsum += abs(inputrehihi_h[i][j] - inputrehihi_d[i][j])
                    + abs(inputrelohi_h[i][j] - inputrelohi_d[i][j])
                    + abs(inputrehilo_h[i][j] - inputrehilo_d[i][j])
                    + abs(inputrelolo_h[i][j] - inputrelolo_d[i][j])
                    + abs(inputimhihi_h[i][j] - inputimhihi_d[i][j])
                    + abs(inputimlohi_h[i][j] - inputimlohi_d[i][j])
                    + abs(inputimhilo_h[i][j] - inputimhilo_d[i][j])
                    + abs(inputimlolo_h[i][j] - inputimlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl4_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhihi_h = new double*[dim];
   double **inputlohi_h = new double*[dim];
   double **inputhilo_h = new double*[dim];
   double **inputlolo_h = new double*[dim];
   double **inputhihi_d = new double*[dim];
   double **inputlohi_d = new double*[dim];
   double **inputhilo_d = new double*[dim];
   double **inputlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihi_h[i] = new double[degp1];
      inputlohi_h[i] = new double[degp1];
      inputhilo_h[i] = new double[degp1];
      inputlolo_h[i] = new double[degp1];
      inputhihi_d[i] = new double[degp1];
      inputlohi_d[i] = new double[degp1];
      inputhilo_d[i] = new double[degp1];
      inputlolo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acchihi = new double[degp1]; // accumulated power series
   double *acclohi = new double[degp1];
   double *acchilo = new double[degp1];
   double *acclolo = new double[degp1];
   double **cffhihi = new double*[dim]; // the coefficients of monomials
   double **cfflohi = new double*[dim];
   double **cffhilo = new double*[dim];
   double **cfflolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihi[i] = new double[degp1];
      cfflohi[i] = new double[degp1];
      cffhilo[i] = new double[degp1];
      cfflolo[i] = new double[degp1];
   }
   double ***outputhihi_h = new double**[dim];
   double ***outputlohi_h = new double**[dim];
   double ***outputhilo_h = new double**[dim];
   double ***outputlolo_h = new double**[dim];
   double ***outputhihi_d = new double**[dim];
   double ***outputlohi_d = new double**[dim];
   double ***outputhilo_d = new double**[dim];
   double ***outputlolo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihi_h[i] = new double*[dim+1];
      outputlohi_h[i] = new double*[dim+1];
      outputhilo_h[i] = new double*[dim+1];
      outputlolo_h[i] = new double*[dim+1];
      outputhihi_d[i] = new double*[dim+1];
      outputlohi_d[i] = new double*[dim+1];
      outputhilo_d[i] = new double*[dim+1];
      outputlolo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihi_h[i][j] = new double[degp1];
         outputlohi_h[i][j] = new double[degp1];
         outputhilo_h[i][j] = new double[degp1];
         outputlolo_h[i][j] = new double[degp1];
         outputhihi_d[i][j] = new double[degp1];
         outputlohi_d[i][j] = new double[degp1];
         outputhilo_d[i][j] = new double[degp1];
         outputlolo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihi_h = new double*[dim];
   double **funvallohi_h = new double*[dim];
   double **funvalhilo_h = new double*[dim];
   double **funvallolo_h = new double*[dim];
   double **funvalhihi_d = new double*[dim];
   double **funvallohi_d = new double*[dim];
   double **funvalhilo_d = new double*[dim];
   double **funvallolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhihi_h[i] = new double[degp1];
      funvallohi_h[i] = new double[degp1];
      funvalhilo_h[i] = new double[degp1];
      funvallolo_h[i] = new double[degp1];
      funvalhihi_d[i] = new double[degp1];
      funvallohi_d[i] = new double[degp1];
      funvalhilo_d[i] = new double[degp1];
      funvallolo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihi_h = new double**[degp1];
   double ***jacvallohi_h = new double**[degp1];
   double ***jacvalhilo_h = new double**[degp1];
   double ***jacvallolo_h = new double**[degp1];
   double ***jacvalhihi_d = new double**[degp1];
   double ***jacvallohi_d = new double**[degp1];
   double ***jacvalhilo_d = new double**[degp1];
   double ***jacvallolo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhihi_h[i] = new double*[dim];
      jacvallohi_h[i] = new double*[dim];
      jacvalhilo_h[i] = new double*[dim];
      jacvallolo_h[i] = new double*[dim];
      jacvalhihi_d[i] = new double*[dim];
      jacvallohi_d[i] = new double*[dim];
      jacvalhilo_d[i] = new double*[dim];
      jacvallolo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhihi_h[i][j] = new double[dim];
         jacvallohi_h[i][j] = new double[dim];
         jacvalhilo_h[i][j] = new double[dim];
         jacvallolo_h[i][j] = new double[dim];
         jacvalhihi_d[i][j] = new double[dim];
         jacvallohi_d[i][j] = new double[dim];
         jacvalhilo_d[i][j] = new double[dim];
         jacvallolo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihi_h = new double*[degp1];
   double **sollohi_h = new double*[degp1];
   double **solhilo_h = new double*[degp1];
   double **sollolo_h = new double*[degp1];
   double **solhihi_d = new double*[degp1];
   double **sollohi_d = new double*[degp1];
   double **solhilo_d = new double*[degp1];
   double **sollolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihi_h[i] = new double[dim];
      sollohi_h[i] = new double[dim];
      solhilo_h[i] = new double[dim];
      sollolo_h[i] = new double[dim];
      solhihi_d[i] = new double[dim];
      sollohi_d[i] = new double[dim];
      solhilo_d[i] = new double[dim];
      sollolo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihi_h = new double*[degp1];
   double **rhslohi_h = new double*[degp1];
   double **rhshilo_h = new double*[degp1];
   double **rhslolo_h = new double*[degp1];
   double **rhshihi_d = new double*[degp1];
   double **rhslohi_d = new double*[degp1];
   double **rhshilo_d = new double*[degp1];
   double **rhslolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihi_h[i] = new double[dim];
      rhslohi_h[i] = new double[dim];
      rhshilo_h[i] = new double[dim];
      rhslolo_h[i] = new double[dim];
      rhshihi_d[i] = new double[dim];
      rhslohi_d[i] = new double[dim];
      rhshilo_d[i] = new double[dim];
      rhslolo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihi_h = new double*[degp1];
   double **urhslohi_h = new double*[degp1];
   double **urhshilo_h = new double*[degp1];
   double **urhslolo_h = new double*[degp1];
   double **urhshihi_d = new double*[degp1];
   double **urhslohi_d = new double*[degp1];
   double **urhshilo_d = new double*[degp1];
   double **urhslolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhshihi_h[i] = new double[dim];
      urhslohi_h[i] = new double[dim];
      urhshilo_h[i] = new double[dim];
      urhslolo_h[i] = new double[dim];
      urhshihi_d[i] = new double[dim];
      urhslohi_d[i] = new double[dim];
      urhshilo_d[i] = new double[dim];
      urhslolo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihi = new double*[dim];
   double **workmatlohi = new double*[dim];
   double **workmathilo = new double*[dim];
   double **workmatlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihi[i] = new double[dim];
      workmatlohi[i] = new double[dim];
      workmathilo[i] = new double[dim];
      workmatlolo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechihi = new double[dim];
   double *workveclohi = new double[dim];
   double *workvechilo = new double[dim];
   double *workveclolo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihi = new double*[degp1];
   double **workrhslohi = new double*[degp1];
   double **workrhshilo = new double*[degp1];
   double **workrhslolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihi[i] = new double[dim];
      workrhslohi[i] = new double[dim];
      workrhshilo[i] = new double[dim];
      workrhslolo[i] = new double[dim];
   }
   double **resvechihi = new double*[degp1];
   double **resveclohi = new double*[degp1];
   double **resvechilo = new double*[degp1];
   double **resveclolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihi[i] = new double[dim];
      resveclohi[i] = new double[dim];
      resvechilo[i] = new double[dim];
      resveclolo[i] = new double[dim];
   }
   double resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo;
   double **Qhihi_h = new double*[dim];
   double **Qlohi_h = new double*[dim];
   double **Qhilo_h = new double*[dim];
   double **Qlolo_h = new double*[dim];
   double **Qhihi_d = new double*[dim];
   double **Qlohi_d = new double*[dim];
   double **Qhilo_d = new double*[dim];
   double **Qlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qhihi_h[i] = new double[dim];
      Qlohi_h[i] = new double[dim];
      Qhilo_h[i] = new double[dim];
      Qlolo_h[i] = new double[dim];
      Qhihi_d[i] = new double[dim];
      Qlohi_d[i] = new double[dim];
      Qhilo_d[i] = new double[dim];
      Qlolo_d[i] = new double[dim];
   }
   double **Rhihi_h = new double*[dim];
   double **Rlohi_h = new double*[dim];
   double **Rhilo_h = new double*[dim];
   double **Rlolo_h = new double*[dim];
   double **Rhihi_d = new double*[dim];
   double **Rlohi_d = new double*[dim];
   double **Rhilo_d = new double*[dim];
   double **Rlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhihi_h[i] = new double[dim];
      Rlohi_h[i] = new double[dim];
      Rhilo_h[i] = new double[dim];
      Rlolo_h[i] = new double[dim];
      Rhihi_d[i] = new double[dim];
      Rlohi_d[i] = new double[dim];
      Rhilo_d[i] = new double[dim];
      Rlolo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl4_unit_series_vector
      (dim,deg,inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h);
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihi_d[i][j] = inputhihi_h[i][j];
         inputlohi_d[i][j] = inputlohi_h[i][j];
         inputhilo_d[i][j] = inputhilo_h[i][j];
         inputlolo_d[i][j] = inputlolo_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihi_h[i][0] << "  "
                            << inputlohi_h[i][0] << endl;
         cout << "     " << inputhilo_h[i][0] << "  "
                         << inputlolo_h[i][0] << endl;
      }
   }
   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;
/*
      dbl4_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
           inputhihi, inputlohi, inputhilo, inputlolo,
          outputhihi,outputlohi,outputhilo,outputlolo,
          funvalhihi,funvallohi,funvalhilo,funvallolo,
          jacvalhihi,jacvallohi,jacvalhilo,jacvallolo,
          rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          workvechihi,workveclohi,workvechilo,workveclolo,
          workrhshihi,workrhslohi,workrhshilo,workrhslolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,ipvt,vrblvl);
 */
      dbl4_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
          inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
          outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
          funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
          funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
          jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
          jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
          rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
          rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
          urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
          urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
          solhihi_h,sollohi_h,solhilo_h,sollolo_h,
          solhihi_d,sollohi_d,solhilo_d,sollolo_d,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
          Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          workvechihi,workveclohi,workvechilo,workveclolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,vrblvl,mode);
   }
   return 0;
}

int test_dbl4_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputrehihi_h = new double*[dim];
   double **inputrelohi_h = new double*[dim];
   double **inputrehilo_h = new double*[dim];
   double **inputrelolo_h = new double*[dim];
   double **inputimhihi_h = new double*[dim];
   double **inputimlohi_h = new double*[dim];
   double **inputimhilo_h = new double*[dim];
   double **inputimlolo_h = new double*[dim];
   double **inputrehihi_d = new double*[dim];
   double **inputrelohi_d = new double*[dim];
   double **inputrehilo_d = new double*[dim];
   double **inputrelolo_d = new double*[dim];
   double **inputimhihi_d = new double*[dim];
   double **inputimlohi_d = new double*[dim];
   double **inputimhilo_d = new double*[dim];
   double **inputimlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputrehihi_h[i] = new double[degp1];
       inputrelohi_h[i] = new double[degp1];
       inputrehilo_h[i] = new double[degp1];
       inputrelolo_h[i] = new double[degp1];
       inputimhihi_h[i] = new double[degp1];
       inputimlohi_h[i] = new double[degp1];
       inputimhilo_h[i] = new double[degp1];
       inputimlolo_h[i] = new double[degp1];
       inputrehihi_d[i] = new double[degp1];
       inputrelohi_d[i] = new double[degp1];
       inputrehilo_d[i] = new double[degp1];
       inputrelolo_d[i] = new double[degp1];
       inputimhihi_d[i] = new double[degp1];
       inputimlohi_d[i] = new double[degp1];
       inputimhilo_d[i] = new double[degp1];
       inputimlolo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *accrehihi = new double[degp1]; // accumulated power series
   double *accrelohi = new double[degp1];
   double *accrehilo = new double[degp1];
   double *accrelolo = new double[degp1];
   double *accimhihi = new double[degp1];
   double *accimlohi = new double[degp1];
   double *accimhilo = new double[degp1];
   double *accimlolo = new double[degp1];
   double **cffrehihi = new double*[dim]; // the coefficients of monomials
   double **cffrelohi = new double*[dim];
   double **cffrehilo = new double*[dim]; 
   double **cffrelolo = new double*[dim];
   double **cffimhihi = new double*[dim]; 
   double **cffimlohi = new double*[dim]; 
   double **cffimhilo = new double*[dim]; 
   double **cffimlolo = new double*[dim]; 

   for(int i=0; i<dim; i++)
   {
      cffrehihi[i] = new double[degp1];
      cffrelohi[i] = new double[degp1];
      cffrehilo[i] = new double[degp1];
      cffrelolo[i] = new double[degp1];
      cffimhihi[i] = new double[degp1];
      cffimlohi[i] = new double[degp1];
      cffimhilo[i] = new double[degp1];
      cffimlolo[i] = new double[degp1];
   }
   double ***outputrehihi_h = new double**[dim];
   double ***outputrelohi_h = new double**[dim];
   double ***outputrehilo_h = new double**[dim];
   double ***outputrelolo_h = new double**[dim];
   double ***outputimhihi_h = new double**[dim];
   double ***outputimlohi_h = new double**[dim];
   double ***outputimhilo_h = new double**[dim];
   double ***outputimlolo_h = new double**[dim];
   double ***outputrehihi_d = new double**[dim];
   double ***outputrelohi_d = new double**[dim];
   double ***outputrehilo_d = new double**[dim];
   double ***outputrelolo_d = new double**[dim];
   double ***outputimhihi_d = new double**[dim];
   double ***outputimlohi_d = new double**[dim];
   double ***outputimhilo_d = new double**[dim];
   double ***outputimlolo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehihi_h[i] = new double*[dim+1];
      outputrelohi_h[i] = new double*[dim+1];
      outputrehilo_h[i] = new double*[dim+1];
      outputrelolo_h[i] = new double*[dim+1];
      outputimhihi_h[i] = new double*[dim+1];
      outputimlohi_h[i] = new double*[dim+1];
      outputimhilo_h[i] = new double*[dim+1];
      outputimlolo_h[i] = new double*[dim+1];
      outputrehihi_d[i] = new double*[dim+1];
      outputrelohi_d[i] = new double*[dim+1];
      outputrehilo_d[i] = new double*[dim+1];
      outputrelolo_d[i] = new double*[dim+1];
      outputimhihi_d[i] = new double*[dim+1];
      outputimlohi_d[i] = new double*[dim+1];
      outputimhilo_d[i] = new double*[dim+1];
      outputimlolo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehihi_h[i][j] = new double[degp1];
         outputrelohi_h[i][j] = new double[degp1];
         outputrehilo_h[i][j] = new double[degp1];
         outputrelolo_h[i][j] = new double[degp1];
         outputimhihi_h[i][j] = new double[degp1];
         outputimlohi_h[i][j] = new double[degp1];
         outputimhilo_h[i][j] = new double[degp1];
         outputimlolo_h[i][j] = new double[degp1];
         outputrehihi_d[i][j] = new double[degp1];
         outputrelohi_d[i][j] = new double[degp1];
         outputrehilo_d[i][j] = new double[degp1];
         outputrelolo_d[i][j] = new double[degp1];
         outputimhihi_d[i][j] = new double[degp1];
         outputimlohi_d[i][j] = new double[degp1];
         outputimhilo_d[i][j] = new double[degp1];
         outputimlolo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalrehihi_h = new double*[dim];
   double **funvalrelohi_h = new double*[dim];
   double **funvalrehilo_h = new double*[dim];
   double **funvalrelolo_h = new double*[dim];
   double **funvalimhihi_h = new double*[dim];
   double **funvalimlohi_h = new double*[dim];
   double **funvalimhilo_h = new double*[dim];
   double **funvalimlolo_h = new double*[dim];
   double **funvalrehihi_d = new double*[dim];
   double **funvalrelohi_d = new double*[dim];
   double **funvalrehilo_d = new double*[dim];
   double **funvalrelolo_d = new double*[dim];
   double **funvalimhihi_d = new double*[dim];
   double **funvalimlohi_d = new double*[dim];
   double **funvalimhilo_d = new double*[dim];
   double **funvalimlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalrehihi_h[i] = new double[degp1];
      funvalrelohi_h[i] = new double[degp1];
      funvalrehilo_h[i] = new double[degp1];
      funvalrelolo_h[i] = new double[degp1];
      funvalimhihi_h[i] = new double[degp1];
      funvalimlohi_h[i] = new double[degp1];
      funvalimhilo_h[i] = new double[degp1];
      funvalimlolo_h[i] = new double[degp1];
      funvalrehihi_d[i] = new double[degp1];
      funvalrelohi_d[i] = new double[degp1];
      funvalrehilo_d[i] = new double[degp1];
      funvalrelolo_d[i] = new double[degp1];
      funvalimhihi_d[i] = new double[degp1];
      funvalimlohi_d[i] = new double[degp1];
      funvalimhilo_d[i] = new double[degp1];
      funvalimlolo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalrehihi_h = new double**[degp1];
   double ***jacvalrelohi_h = new double**[degp1];
   double ***jacvalrehilo_h = new double**[degp1];
   double ***jacvalrelolo_h = new double**[degp1];
   double ***jacvalimhihi_h = new double**[degp1];
   double ***jacvalimlohi_h = new double**[degp1];
   double ***jacvalimhilo_h = new double**[degp1];
   double ***jacvalimlolo_h = new double**[degp1];
   double ***jacvalrehihi_d = new double**[degp1];
   double ***jacvalrelohi_d = new double**[degp1];
   double ***jacvalrehilo_d = new double**[degp1];
   double ***jacvalrelolo_d = new double**[degp1];
   double ***jacvalimhihi_d = new double**[degp1];
   double ***jacvalimlohi_d = new double**[degp1];
   double ***jacvalimhilo_d = new double**[degp1];
   double ***jacvalimlolo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalrehihi_h[i] = new double*[dim];
      jacvalrelohi_h[i] = new double*[dim];
      jacvalrehilo_h[i] = new double*[dim];
      jacvalrelolo_h[i] = new double*[dim];
      jacvalimhihi_h[i] = new double*[dim];
      jacvalimlohi_h[i] = new double*[dim];
      jacvalimhilo_h[i] = new double*[dim];
      jacvalimlolo_h[i] = new double*[dim];
      jacvalrehihi_d[i] = new double*[dim];
      jacvalrelohi_d[i] = new double*[dim];
      jacvalrehilo_d[i] = new double*[dim];
      jacvalrelolo_d[i] = new double*[dim];
      jacvalimhihi_d[i] = new double*[dim];
      jacvalimlohi_d[i] = new double*[dim];
      jacvalimhilo_d[i] = new double*[dim];
      jacvalimlolo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalrehihi_h[i][j] = new double[dim];
         jacvalrelohi_h[i][j] = new double[dim];
         jacvalrehilo_h[i][j] = new double[dim];
         jacvalrelolo_h[i][j] = new double[dim];
         jacvalimhihi_h[i][j] = new double[dim];
         jacvalimlohi_h[i][j] = new double[dim];
         jacvalimhilo_h[i][j] = new double[dim];
         jacvalimlolo_h[i][j] = new double[dim];
         jacvalrehihi_d[i][j] = new double[dim];
         jacvalrelohi_d[i][j] = new double[dim];
         jacvalrehilo_d[i][j] = new double[dim];
         jacvalrelolo_d[i][j] = new double[dim];
         jacvalimhihi_d[i][j] = new double[dim];
         jacvalimlohi_d[i][j] = new double[dim];
         jacvalimhilo_d[i][j] = new double[dim];
         jacvalimlolo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solrehihi_h = new double*[degp1];
   double **solrelohi_h = new double*[degp1];
   double **solrehilo_h = new double*[degp1];
   double **solrelolo_h = new double*[degp1];
   double **solimhihi_h = new double*[degp1];
   double **solimlohi_h = new double*[degp1];
   double **solimhilo_h = new double*[degp1];
   double **solimlolo_h = new double*[degp1];
   double **solrehihi_d = new double*[degp1];
   double **solrelohi_d = new double*[degp1];
   double **solrehilo_d = new double*[degp1];
   double **solrelolo_d = new double*[degp1];
   double **solimhihi_d = new double*[degp1];
   double **solimlohi_d = new double*[degp1];
   double **solimhilo_d = new double*[degp1];
   double **solimlolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++) 
   {
      solrehihi_h[i] = new double[dim];
      solrelohi_h[i] = new double[dim];
      solrehilo_h[i] = new double[dim];
      solrelolo_h[i] = new double[dim];
      solimhihi_h[i] = new double[dim];
      solimlohi_h[i] = new double[dim];
      solimhilo_h[i] = new double[dim];
      solimlolo_h[i] = new double[dim];
      solrehihi_d[i] = new double[dim];
      solrelohi_d[i] = new double[dim];
      solrehilo_d[i] = new double[dim];
      solrelolo_d[i] = new double[dim];
      solimhihi_d[i] = new double[dim];
      solimlohi_d[i] = new double[dim];
      solimhilo_d[i] = new double[dim];
      solimlolo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsrehihi_h = new double*[degp1];
   double **rhsrelohi_h = new double*[degp1];
   double **rhsrehilo_h = new double*[degp1];
   double **rhsrelolo_h = new double*[degp1];
   double **rhsimhihi_h = new double*[degp1];
   double **rhsimlohi_h = new double*[degp1];
   double **rhsimhilo_h = new double*[degp1];
   double **rhsimlolo_h = new double*[degp1];
   double **rhsrehihi_d = new double*[degp1];
   double **rhsrelohi_d = new double*[degp1];
   double **rhsrehilo_d = new double*[degp1];
   double **rhsrelolo_d = new double*[degp1];
   double **rhsimhihi_d = new double*[degp1];
   double **rhsimlohi_d = new double*[degp1];
   double **rhsimhilo_d = new double*[degp1];
   double **rhsimlolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhsrehihi_h[i] = new double[dim];
      rhsrelohi_h[i] = new double[dim];
      rhsrehilo_h[i] = new double[dim];
      rhsrelolo_h[i] = new double[dim];
      rhsimhihi_h[i] = new double[dim];
      rhsimlohi_h[i] = new double[dim];
      rhsimhilo_h[i] = new double[dim];
      rhsimlolo_h[i] = new double[dim];
      rhsrehihi_d[i] = new double[dim];
      rhsrelohi_d[i] = new double[dim];
      rhsrehilo_d[i] = new double[dim];
      rhsrelolo_d[i] = new double[dim];
      rhsimhihi_d[i] = new double[dim];
      rhsimlohi_d[i] = new double[dim];
      rhsimhilo_d[i] = new double[dim];
      rhsimlolo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmatrehihi = new double*[dim];
   double **workmatrelohi = new double*[dim];
   double **workmatrehilo = new double*[dim];
   double **workmatrelolo = new double*[dim];
   double **workmatimhihi = new double*[dim];
   double **workmatimlohi = new double*[dim];
   double **workmatimhilo = new double*[dim];
   double **workmatimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmatrehihi[i] = new double[dim];
      workmatrelohi[i] = new double[dim];
      workmatrehilo[i] = new double[dim];
      workmatrelolo[i] = new double[dim];
      workmatimhihi[i] = new double[dim];
      workmatimlohi[i] = new double[dim];
      workmatimhilo[i] = new double[dim];
      workmatimlolo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvecrehihi = new double[dim];
   double *workvecrelohi = new double[dim];
   double *workvecrehilo = new double[dim];
   double *workvecrelolo = new double[dim];
   double *workvecimhihi = new double[dim];
   double *workvecimlohi = new double[dim];
   double *workvecimhilo = new double[dim];
   double *workvecimlolo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **urhsrehihi_h = new double*[degp1];
   double **urhsrelohi_h = new double*[degp1];
   double **urhsrehilo_h = new double*[degp1];
   double **urhsrelolo_h = new double*[degp1];
   double **urhsimhihi_h = new double*[degp1];
   double **urhsimlohi_h = new double*[degp1];
   double **urhsimhilo_h = new double*[degp1];
   double **urhsimlolo_h = new double*[degp1];
   double **urhsrehihi_d = new double*[degp1];
   double **urhsrelohi_d = new double*[degp1];
   double **urhsrehilo_d = new double*[degp1];
   double **urhsrelolo_d = new double*[degp1];
   double **urhsimhihi_d = new double*[degp1];
   double **urhsimlohi_d = new double*[degp1];
   double **urhsimhilo_d = new double*[degp1];
   double **urhsimlolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhsrehihi_h[i] = new double[dim];
      urhsrelohi_h[i] = new double[dim];
      urhsrehilo_h[i] = new double[dim];
      urhsrelolo_h[i] = new double[dim];
      urhsimhihi_h[i] = new double[dim];
      urhsimlohi_h[i] = new double[dim];
      urhsimhilo_h[i] = new double[dim];
      urhsimlolo_h[i] = new double[dim];
      urhsrehihi_d[i] = new double[dim];
      urhsrelohi_d[i] = new double[dim];
      urhsrehilo_d[i] = new double[dim];
      urhsrelolo_d[i] = new double[dim];
      urhsimhihi_d[i] = new double[dim];
      urhsimlohi_d[i] = new double[dim];
      urhsimhilo_d[i] = new double[dim];
      urhsimlolo_d[i] = new double[dim];
   }
   double **resvecrehihi = new double*[degp1];
   double **resvecrelohi = new double*[degp1];
   double **resvecrehilo = new double*[degp1];
   double **resvecrelolo = new double*[degp1];
   double **resvecimhihi = new double*[degp1];
   double **resvecimlohi = new double*[degp1];
   double **resvecimhilo = new double*[degp1];
   double **resvecimlolo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecrehihi[i] = new double[dim];
      resvecrelohi[i] = new double[dim];
      resvecrehilo[i] = new double[dim];
      resvecrelolo[i] = new double[dim];
      resvecimhihi[i] = new double[dim];
      resvecimlohi[i] = new double[dim];
      resvecimhilo[i] = new double[dim];
      resvecimlolo[i] = new double[dim];
   }
   double resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo;
   double **Qrehihi_h = new double*[dim];
   double **Qrelohi_h = new double*[dim];
   double **Qrehilo_h = new double*[dim];
   double **Qrelolo_h = new double*[dim];
   double **Qimhihi_h = new double*[dim];
   double **Qimlohi_h = new double*[dim];
   double **Qimhilo_h = new double*[dim];
   double **Qimlolo_h = new double*[dim];
   double **Qrehihi_d = new double*[dim];
   double **Qrelohi_d = new double*[dim];
   double **Qrehilo_d = new double*[dim];
   double **Qrelolo_d = new double*[dim];
   double **Qimhihi_d = new double*[dim];
   double **Qimlohi_d = new double*[dim];
   double **Qimhilo_d = new double*[dim];
   double **Qimlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qrehihi_h[i] = new double[dim];
      Qrelohi_h[i] = new double[dim];
      Qrehilo_h[i] = new double[dim];
      Qrelolo_h[i] = new double[dim];
      Qimhihi_h[i] = new double[dim];
      Qimlohi_h[i] = new double[dim];
      Qimhilo_h[i] = new double[dim];
      Qimlolo_h[i] = new double[dim];
      Qrehihi_d[i] = new double[dim];
      Qrelohi_d[i] = new double[dim];
      Qrehilo_d[i] = new double[dim];
      Qrelolo_d[i] = new double[dim];
      Qimhihi_d[i] = new double[dim];
      Qimlohi_d[i] = new double[dim];
      Qimhilo_d[i] = new double[dim];
      Qimlolo_d[i] = new double[dim];
   }
   double **Rrehihi_h = new double*[dim];
   double **Rrelohi_h = new double*[dim];
   double **Rrehilo_h = new double*[dim];
   double **Rrelolo_h = new double*[dim];
   double **Rimhihi_h = new double*[dim];
   double **Rimlohi_h = new double*[dim];
   double **Rimhilo_h = new double*[dim];
   double **Rimlolo_h = new double*[dim];
   double **Rrehihi_d = new double*[dim];
   double **Rrelohi_d = new double*[dim];
   double **Rrehilo_d = new double*[dim];
   double **Rrelolo_d = new double*[dim];
   double **Rimhihi_d = new double*[dim];
   double **Rimlohi_d = new double*[dim];
   double **Rimhilo_d = new double*[dim];
   double **Rimlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rrehihi_h[i] = new double[dim];
      Rrelohi_h[i] = new double[dim];
      Rrehilo_h[i] = new double[dim];
      Rrelolo_h[i] = new double[dim];
      Rimhihi_h[i] = new double[dim];
      Rimlohi_h[i] = new double[dim];
      Rimhilo_h[i] = new double[dim];
      Rimlolo_h[i] = new double[dim];
      Rrehihi_d[i] = new double[dim];
      Rrelohi_d[i] = new double[dim];
      Rrehilo_d[i] = new double[dim];
      Rrelolo_d[i] = new double[dim];
      Rimhihi_d[i] = new double[dim];
      Rimlohi_d[i] = new double[dim];
      Rimhilo_d[i] = new double[dim];
      Rimlolo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   cmplx4_unit_series_vector
      (dim,deg,inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
               inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehihi_d[i][j] = inputrehihi_h[i][j];
         inputrelohi_d[i][j] = inputrelohi_h[i][j];
         inputrehilo_d[i][j] = inputrehilo_h[i][j];
         inputrelolo_d[i][j] = inputrelolo_h[i][j];
         inputimhihi_d[i][j] = inputimhihi_h[i][j];
         inputimlohi_d[i][j] = inputimlohi_h[i][j];
         inputimhilo_d[i][j] = inputimhilo_h[i][j];
         inputimlolo_d[i][j] = inputimlolo_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << inputrehihi_h[i][0] << "  "
              << inputrelohi_h[i][0] << endl << "  "
              << inputrehilo_h[i][0] << "  "
              << inputrelolo_h[i][0] << endl << "  "
              << inputimhihi_h[i][0] << "  "
              << inputimlohi_h[i][0] << endl << "  "
              << inputimhilo_h[i][0] << "  "
              << inputimlolo_h[i][0] << endl;
   }
   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;

      cmplx4_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,
          accrehihi,accrelohi,accrehilo,accrelolo,
          accimhihi,accimlohi,accimhilo,accimlolo,
          inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
          inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
          outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
          outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
          outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
          outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
          funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
          funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
          funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
          funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
          jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
          jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
          jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
          jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
          rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
          rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
          rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
          rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
          urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
          urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
          urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
          urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
          solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
          solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
          solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
          solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
          Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
          Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
          Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
          Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
          Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
          Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
          Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
          Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
          workmatrehihi,workmatrelohi,workmatrehilo,workmatrelolo,
          workmatimhihi,workmatimlohi,workmatimhilo,workmatimlolo,
          workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
          workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,
          resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
          resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,vrblvl,mode);
   }
   return 0;
}
