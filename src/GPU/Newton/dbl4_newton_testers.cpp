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

   if(vrblvl > 0)
      cout << "maximum residual : "
           << *resmaxhihi << "  " << *resmaxlohi << endl
           << *resmaxhilo << "  " << *resmaxlolo << endl;

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
            cout << "rhs_h[" << i << "][" << j << "] : "
                 << rhshihi_h[i][j] << "  " << rhslohi_h[i][j] << endl << "  "
                 << rhshilo_h[i][j] << "  " << rhslolo_h[i][j] << endl;
            cout << "rhs_d[" << i << "][" << j << "] : "
                 << rhshihi_d[i][j] << "  " << rhslohi_d[i][j] << endl << "  "
                 << rhshilo_d[i][j] << "  " << rhslolo_d[i][j] << endl;
            errsum += abs(rhshihi_h[i][j] - rhshihi_d[i][j])
                    + abs(rhslohi_h[i][j] - rhslohi_d[i][j])
                    + abs(rhshilo_h[i][j] - rhshilo_d[i][j])
                    + abs(rhslolo_h[i][j] - rhslolo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }

   if(vrblvl > 0) cout << "saving the original rhs ..." << endl;

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

   if(vrblvl > 0) cout << "initializing the solution ..." << endl;

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

         cout << "maximum residual : "
              << *resmaxhihi << "  " << *resmaxlohi << "  "
              << *resmaxhilo << "  " << *resmaxlolo << endl;
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

         cout << "maximum residual : "
              << *resmaxhihi << "  " << *resmaxlohi << "  "
              << *resmaxhilo << "  " << *resmaxlolo << endl;
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
             cout << "Q_h[" << i << "][" << j << "] : "
                  << Qhihi_h[i][j] << "  " << Qlohi_h[i][j] << endl << "  "
                  << Qhilo_h[i][j] << "  " << Qlolo_h[i][j] << endl;
             cout << "Q_d[" << i << "][" << j << "] : "
                  << Qhihi_d[i][j] << "  " << Qlohi_d[i][j] << endl << "  "
                  << Qhilo_d[i][j] << "  " << Qlolo_d[i][j] << endl;
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
             cout << "R_h[" << i << "][" << j << "] : "
                  << Rhihi_h[i][j] << "  " << Rlohi_h[i][j] << endl << "  "
                  << Rhilo_h[i][j] << "  " << Rlolo_h[i][j] << endl;
             cout << "R_d[" << i << "][" << j << "] : "
                  << Rhihi_d[i][j] << "  " << Rlohi_d[i][j] << endl << "  "
                  << Rhilo_d[i][j] << "  " << Rlolo_d[i][j] << endl;
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
             cout << "urhs_h[" << i << "][" << j << "] : "
                  << urhshihi_h[i][j] << "  " << urhslohi_h[i][j]
                  << endl << "  "
                  << urhshilo_h[i][j] << "  " << urhslolo_h[i][j] << endl;
             cout << "urhs_d[" << i << "][" << j << "] : "
                  << urhshihi_d[i][j] << "  " << urhslohi_d[i][j]
                  << endl << "  "
                  << urhshilo_d[i][j] << "  " << urhslolo_d[i][j] << endl;
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
             cout << "sol_h[" << i << "][" << j << "] : "
                  << solhihi_h[i][j] << "  " << sollohi_h[i][j] << endl << "  "
                  << solhilo_h[i][j] << "  " << sollolo_h[i][j] << endl;
             cout << "sol_d[" << i << "][" << j << "] : "
                  << solhihi_d[i][j] << "  " << sollohi_d[i][j] << endl << "  "
                  << solhilo_d[i][j] << "  " << sollolo_d[i][j] << endl;
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
             errsum += abs(inputhihi_h[i][j] - inputhihi_d[i][j])
                     + abs(inputlohi_h[i][j] - inputlohi_d[i][j])
                     + abs(inputhilo_h[i][j] - inputhilo_d[i][j])
                     + abs(inputlolo_h[i][j] - inputlolo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}
