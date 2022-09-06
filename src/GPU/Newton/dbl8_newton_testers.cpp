// The file dbl8_newton_testers.cpp defines the functions with prototypes in
// the file dbl8_newton_testers.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#include "unimodular_matrices.h"
#include "random_monomials.h"
#include "octo_double_functions.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"
#include "dbl8_factorizations.h"
#include "dbl8_bals_host.h"
#include "dbl8_bals_kernels.h"
#include "dbl8_systems_host.h"
#include "dbl8_systems_kernels.h"

using namespace std;

void dbl8_unit_series_vector
 ( int dim, int deg,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo )
{
   for(int i=0; i<dim; i++)
   {
      cffhihihi[i][0] = 1.0; cfflohihi[i][0] = 0.0;
      cffhilohi[i][0] = 0.0; cfflolohi[i][0] = 0.0;
      cffhihilo[i][0] = 0.0; cfflohilo[i][0] = 0.0;
      cffhilolo[i][0] = 0.0; cfflololo[i][0] = 0.0;

      for(int j=1; j<=deg; j++)
      {
         cffhihihi[i][j] = 0.0; cfflohihi[i][j] = 0.0;
         cffhilohi[i][j] = 0.0; cfflolohi[i][j] = 0.0;
         cffhihilo[i][j] = 0.0; cfflohilo[i][j] = 0.0;
         cffhilolo[i][j] = 0.0; cfflololo[i][j] = 0.0;
      }
   }
}

void dbl8_update_series
 ( int dim, int degp1,
   double **xhihihi, double **xlohihi, double **xhilohi, double **xlolohi,
   double **xhihilo, double **xlohilo, double **xhilolo, double **xlololo,
   double **dxhihihi, double **dxlohihi, double **dxhilohi, double **dxlolohi,
   double **dxhihilo, double **dxlohilo, double **dxhilolo, double **dxlololo,
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
            cout << xhihihi[i][j] << "  " << xlohihi[i][j] << endl;
            cout << xhilohi[i][j] << "  " << xlolohi[i][j] << endl;
            cout << xhihilo[i][j] << "  " << xlohilo[i][j] << endl;
            cout << xhilolo[i][j] << "  " << xlololo[i][j] << endl;
         }
      }
   }
   // The update dx is linearized, the series x is not.
   for(int j=0; j<degp1; j++) 
      for(int i=0; i<dim; i++) // x[i][j] = x[i][j] + dx[j][i];
      {
         odf_inc(&xhihihi[i][j],&xlohihi[i][j],&xhilohi[i][j],&xlolohi[i][j],
                 &xhihilo[i][j],&xlohilo[i][j],&xhilolo[i][j],&xlololo[i][j],
                 dxhihihi[j][i],dxlohihi[j][i],dxhilohi[j][i],dxlolohi[j][i],
                 dxhihilo[j][i],dxlohilo[j][i],dxhilolo[j][i],dxlololo[j][i]);
      }

   if(vrblvl > 0)
   {
      cout << "The series after the update : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            cout << xhihihi[i][j] << "  " << xlohihi[i][j] << endl;
            cout << xhilohi[i][j] << "  " << xlolohi[i][j] << endl;
            cout << xhihilo[i][j] << "  " << xlohilo[i][j] << endl;
            cout << xhilolo[i][j] << "  " << xlololo[i][j] << endl;
         }
      }
   }
}

void dbl8_newton_lustep
 ( int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi,
   double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo,
   double *acchilolo, double *acclololo,
   double **inputhihihi, double **inputlohihi,
   double **inputhilohi, double **inputlolohi,
   double **inputhihilo, double **inputlohilo,
   double **inputhilolo, double **inputlololo,
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
   double ***jacvalhilolo, double ***jacvallololo,
   double **rhshihihi, double **rhslohihi,
   double **rhshilohi, double **rhslolohi,
   double **rhshihilo, double **rhslohilo,
   double **rhshilolo, double **rhslololo,
   double **solhihihi, double **sollohihi,
   double **solhilohi, double **sollolohi,
   double **solhihilo, double **sollohilo,
   double **solhilolo, double **sollololo,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **workrhshihihi, double **workrhslohihi,
   double **workrhshilohi, double **workrhslolohi,
   double **workrhshihilo, double **workrhslohilo,
   double **workrhshilolo, double **workrhslololo,
   double **resvechihihi, double **resveclohihi,
   double **resvechilohi, double **resveclolohi,
   double **resvechihilo, double **resveclohilo,
   double **resvechilolo, double **resveclololo,
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int *ipvt, int vrblvl )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   for(int i=0; i<dim; i++)
   {
      cffhihihi[i][0] = 1.0; cfflohihi[i][0] = 0.0;
      cffhilohi[i][0] = 0.0; cfflolohi[i][0] = 0.0;
      cffhihilo[i][0] = 0.0; cfflohilo[i][0] = 0.0;
      cffhilolo[i][0] = 0.0; cfflololo[i][0] = 0.0;

      for(int j=1; j<degp1; j++)
      {
         cffhihihi[i][j] = 0.0; cfflohihi[i][j] = 0.0;
         cffhilohi[i][j] = 0.0; cfflolohi[i][j] = 0.0;
         cffhihilo[i][j] = 0.0; cfflohilo[i][j] = 0.0;
         cffhilolo[i][j] = 0.0; cfflololo[i][j] = 0.0;
      }
   }
   CPU_dbl8_evaluate_monomials
      (dim,deg,nvr,idx,exp,nbrfac,expfac,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       acchihihi,acclohihi,acchilohi,acclolohi,
       acchihilo,acclohilo,acchilolo,acclololo,
       inputhihihi, inputlohihi, inputhilohi, inputlolohi,
       inputhihilo, inputlohilo, inputhilolo, inputlololo,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,vrblvl);

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhihihi[i][j][k] = 0.0; jacvallohihi[i][j][k] = 0.0;
            jacvalhilohi[i][j][k] = 0.0; jacvallolohi[i][j][k] = 0.0;
            jacvalhihilo[i][j][k] = 0.0; jacvallohilo[i][j][k] = 0.0;
            jacvalhilolo[i][j][k] = 0.0; jacvallololo[i][j][k] = 0.0;
         }

   dbl8_linearize_evaldiff_output
      (dim,degp1,nvr,idx,
       outputhihihi,outputlohihi,outputhilohi,outputlolohi,
       outputhihilo,outputlohilo,outputhilolo,outputlololo,
       funvalhihihi,funvallohihi,funvalhilohi,funvallolohi,
       funvalhihilo,funvallohilo,funvalhilolo,funvallololo,
       rhshihihi,rhslohihi,rhshilohi,rhslolohi,
       rhshihilo,rhslohilo,rhshilolo,rhslololo,
       jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
       jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,vrblvl);

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         workrhshihihi[i][j] = rhshihihi[i][j];
         workrhslohihi[i][j] = rhslohihi[i][j];
         workrhshilohi[i][j] = rhshilohi[i][j];
         workrhslolohi[i][j] = rhslolohi[i][j];
         workrhshihilo[i][j] = rhshihilo[i][j];
         workrhslohilo[i][j] = rhslohilo[i][j];
         workrhshilolo[i][j] = rhshilolo[i][j];
         workrhslololo[i][j] = rhslololo[i][j];
      }
   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhihihi[i][j] = 0.0; sollohihi[i][j] = 0.0;
         solhilohi[i][j] = 0.0; sollolohi[i][j] = 0.0;
         solhihilo[i][j] = 0.0; sollohilo[i][j] = 0.0;
         solhilolo[i][j] = 0.0; sollololo[i][j] = 0.0;
      }
 
   CPU_dbl8_lusb_solve
      (dim,degp1,
       jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
       jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
       workrhshihihi,workrhslohihi,workrhshilohi,workrhslolohi,
       workrhshihilo,workrhslohilo,workrhshilolo,workrhslololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,
       workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
       workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
       workvechihihi,workveclohihi,workvechilohi,workveclolohi,
       workvechihilo,workveclohilo,workvechilolo,workveclololo,
       ipvt,0); // vrblvl);

   CPU_dbl8_linear_residue
      (dim,degp1,
       jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
       jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
       rhshihihi,rhslohihi,rhshilohi,rhslolohi,
       rhshihilo,rhslohilo,rhshilolo,rhslololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,
       resvechihihi,resveclohihi,resvechilohi,resveclolohi,
       resvechihilo,resveclohilo,resvechilolo,resveclololo,
       resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
       resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);

   if(vrblvl > 0)
      cout << "maximum residual : "
           << *resmaxhihihi << "  " << *resmaxlohihi << endl
           << *resmaxhilohi << "  " << *resmaxlolohi << endl
           << *resmaxhihilo << "  " << *resmaxlohilo << endl
           << *resmaxhilolo << "  " << *resmaxlololo << endl;

   dbl8_update_series
      (dim,degp1,
       inputhihihi,inputlohihi,inputhilohi,inputlolohi,
       inputhihilo,inputlohilo,inputhilolo,inputlololo,
       solhihihi,sollohihi,solhilohi,sollolohi,
       solhihilo,sollohilo,solhilolo,sollololo,vrblvl);
}

void dbl8_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double **cffhihihi, double **cfflohihi,
   double **cffhilohi, double **cfflolohi,
   double **cffhihilo, double **cfflohilo,
   double **cffhilolo, double **cfflololo,
   double *acchihihi, double *acclohihi, double *acchilohi, double *acclolohi,
   double *acchihilo, double *acclohilo, double *acchilolo, double *acclololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
   double ***outputhihihi_h, double ***outputlohihi_h,
   double ***outputhilohi_h, double ***outputlolohi_h,
   double ***outputhihilo_h, double ***outputlohilo_h,
   double ***outputhilolo_h, double ***outputlololo_h,
   double ***outputhihihi_d, double ***outputlohihi_d,
   double ***outputhilohi_d, double ***outputlolohi_d,
   double ***outputhihilo_d, double ***outputlohilo_d,
   double ***outputhilolo_d, double ***outputlololo_d,
   double **funvalhihihi_h, double **funvallohihi_h,
   double **funvalhilohi_h, double **funvallolohi_h,
   double **funvalhihilo_h, double **funvallohilo_h,
   double **funvalhilolo_h, double **funvallololo_h,
   double **funvalhihihi_d, double **funvallohihi_d,
   double **funvalhilohi_d, double **funvallolohi_d,
   double **funvalhihilo_d, double **funvallohilo_d,
   double **funvalhilolo_d, double **funvallololo_d,
   double ***jacvalhihihi_h, double ***jacvallohihi_h,
   double ***jacvalhilohi_h, double ***jacvallolohi_h,
   double ***jacvalhihilo_h, double ***jacvallohilo_h,
   double ***jacvalhilolo_h, double ***jacvallololo_h,
   double ***jacvalhihihi_d, double ***jacvallohihi_d,
   double ***jacvalhilohi_d, double ***jacvallolohi_d,
   double ***jacvalhihilo_d, double ***jacvallohilo_d,
   double ***jacvalhilolo_d, double ***jacvallololo_d,
   double **rhshihihi_h, double **rhslohihi_h,
   double **rhshilohi_h, double **rhslolohi_h,
   double **rhshihilo_h, double **rhslohilo_h,
   double **rhshilolo_h, double **rhslololo_h,
   double **rhshihihi_d, double **rhslohihi_d,
   double **rhshilohi_d, double **rhslolohi_d,
   double **rhshihilo_d, double **rhslohilo_d,
   double **rhshilolo_d, double **rhslololo_d,
   double **urhshihihi_h, double **urhslohihi_h,
   double **urhshilohi_h, double **urhslolohi_h,
   double **urhshihilo_h, double **urhslohilo_h,
   double **urhshilolo_h, double **urhslololo_h,
   double **urhshihihi_d, double **urhslohihi_d,
   double **urhshilohi_d, double **urhslolohi_d,
   double **urhshihilo_d, double **urhslohilo_d,
   double **urhshilolo_d, double **urhslololo_d,
   double **solhihihi_h, double **sollohihi_h,
   double **solhilohi_h, double **sollolohi_h,
   double **solhihilo_h, double **sollohilo_h,
   double **solhilolo_h, double **sollololo_h,
   double **solhihihi_d, double **sollohihi_d,
   double **solhilohi_d, double **sollolohi_d,
   double **solhihilo_d, double **sollohilo_d,
   double **solhilolo_d, double **sollololo_d,
   double **Qhihihi_h, double **Qlohihi_h,
   double **Qhilohi_h, double **Qlolohi_h,
   double **Qhihilo_h, double **Qlohilo_h,
   double **Qhilolo_h, double **Qlololo_h,
   double **Qhihihi_d, double **Qlohihi_d,
   double **Qhilohi_d, double **Qlolohi_d,
   double **Qhihilo_d, double **Qlohilo_d,
   double **Qhilolo_d, double **Qlololo_d,
   double **Rhihihi_h, double **Rlohihi_h,
   double **Rhilohi_h, double **Rlolohi_h,
   double **Rhihilo_h, double **Rlohilo_h,
   double **Rhilolo_h, double **Rlololo_h,
   double **Rhihihi_d, double **Rlohihi_d,
   double **Rhilohi_d, double **Rlolohi_d,
   double **Rhihilo_d, double **Rlohilo_d,
   double **Rhilolo_d, double **Rlololo_d,
   double **workmathihihi, double **workmatlohihi,
   double **workmathilohi, double **workmatlolohi,
   double **workmathihilo, double **workmatlohilo,
   double **workmathilolo, double **workmatlololo,
   double *workvechihihi, double *workveclohihi,
   double *workvechilohi, double *workveclolohi,
   double *workvechihilo, double *workveclohilo,
   double *workvechilolo, double *workveclololo,
   double **resvechihihi, double **resveclohihi, 
   double **resvechilohi, double **resveclolohi, 
   double **resvechihilo, double **resveclohilo, 
   double **resvechilolo, double **resveclololo, 
   double *resmaxhihihi, double *resmaxlohihi,
   double *resmaxhilohi, double *resmaxlolohi,
   double *resmaxhihilo, double *resmaxlohilo,
   double *resmaxhilolo, double *resmaxlololo, int vrblvl, int mode )
{
   const int degp1 = deg+1;

   // The series coefficients accumulate common factors,
   // initially the coefficients are set to one.
   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         cffhihihi[i][0] = 1.0; cfflohihi[i][0] = 0.0;
         cffhilohi[i][0] = 0.0; cfflolohi[i][0] = 0.0;
         cffhihilo[i][0] = 0.0; cfflohilo[i][0] = 0.0;
         cffhilolo[i][0] = 0.0; cfflololo[i][0] = 0.0;

         for(int j=1; j<degp1; j++)
         {
            cffhihihi[i][j] = 0.0; cfflohihi[i][j] = 0.0;
            cffhilohi[i][j] = 0.0; cfflolohi[i][j] = 0.0;
            cffhihilo[i][j] = 0.0; cfflohilo[i][j] = 0.0;
            cffhilolo[i][j] = 0.0; cfflololo[i][j] = 0.0;
         }
      }
      CPU_dbl8_evaluate_monomials
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         cffhihihi[i][0] = 1.0; cfflohihi[i][0] = 0.0;
         cffhilohi[i][0] = 0.0; cfflolohi[i][0] = 0.0;
         cffhihilo[i][0] = 0.0; cfflohilo[i][0] = 0.0;
         cffhilolo[i][0] = 0.0; cfflololo[i][0] = 0.0;

         for(int j=1; j<degp1; j++)
         {
            cffhihihi[i][j] = 0.0; cfflohihi[i][j] = 0.0;
            cffhilohi[i][j] = 0.0; cfflolohi[i][j] = 0.0;
            cffhihilo[i][j] = 0.0; cfflohilo[i][j] = 0.0;
            cffhilolo[i][j] = 0.0; cfflololo[i][j] = 0.0;
         }
      }
      GPU_dbl8_evaluate_monomials
         (dim,deg,szt,nbt,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
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
             cout << "output_h[" << k << "][" << i << "][" << j << "] : "
                  << outputhihihi_h[k][i][j] << "  "
                  << outputlohihi_h[k][i][j] << endl << "  "
                  << outputhilohi_h[k][i][j] << "  "
                  << outputlolohi_h[k][i][j] << endl << "  "
                  << outputhihilo_h[k][i][j] << "  "
                  << outputlohilo_h[k][i][j] << endl << "  "
                  << outputhilolo_h[k][i][j] << "  "
                  << outputlololo_h[k][i][j] << endl;
             cout << "output_d[" << k << "][" << i << "][" << j << "] : "
                  << outputhihihi_d[k][i][j] << "  "
                  << outputlohihi_d[k][i][j] << endl << "  "
                  << outputhilohi_d[k][i][j] << "  "
                  << outputlolohi_d[k][i][j] << endl << "  "
                  << outputhihilo_d[k][i][j] << "  "
                  << outputlohilo_d[k][i][j] << endl << "  "
                  << outputhilolo_d[k][i][j] << "  "
                  << outputlololo_d[k][i][j] << endl;
             errsum += abs(outputhihihi_h[k][i][j] - outputhihihi_d[k][i][j])
                     + abs(outputlohihi_h[k][i][j] - outputlohihi_d[k][i][j])
                     + abs(outputhilohi_h[k][i][j] - outputhilohi_d[k][i][j])
                     + abs(outputlolohi_h[k][i][j] - outputlolohi_d[k][i][j])
                     + abs(outputhihilo_h[k][i][j] - outputhihilo_d[k][i][j])
                     + abs(outputlohilo_h[k][i][j] - outputlohilo_d[k][i][j])
                     + abs(outputhilolo_h[k][i][j] - outputhilolo_d[k][i][j])
                     + abs(outputlololo_h[k][i][j] - outputlololo_d[k][i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }

   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
      for(int j=0; j<dim; j++) 
         for(int k=0; k<dim; k++)
         {
            jacvalhihihi_h[i][j][k] = 0.0; jacvallohihi_h[i][j][k] = 0.0;
            jacvalhilohi_h[i][j][k] = 0.0; jacvallolohi_h[i][j][k] = 0.0;
            jacvalhihilo_h[i][j][k] = 0.0; jacvallohilo_h[i][j][k] = 0.0;
            jacvalhilolo_h[i][j][k] = 0.0; jacvallololo_h[i][j][k] = 0.0;
            jacvalhihihi_d[i][j][k] = 0.0; jacvallohihi_d[i][j][k] = 0.0;
            jacvalhilohi_d[i][j][k] = 0.0; jacvallolohi_d[i][j][k] = 0.0;
            jacvalhihilo_d[i][j][k] = 0.0; jacvallohilo_d[i][j][k] = 0.0;
            jacvalhilolo_d[i][j][k] = 0.0; jacvallololo_d[i][j][k] = 0.0;
         }

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
      dbl8_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
          funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
          rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
          rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,vrblvl);

   if((mode == 1) || (mode == 2))
      dbl8_linearize_evaldiff_output
         (dim,degp1,nvr,idx,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
          funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
          funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
          rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
          rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,vrblvl);

   if((vrblvl > 0) && (mode == 2))
   {
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
         {
            cout << "funval_h[" << i << "][" << j << "] : "
                 << funvalhihihi_h[i][j] << "  "
                 << funvallohihi_h[i][j] << endl << "  "
                 << funvalhilohi_h[i][j] << "  "
                 << funvallolohi_h[i][j] << endl << "  "
                 << funvalhihilo_h[i][j] << "  "
                 << funvallohilo_h[i][j] << endl << "  "
                 << funvalhilolo_h[i][j] << "  "
                 << funvallololo_h[i][j] << endl;
            cout << "funval_d[" << i << "][" << j << "] : "
                 << funvalhihihi_d[i][j] << "  "
                 << funvallohihi_d[i][j] << endl << "  "
                 << funvalhilohi_d[i][j] << "  "
                 << funvallolohi_d[i][j] << endl << "  "
                 << funvalhihilo_d[i][j] << "  "
                 << funvallohilo_d[i][j] << endl << "  "
                 << funvalhilolo_d[i][j] << "  "
                 << funvallololo_d[i][j] << endl;
            errsum += abs(funvalhihihi_h[i][j] - funvalhihihi_d[i][j])
                    + abs(funvallohihi_h[i][j] - funvallohihi_d[i][j])
                    + abs(funvalhilohi_h[i][j] - funvalhilohi_d[i][j])
                    + abs(funvallolohi_h[i][j] - funvallolohi_d[i][j])
                    + abs(funvalhihilo_h[i][j] - funvalhihilo_d[i][j])
                    + abs(funvallohilo_h[i][j] - funvallohilo_d[i][j])
                    + abs(funvalhilolo_h[i][j] - funvalhilolo_d[i][j])
                    + abs(funvallololo_h[i][j] - funvallololo_d[i][j]);
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
                    << jacvalhihihi_h[i][j][k] << "  "
                    << jacvallohihi_h[i][j][k] << endl << "  "
                    << jacvalhilohi_h[i][j][k] << "  "
                    << jacvallolohi_h[i][j][k] << endl << "  "
                    << jacvalhihilo_h[i][j][k] << "  "
                    << jacvallohilo_h[i][j][k] << endl << "  "
                    << jacvalhilolo_h[i][j][k] << "  "
                    << jacvallololo_h[i][j][k] << endl;
               cout << "jacval_d[" << i << "][" << j << "][" << k << "] : "
                    << jacvalhihihi_d[i][j][k] << "  "
                    << jacvallohihi_d[i][j][k] << endl << "  "
                    << jacvalhilohi_d[i][j][k] << "  "
                    << jacvallolohi_d[i][j][k] << endl << "  "
                    << jacvalhihilo_d[i][j][k] << "  "
                    << jacvallohilo_d[i][j][k] << endl << "  "
                    << jacvalhilolo_d[i][j][k] << "  "
                    << jacvallololo_d[i][j][k] << endl;
               errsum += abs(jacvalhihihi_h[i][j][k] - jacvalhihihi_d[i][j][k])
                       + abs(jacvallohihi_h[i][j][k] - jacvallohihi_d[i][j][k])
                       + abs(jacvalhilohi_h[i][j][k] - jacvalhilohi_d[i][j][k])
                       + abs(jacvallolohi_h[i][j][k] - jacvallolohi_d[i][j][k])
                       + abs(jacvalhihilo_h[i][j][k] - jacvalhihilo_d[i][j][k])
                       + abs(jacvallohilo_h[i][j][k] - jacvallohilo_d[i][j][k])
                       + abs(jacvalhilolo_h[i][j][k] - jacvalhilolo_d[i][j][k])
                       + abs(jacvallololo_h[i][j][k] - jacvallololo_d[i][j][k]);
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
                 << rhshihihi_h[i][j] << "  "
                 << rhslohihi_h[i][j] << endl << "  "
                 << rhshilohi_h[i][j] << "  "
                 << rhslolohi_h[i][j] << endl << "  "
                 << rhshihilo_h[i][j] << "  "
                 << rhslohilo_h[i][j] << endl << "  "
                 << rhshilolo_h[i][j] << "  "
                 << rhslololo_h[i][j] << endl;
            cout << "rhs_d[" << i << "][" << j << "] : "
                 << rhshihihi_d[i][j] << "  "
                 << rhslohihi_d[i][j] << endl << "  "
                 << rhshilohi_d[i][j] << "  "
                 << rhslolohi_d[i][j] << endl << "  "
                 << rhshihilo_d[i][j] << "  "
                 << rhslohilo_d[i][j] << endl << "  "
                 << rhshilolo_d[i][j] << "  "
                 << rhslololo_d[i][j] << endl;
            errsum += abs(rhshihihi_h[i][j] - rhshihihi_d[i][j])
                    + abs(rhslohihi_h[i][j] - rhslohihi_d[i][j])
                    + abs(rhshilohi_h[i][j] - rhshilohi_d[i][j])
                    + abs(rhslolohi_h[i][j] - rhslolohi_d[i][j])
                    + abs(rhshihilo_h[i][j] - rhshihilo_d[i][j])
                    + abs(rhslohilo_h[i][j] - rhslohilo_d[i][j])
                    + abs(rhshilolo_h[i][j] - rhshilolo_d[i][j])
                    + abs(rhslololo_h[i][j] - rhslololo_d[i][j]);
         }
      }
      cout << "sum of errors : " << errsum << endl;
   }

   if(vrblvl > 0) cout << "saving the original rhs ..." << endl;

   for(int i=0; i<degp1; i++) // save original rhs for residual
      for(int j=0; j<dim; j++)
      {
         urhshihihi_h[i][j] = rhshihihi_h[i][j];
         urhslohihi_h[i][j] = rhslohihi_h[i][j];
         urhshilohi_h[i][j] = rhshilohi_h[i][j];
         urhslolohi_h[i][j] = rhslolohi_h[i][j];
         urhshihilo_h[i][j] = rhshihilo_h[i][j];
         urhslohilo_h[i][j] = rhslohilo_h[i][j];
         urhshilolo_h[i][j] = rhshilolo_h[i][j];
         urhslololo_h[i][j] = rhslololo_h[i][j];
         urhshihihi_d[i][j] = rhshihihi_d[i][j];
         urhslohihi_d[i][j] = rhslohihi_d[i][j];
         urhshilohi_d[i][j] = rhshilohi_d[i][j];
         urhslolohi_d[i][j] = rhslolohi_d[i][j];
         urhshihilo_d[i][j] = rhshihilo_d[i][j];
         urhslohilo_d[i][j] = rhslohilo_d[i][j];
         urhshilolo_d[i][j] = rhshilolo_d[i][j];
         urhslololo_d[i][j] = rhslololo_d[i][j];
      }

   if(vrblvl > 0) cout << "initializing the solution ..." << endl;

   for(int i=0; i<degp1; i++) // initialize the solution to zero
      for(int j=0; j<dim; j++)
      {
         solhihihi_h[i][j] = 0.0; sollohihi_h[i][j] = 0.0;
         solhilohi_h[i][j] = 0.0; sollolohi_h[i][j] = 0.0;
         solhihilo_h[i][j] = 0.0; sollohilo_h[i][j] = 0.0;
         solhilolo_h[i][j] = 0.0; sollololo_h[i][j] = 0.0;
         solhihihi_d[i][j] = 0.0; sollohihi_d[i][j] = 0.0;
         solhilohi_d[i][j] = 0.0; sollolohi_d[i][j] = 0.0;
         solhihilo_d[i][j] = 0.0; sollohilo_d[i][j] = 0.0;
         solhilolo_d[i][j] = 0.0; sollololo_d[i][j] = 0.0;
      }
 
   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling CPU_dbl8_qrbs_solve ..." << endl;
      CPU_dbl8_qrbs_solve
         (dim,degp1,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
          urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
          urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl4_linear_residue ..." << endl;

         CPU_dbl8_linear_residue
            (dim,degp1,
             jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
             jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
             rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
             rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
             solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
             solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
             resvechihihi,resveclohihi,resvechilohi,resveclolohi,
             resvechihilo,resveclohilo,resvechilolo,resveclololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);

         cout << "maximum residual : "
              << *resmaxhihihi << "  " << *resmaxlohihi << "  "
              << *resmaxhilohi << "  " << *resmaxlolohi << "  "
              << *resmaxhihilo << "  " << *resmaxlohilo << "  "
              << *resmaxhilolo << "  " << *resmaxlololo << endl;
      }
      dbl8_update_series
         (dim,degp1,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0) cout << "calling GPU_dbl4_bals_solve ..." << endl;

      GPU_dbl8_bals_solve
         (dim,degp1,szt,nbt,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
          urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl4_linear_residue ..." << endl;

         CPU_dbl8_linear_residue
            (dim,degp1,
             jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
             jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
             rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
             rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
             solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
             solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
             resvechihihi,resveclohihi,resvechilohi,resveclolohi,
             resvechihilo,resveclohilo,resvechilolo,resveclololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,vrblvl);

         cout << "maximum residual : "
              << *resmaxhihihi << "  " << *resmaxlohihi << "  "
              << *resmaxhilohi << "  " << *resmaxlolohi << "  "
              << *resmaxhihilo << "  " << *resmaxlohilo << "  "
              << *resmaxhilolo << "  " << *resmaxlololo << endl;
      }
      dbl8_update_series
         (dim,degp1,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "Q_h[" << i << "][" << j << "] : "
                  << Qhihihi_h[i][j] << "  "
                  << Qlohihi_h[i][j] << endl << "  "
                  << Qhilohi_h[i][j] << "  "
                  << Qlolohi_h[i][j] << endl << "  "
                  << Qhihilo_h[i][j] << "  "
                  << Qlohilo_h[i][j] << endl << "  "
                  << Qhilolo_h[i][j] << "  "
                  << Qlololo_h[i][j] << endl;
             cout << "Q_d[" << i << "][" << j << "] : "
                  << Qhihihi_d[i][j] << "  "
                  << Qlohihi_d[i][j] << endl << "  "
                  << Qhilohi_d[i][j] << "  "
                  << Qlolohi_d[i][j] << endl << "  "
                  << Qhihilo_d[i][j] << "  "
                  << Qlohilo_d[i][j] << endl << "  "
                  << Qhilolo_d[i][j] << "  "
                  << Qlololo_d[i][j] << endl;
             errsum += abs(Qhihihi_h[i][j] - Qhihihi_d[i][j])
                     + abs(Qlohihi_h[i][j] - Qlohihi_d[i][j])
                     + abs(Qhilohi_h[i][j] - Qhilohi_d[i][j])
                     + abs(Qlolohi_h[i][j] - Qlolohi_d[i][j])
                     + abs(Qhihilo_h[i][j] - Qhihilo_d[i][j])
                     + abs(Qlohilo_h[i][j] - Qlohilo_d[i][j])
                     + abs(Qhilolo_h[i][j] - Qhilolo_d[i][j])
                     + abs(Qlololo_h[i][j] - Qlololo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "R_h[" << i << "][" << j << "] : "
                  << Rhihihi_h[i][j] << "  "
                  << Rlohihi_h[i][j] << endl << "  "
                  << Rhilohi_h[i][j] << "  "
                  << Rlolohi_h[i][j] << endl << "  "
                  << Rhihilo_h[i][j] << "  "
                  << Rlohilo_h[i][j] << endl << "  "
                  << Rhilolo_h[i][j] << "  "
                  << Rlololo_h[i][j] << endl;
             cout << "R_d[" << i << "][" << j << "] : "
                  << Rhihihi_d[i][j] << "  "
                  << Rlohihi_d[i][j] << endl << "  "
                  << Rhilohi_d[i][j] << "  "
                  << Rlolohi_d[i][j] << endl << "  "
                  << Rhihilo_d[i][j] << "  "
                  << Rlohilo_d[i][j] << endl << "  "
                  << Rhilolo_d[i][j] << "  "
                  << Rlololo_d[i][j] << endl;
             errsum += abs(Rhihihi_h[i][j] - Rhihihi_d[i][j])
                     + abs(Rlohihi_h[i][j] - Rlohihi_d[i][j])
                     + abs(Rhilohi_h[i][j] - Rhilohi_d[i][j])
                     + abs(Rlolohi_h[i][j] - Rlolohi_d[i][j])
                     + abs(Rhihilo_h[i][j] - Rhihilo_d[i][j])
                     + abs(Rlohilo_h[i][j] - Rlohilo_d[i][j])
                     + abs(Rhilolo_h[i][j] - Rhilolo_d[i][j])
                     + abs(Rlololo_h[i][j] - Rlololo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "urhs_h[" << i << "][" << j << "] : "
                  << urhshihihi_h[i][j] << "  "
                  << urhslohihi_h[i][j] << endl << "  "
                  << urhshilohi_h[i][j] << "  "
                  << urhslolohi_h[i][j] << endl << "  "
                  << urhshihilo_h[i][j] << "  "
                  << urhslohilo_h[i][j] << endl << "  "
                  << urhshilolo_h[i][j] << "  "
                  << urhslololo_h[i][j] << endl;
             cout << "urhs_d[" << i << "][" << j << "] : "
                  << urhshihihi_d[i][j] << "  "
                  << urhslohihi_d[i][j] << endl << "  "
                  << urhshilohi_d[i][j] << "  "
                  << urhslolohi_d[i][j] << endl << "  "
                  << urhshihilo_d[i][j] << "  "
                  << urhslohilo_d[i][j] << endl << "  "
                  << urhshilolo_d[i][j] << "  "
                  << urhslololo_d[i][j] << endl;
             errsum += abs(urhshihihi_h[i][j] - urhshihihi_d[i][j])
                     + abs(urhslohihi_h[i][j] - urhslohihi_d[i][j])
                     + abs(urhshilohi_h[i][j] - urhshilohi_d[i][j])
                     + abs(urhslolohi_h[i][j] - urhslolohi_d[i][j])
                     + abs(urhshihilo_h[i][j] - urhshihilo_d[i][j])
                     + abs(urhslohilo_h[i][j] - urhslohilo_d[i][j])
                     + abs(urhshilolo_h[i][j] - urhshilolo_d[i][j])
                     + abs(urhslololo_h[i][j] - urhslololo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
             cout << "sol_h[" << i << "][" << j << "] : "
                  << solhihihi_h[i][j] << "  "
                  << sollohihi_h[i][j] << endl << "  "
                  << solhilohi_h[i][j] << "  "
                  << sollolohi_h[i][j] << endl << "  "
                  << solhihilo_h[i][j] << "  "
                  << sollohilo_h[i][j] << endl << "  "
                  << solhilolo_h[i][j] << "  "
                  << sollololo_h[i][j] << endl;
             cout << "sol_d[" << i << "][" << j << "] : "
                  << solhihihi_d[i][j] << "  "
                  << sollohihi_d[i][j] << endl << "  "
                  << solhilohi_d[i][j] << "  "
                  << sollolohi_d[i][j] << endl << "  "
                  << solhihilo_d[i][j] << "  "
                  << sollohilo_d[i][j] << endl << "  "
                  << solhilolo_d[i][j] << "  "
                  << sollololo_d[i][j] << endl;
             errsum += abs(solhihihi_h[i][j] - solhihihi_d[i][j])
                     + abs(sollohihi_h[i][j] - sollohihi_d[i][j])
                     + abs(solhilohi_h[i][j] - solhilohi_d[i][j])
                     + abs(sollolohi_h[i][j] - sollolohi_d[i][j])
                     + abs(solhihilo_h[i][j] - solhihilo_d[i][j])
                     + abs(sollohilo_h[i][j] - sollohilo_d[i][j])
                     + abs(solhilolo_h[i][j] - solhilolo_d[i][j])
                     + abs(sollololo_h[i][j] - sollololo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
      errsum = 0.0;
      cout << "comparing CPU with GPU series ... " << endl;
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
             cout << "input_h[" << i << "][" << j << "] : "
                  << inputhihihi_h[i][j] << "  "
                  << inputlohihi_h[i][j] << endl << "  "
                  << inputhilohi_h[i][j] << "  "
                  << inputlolohi_h[i][j] << endl << "  "
                  << inputhihilo_h[i][j] << "  "
                  << inputlohilo_h[i][j] << endl << "  "
                  << inputhilolo_h[i][j] << "  "
                  << inputlololo_h[i][j] << endl;
             cout << "input_d[" << i << "][" << j << "] : "
                  << inputhihihi_d[i][j] << "  "
                  << inputlohihi_d[i][j] << endl << "  "
                  << inputhilohi_d[i][j] << "  "
                  << inputlolohi_d[i][j] << endl << "  "
                  << inputhihilo_d[i][j] << "  "
                  << inputlohilo_d[i][j] << endl << "  "
                  << inputhilolo_d[i][j] << "  "
                  << inputlololo_d[i][j] << endl;
             errsum += abs(inputhihihi_h[i][j] - inputhihihi_d[i][j])
                     + abs(inputlohihi_h[i][j] - inputlohihi_d[i][j])
                     + abs(inputhilohi_h[i][j] - inputhilohi_d[i][j])
                     + abs(inputlolohi_h[i][j] - inputlolohi_d[i][j])
                     + abs(inputhihilo_h[i][j] - inputhihilo_d[i][j])
                     + abs(inputlohilo_h[i][j] - inputlohilo_d[i][j])
                     + abs(inputhilolo_h[i][j] - inputhilolo_d[i][j])
                     + abs(inputlololo_h[i][j] - inputlololo_d[i][j]);
         }
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl8_real_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhihihi_h = new double*[dim];
   double **inputlohihi_h = new double*[dim];
   double **inputhilohi_h = new double*[dim];
   double **inputlolohi_h = new double*[dim];
   double **inputhihilo_h = new double*[dim];
   double **inputlohilo_h = new double*[dim];
   double **inputhilolo_h = new double*[dim];
   double **inputlololo_h = new double*[dim];
   double **inputhihihi_d = new double*[dim];
   double **inputlohihi_d = new double*[dim];
   double **inputhilohi_d = new double*[dim];
   double **inputlolohi_d = new double*[dim];
   double **inputhihilo_d = new double*[dim];
   double **inputlohilo_d = new double*[dim];
   double **inputhilolo_d = new double*[dim];
   double **inputlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhihihi_h[i] = new double[degp1];
      inputlohihi_h[i] = new double[degp1];
      inputhilohi_h[i] = new double[degp1];
      inputlolohi_h[i] = new double[degp1];
      inputhihilo_h[i] = new double[degp1];
      inputlohilo_h[i] = new double[degp1];
      inputhilolo_h[i] = new double[degp1];
      inputlololo_h[i] = new double[degp1];
      inputhihihi_d[i] = new double[degp1];
      inputlohihi_d[i] = new double[degp1];
      inputhilohi_d[i] = new double[degp1];
      inputlolohi_d[i] = new double[degp1];
      inputhihilo_d[i] = new double[degp1];
      inputlohilo_d[i] = new double[degp1];
      inputhilolo_d[i] = new double[degp1];
      inputlololo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double *acchihihi = new double[degp1]; // accumulated power series
   double *acclohihi = new double[degp1];
   double *acchilohi = new double[degp1];
   double *acclolohi = new double[degp1];
   double *acchihilo = new double[degp1];
   double *acclohilo = new double[degp1];
   double *acchilolo = new double[degp1];
   double *acclololo = new double[degp1];
   double **cffhihihi = new double*[dim]; // the coefficients of monomials
   double **cfflohihi = new double*[dim];
   double **cffhilohi = new double*[dim];
   double **cfflolohi = new double*[dim];
   double **cffhihilo = new double*[dim];
   double **cfflohilo = new double*[dim];
   double **cffhilolo = new double*[dim];
   double **cfflololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double[degp1];
      cfflohihi[i] = new double[degp1];
      cffhilohi[i] = new double[degp1];
      cfflolohi[i] = new double[degp1];
      cffhihilo[i] = new double[degp1];
      cfflohilo[i] = new double[degp1];
      cffhilolo[i] = new double[degp1];
      cfflololo[i] = new double[degp1];
   }
   double ***outputhihihi_h = new double**[dim];
   double ***outputlohihi_h = new double**[dim];
   double ***outputhilohi_h = new double**[dim];
   double ***outputlolohi_h = new double**[dim];
   double ***outputhihilo_h = new double**[dim];
   double ***outputlohilo_h = new double**[dim];
   double ***outputhilolo_h = new double**[dim];
   double ***outputlololo_h = new double**[dim];
   double ***outputhihihi_d = new double**[dim];
   double ***outputlohihi_d = new double**[dim];
   double ***outputhilohi_d = new double**[dim];
   double ***outputlolohi_d = new double**[dim];
   double ***outputhihilo_d = new double**[dim];
   double ***outputlohilo_d = new double**[dim];
   double ***outputhilolo_d = new double**[dim];
   double ***outputlololo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihihi_h[i] = new double*[dim+1];
      outputlohihi_h[i] = new double*[dim+1];
      outputhilohi_h[i] = new double*[dim+1];
      outputlolohi_h[i] = new double*[dim+1];
      outputhihilo_h[i] = new double*[dim+1];
      outputlohilo_h[i] = new double*[dim+1];
      outputhilolo_h[i] = new double*[dim+1];
      outputlololo_h[i] = new double*[dim+1];
      outputhihihi_d[i] = new double*[dim+1];
      outputlohihi_d[i] = new double*[dim+1];
      outputhilohi_d[i] = new double*[dim+1];
      outputlolohi_d[i] = new double*[dim+1];
      outputhihilo_d[i] = new double*[dim+1];
      outputlohilo_d[i] = new double*[dim+1];
      outputhilolo_d[i] = new double*[dim+1];
      outputlololo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihihi_h[i][j] = new double[degp1];
         outputlohihi_h[i][j] = new double[degp1];
         outputhilohi_h[i][j] = new double[degp1];
         outputlolohi_h[i][j] = new double[degp1];
         outputhihilo_h[i][j] = new double[degp1];
         outputlohilo_h[i][j] = new double[degp1];
         outputhilolo_h[i][j] = new double[degp1];
         outputlololo_h[i][j] = new double[degp1];
         outputhihihi_d[i][j] = new double[degp1];
         outputlohihi_d[i][j] = new double[degp1];
         outputhilohi_d[i][j] = new double[degp1];
         outputlolohi_d[i][j] = new double[degp1];
         outputhihilo_d[i][j] = new double[degp1];
         outputlohilo_d[i][j] = new double[degp1];
         outputhilolo_d[i][j] = new double[degp1];
         outputlololo_d[i][j] = new double[degp1];
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhihihi_h = new double*[dim];
   double **funvallohihi_h = new double*[dim];
   double **funvalhilohi_h = new double*[dim];
   double **funvallolohi_h = new double*[dim];
   double **funvalhihilo_h = new double*[dim];
   double **funvallohilo_h = new double*[dim];
   double **funvalhilolo_h = new double*[dim];
   double **funvallololo_h = new double*[dim];
   double **funvalhihihi_d = new double*[dim];
   double **funvallohihi_d = new double*[dim];
   double **funvalhilohi_d = new double*[dim];
   double **funvallolohi_d = new double*[dim];
   double **funvalhihilo_d = new double*[dim];
   double **funvallohilo_d = new double*[dim];
   double **funvalhilolo_d = new double*[dim];
   double **funvallololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhihihi_h[i] = new double[degp1];
      funvallohihi_h[i] = new double[degp1];
      funvalhilohi_h[i] = new double[degp1];
      funvallolohi_h[i] = new double[degp1];
      funvalhihilo_h[i] = new double[degp1];
      funvallohilo_h[i] = new double[degp1];
      funvalhilolo_h[i] = new double[degp1];
      funvallololo_h[i] = new double[degp1];
      funvalhihihi_d[i] = new double[degp1];
      funvallohihi_d[i] = new double[degp1];
      funvalhilohi_d[i] = new double[degp1];
      funvallolohi_d[i] = new double[degp1];
      funvalhihilo_d[i] = new double[degp1];
      funvallohilo_d[i] = new double[degp1];
      funvalhilolo_d[i] = new double[degp1];
      funvallololo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihihi_h = new double**[degp1];
   double ***jacvallohihi_h = new double**[degp1];
   double ***jacvalhilohi_h = new double**[degp1];
   double ***jacvallolohi_h = new double**[degp1];
   double ***jacvalhihilo_h = new double**[degp1];
   double ***jacvallohilo_h = new double**[degp1];
   double ***jacvalhilolo_h = new double**[degp1];
   double ***jacvallololo_h = new double**[degp1];
   double ***jacvalhihihi_d = new double**[degp1];
   double ***jacvallohihi_d = new double**[degp1];
   double ***jacvalhilohi_d = new double**[degp1];
   double ***jacvallolohi_d = new double**[degp1];
   double ***jacvalhihilo_d = new double**[degp1];
   double ***jacvallohilo_d = new double**[degp1];
   double ***jacvalhilolo_d = new double**[degp1];
   double ***jacvallololo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhihihi_h[i] = new double*[dim];
      jacvallohihi_h[i] = new double*[dim];
      jacvalhilohi_h[i] = new double*[dim];
      jacvallolohi_h[i] = new double*[dim];
      jacvalhihilo_h[i] = new double*[dim];
      jacvallohilo_h[i] = new double*[dim];
      jacvalhilolo_h[i] = new double*[dim];
      jacvallololo_h[i] = new double*[dim];
      jacvalhihihi_d[i] = new double*[dim];
      jacvallohihi_d[i] = new double*[dim];
      jacvalhilohi_d[i] = new double*[dim];
      jacvallolohi_d[i] = new double*[dim];
      jacvalhihilo_d[i] = new double*[dim];
      jacvallohilo_d[i] = new double*[dim];
      jacvalhilolo_d[i] = new double*[dim];
      jacvallololo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhihihi_h[i][j] = new double[dim];
         jacvallohihi_h[i][j] = new double[dim];
         jacvalhilohi_h[i][j] = new double[dim];
         jacvallolohi_h[i][j] = new double[dim];
         jacvalhihilo_h[i][j] = new double[dim];
         jacvallohilo_h[i][j] = new double[dim];
         jacvalhilolo_h[i][j] = new double[dim];
         jacvallololo_h[i][j] = new double[dim];
         jacvalhihihi_d[i][j] = new double[dim];
         jacvallohihi_d[i][j] = new double[dim];
         jacvalhilohi_d[i][j] = new double[dim];
         jacvallolohi_d[i][j] = new double[dim];
         jacvalhihilo_d[i][j] = new double[dim];
         jacvallohilo_d[i][j] = new double[dim];
         jacvalhilolo_d[i][j] = new double[dim];
         jacvallololo_d[i][j] = new double[dim];
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihihi_h = new double*[degp1];
   double **sollohihi_h = new double*[degp1];
   double **solhilohi_h = new double*[degp1];
   double **sollolohi_h = new double*[degp1];
   double **solhihilo_h = new double*[degp1];
   double **sollohilo_h = new double*[degp1];
   double **solhilolo_h = new double*[degp1];
   double **sollololo_h = new double*[degp1];
   double **solhihihi_d = new double*[degp1];
   double **sollohihi_d = new double*[degp1];
   double **solhilohi_d = new double*[degp1];
   double **sollolohi_d = new double*[degp1];
   double **solhihilo_d = new double*[degp1];
   double **sollohilo_d = new double*[degp1];
   double **solhilolo_d = new double*[degp1];
   double **sollololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihihi_h[i] = new double[dim];
      sollohihi_h[i] = new double[dim];
      solhilohi_h[i] = new double[dim];
      sollolohi_h[i] = new double[dim];
      solhihilo_h[i] = new double[dim];
      sollohilo_h[i] = new double[dim];
      solhilolo_h[i] = new double[dim];
      sollololo_h[i] = new double[dim];
      solhihihi_d[i] = new double[dim];
      sollohihi_d[i] = new double[dim];
      solhilohi_d[i] = new double[dim];
      sollolohi_d[i] = new double[dim];
      solhihilo_d[i] = new double[dim];
      sollohilo_d[i] = new double[dim];
      solhilolo_d[i] = new double[dim];
      sollololo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihihi_h = new double*[degp1];
   double **rhslohihi_h = new double*[degp1];
   double **rhshilohi_h = new double*[degp1];
   double **rhslolohi_h = new double*[degp1];
   double **rhshihilo_h = new double*[degp1];
   double **rhslohilo_h = new double*[degp1];
   double **rhshilolo_h = new double*[degp1];
   double **rhslololo_h = new double*[degp1];
   double **rhshihihi_d = new double*[degp1];
   double **rhslohihi_d = new double*[degp1];
   double **rhshilohi_d = new double*[degp1];
   double **rhslolohi_d = new double*[degp1];
   double **rhshihilo_d = new double*[degp1];
   double **rhslohilo_d = new double*[degp1];
   double **rhshilolo_d = new double*[degp1];
   double **rhslololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihihi_h[i] = new double[dim];
      rhslohihi_h[i] = new double[dim];
      rhshilohi_h[i] = new double[dim];
      rhslolohi_h[i] = new double[dim];
      rhshihilo_h[i] = new double[dim];
      rhslohilo_h[i] = new double[dim];
      rhshilolo_h[i] = new double[dim];
      rhslololo_h[i] = new double[dim];
      rhshihihi_d[i] = new double[dim];
      rhslohihi_d[i] = new double[dim];
      rhshilohi_d[i] = new double[dim];
      rhslolohi_d[i] = new double[dim];
      rhshihilo_d[i] = new double[dim];
      rhslohilo_d[i] = new double[dim];
      rhshilolo_d[i] = new double[dim];
      rhslololo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihihi_h = new double*[degp1];
   double **urhslohihi_h = new double*[degp1];
   double **urhshilohi_h = new double*[degp1];
   double **urhslolohi_h = new double*[degp1];
   double **urhshihilo_h = new double*[degp1];
   double **urhslohilo_h = new double*[degp1];
   double **urhshilolo_h = new double*[degp1];
   double **urhslololo_h = new double*[degp1];
   double **urhshihihi_d = new double*[degp1];
   double **urhslohihi_d = new double*[degp1];
   double **urhshilohi_d = new double*[degp1];
   double **urhslolohi_d = new double*[degp1];
   double **urhshihilo_d = new double*[degp1];
   double **urhslohilo_d = new double*[degp1];
   double **urhshilolo_d = new double*[degp1];
   double **urhslololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhshihihi_h[i] = new double[dim];
      urhslohihi_h[i] = new double[dim];
      urhshilohi_h[i] = new double[dim];
      urhslolohi_h[i] = new double[dim];
      urhshihilo_h[i] = new double[dim];
      urhslohilo_h[i] = new double[dim];
      urhshilolo_h[i] = new double[dim];
      urhslololo_h[i] = new double[dim];
      urhshihihi_d[i] = new double[dim];
      urhslohihi_d[i] = new double[dim];
      urhshilohi_d[i] = new double[dim];
      urhslolohi_d[i] = new double[dim];
      urhshihilo_d[i] = new double[dim];
      urhslohilo_d[i] = new double[dim];
      urhshilolo_d[i] = new double[dim];
      urhslololo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihihi = new double*[dim];
   double **workmatlohihi = new double*[dim];
   double **workmathilohi = new double*[dim];
   double **workmatlolohi = new double*[dim];
   double **workmathihilo = new double*[dim];
   double **workmatlohilo = new double*[dim];
   double **workmathilolo = new double*[dim];
   double **workmatlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihihi[i] = new double[dim];
      workmatlohihi[i] = new double[dim];
      workmathilohi[i] = new double[dim];
      workmatlolohi[i] = new double[dim];
      workmathihilo[i] = new double[dim];
      workmatlohilo[i] = new double[dim];
      workmathilolo[i] = new double[dim];
      workmatlololo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
   double *workvechihihi = new double[dim];
   double *workveclohihi = new double[dim];
   double *workvechilohi = new double[dim];
   double *workveclolohi = new double[dim];
   double *workvechihilo = new double[dim];
   double *workveclohilo = new double[dim];
   double *workvechilolo = new double[dim];
   double *workveclololo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshihihi = new double*[degp1];
   double **workrhslohihi = new double*[degp1];
   double **workrhshilohi = new double*[degp1];
   double **workrhslolohi = new double*[degp1];
   double **workrhshihilo = new double*[degp1];
   double **workrhslohilo = new double*[degp1];
   double **workrhshilolo = new double*[degp1];
   double **workrhslololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshihihi[i] = new double[dim];
      workrhslohihi[i] = new double[dim];
      workrhshilohi[i] = new double[dim];
      workrhslolohi[i] = new double[dim];
      workrhshihilo[i] = new double[dim];
      workrhslohilo[i] = new double[dim];
      workrhshilolo[i] = new double[dim];
      workrhslololo[i] = new double[dim];
   }
   double **resvechihihi = new double*[degp1];
   double **resveclohihi = new double*[degp1];
   double **resvechilohi = new double*[degp1];
   double **resveclolohi = new double*[degp1];
   double **resvechihilo = new double*[degp1];
   double **resveclohilo = new double*[degp1];
   double **resvechilolo = new double*[degp1];
   double **resveclololo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechihihi[i] = new double[dim];
      resveclohihi[i] = new double[dim];
      resvechilohi[i] = new double[dim];
      resveclolohi[i] = new double[dim];
      resvechihilo[i] = new double[dim];
      resveclohilo[i] = new double[dim];
      resvechilolo[i] = new double[dim];
      resveclololo[i] = new double[dim];
   }
   double resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi;
   double resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo;

   double **Qhihihi_h = new double*[dim];
   double **Qlohihi_h = new double*[dim];
   double **Qhilohi_h = new double*[dim];
   double **Qlolohi_h = new double*[dim];
   double **Qhihilo_h = new double*[dim];
   double **Qlohilo_h = new double*[dim];
   double **Qhilolo_h = new double*[dim];
   double **Qlololo_h = new double*[dim];
   double **Qhihihi_d = new double*[dim];
   double **Qlohihi_d = new double*[dim];
   double **Qhilohi_d = new double*[dim];
   double **Qlolohi_d = new double*[dim];
   double **Qhihilo_d = new double*[dim];
   double **Qlohilo_d = new double*[dim];
   double **Qhilolo_d = new double*[dim];
   double **Qlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qhihihi_h[i] = new double[dim];
      Qlohihi_h[i] = new double[dim];
      Qhilohi_h[i] = new double[dim];
      Qlolohi_h[i] = new double[dim];
      Qhihilo_h[i] = new double[dim];
      Qlohilo_h[i] = new double[dim];
      Qhilolo_h[i] = new double[dim];
      Qlololo_h[i] = new double[dim];
      Qhihihi_d[i] = new double[dim];
      Qlohihi_d[i] = new double[dim];
      Qhilohi_d[i] = new double[dim];
      Qlolohi_d[i] = new double[dim];
      Qhihilo_d[i] = new double[dim];
      Qlohilo_d[i] = new double[dim];
      Qhilolo_d[i] = new double[dim];
      Qlololo_d[i] = new double[dim];
   }
   double **Rhihihi_h = new double*[dim];
   double **Rlohihi_h = new double*[dim];
   double **Rhilohi_h = new double*[dim];
   double **Rlolohi_h = new double*[dim];
   double **Rhihilo_h = new double*[dim];
   double **Rlohilo_h = new double*[dim];
   double **Rhilolo_h = new double*[dim];
   double **Rlololo_h = new double*[dim];
   double **Rhihihi_d = new double*[dim];
   double **Rlohihi_d = new double*[dim];
   double **Rhilohi_d = new double*[dim];
   double **Rlolohi_d = new double*[dim];
   double **Rhihilo_d = new double*[dim];
   double **Rlohilo_d = new double*[dim];
   double **Rhilolo_d = new double*[dim];
   double **Rlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhihihi_h[i] = new double[dim];
      Rlohihi_h[i] = new double[dim];
      Rhilohi_h[i] = new double[dim];
      Rlolohi_h[i] = new double[dim];
      Rhihilo_h[i] = new double[dim];
      Rlohilo_h[i] = new double[dim];
      Rhilolo_h[i] = new double[dim];
      Rlololo_h[i] = new double[dim];
      Rhihihi_d[i] = new double[dim];
      Rlohihi_d[i] = new double[dim];
      Rhilohi_d[i] = new double[dim];
      Rlolohi_d[i] = new double[dim];
      Rhihilo_d[i] = new double[dim];
      Rlohilo_d[i] = new double[dim];
      Rhilolo_d[i] = new double[dim];
      Rlololo_d[i] = new double[dim];
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl8_unit_series_vector
      (dim,deg,inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
               inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h);
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihihi_d[i][j] = inputhihihi_h[i][j];
         inputlohihi_d[i][j] = inputlohihi_h[i][j];
         inputhilohi_d[i][j] = inputhilohi_h[i][j];
         inputlolohi_d[i][j] = inputlolohi_h[i][j];
         inputhihilo_d[i][j] = inputhihilo_h[i][j];
         inputlohilo_d[i][j] = inputlohilo_h[i][j];
         inputhilolo_d[i][j] = inputhilolo_h[i][j];
         inputlololo_d[i][j] = inputlololo_h[i][j];
      }

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihihi_h[i][0] << "  "
                            << inputlohihi_h[i][0] << endl;
         cout << "     " << inputhilohi_h[i][0] << "  "
                         << inputlolohi_h[i][0] << endl;
         cout << "     " << inputhihilo_h[i][0] << "  "
                         << inputlohilo_h[i][0] << endl;
         cout << "     " << inputhilolo_h[i][0] << "  "
                         << inputlololo_h[i][0] << endl;
      }
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;
/*
      dbl8_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi, inputlohihi, inputhilohi, inputlolohi,
          inputhihilo, inputlohilo, inputhilolo, inputlololo,
          outputhihihi,outputlohihi,outputhilohi,outputlolohi,
          outputhihilo,outputlohilo,outputhilolo,outputlololo,
          funvalhihihi,funvallohihi,funvalhilohi,funvallolohi,
          funvalhihilo,funvallohilo,funvalhilolo,funvallololo,
          jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
          jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          workrhshihihi,workrhslohihi,workrhshilohi,workrhslolohi,
          workrhshihilo,workrhslohilo,workrhshilolo,workrhslololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          ipvt,vrblvl);
 */
      dbl8_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
          outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
          outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
          outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
          funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
          funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
          funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
          funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
          rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
          rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
          rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
          rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
          urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
          urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
          urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
          urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          vrblvl,mode);
   }
   return 0;
}
