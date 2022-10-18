// The file cmplx4_newton_method.cpp defines the functions with prototypes in
// the file cmplx4_newton_method.h.

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
#include "dbl4_tail_kernels.h"
#include "dbl4_systems_host.h"
#include "dbl4_systems_kernels.h"
#include "dbl4_newton_testers.h"

using namespace std;

void cmplx4_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac, double dpr,
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

   if((mode == 1) || (mode == 2))
   {
      // The series coefficients accumulate common factors,
      // initially the coefficients are set to one.
      cmplx4_unit_series_vector
         (dim,deg,cffrehihi,cffrelohi,cffrehilo,cffrelolo,
                  cffimhihi,cffimlohi,cffimhilo,cffimlolo);

      if(vrblvl > 0)
         cout << "calling CPU_cmplx4_evaluate_monomials ..." << endl;

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
      // reset the coefficients
      cmplx4_unit_series_vector
         (dim,deg,cffrehihi,cffrelohi,cffrehilo,cffrelolo,
                  cffimhihi,cffimlohi,cffimhilo,cffimlolo);

      if(vrblvl > 0)
         cout << "calling GPU_cmplx4_evaluate_monomials ..." << endl;

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
      double errsum = 0.0;

      cout << "comparing CPU with GPU evaluations ... " << endl;

      errsum = cmplx4_error3sum(dim,dim+1,degp1,
                  outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
                  outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
                  outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
                  outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
                  "output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

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

   if(vrblvl > 0) cout << "linearizing the output ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      cmplx4_linearize_evaldiff_output
         (dim,degp1,nvr,idx,dpr,
          outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
          outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
          funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
          funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
          rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
          rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
          jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
          jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
          vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx4_linearize_evaldiff_output
         (dim,degp1,nvr,idx,dpr,
          outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
          outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
          funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
          funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
          rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
          rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
          jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
          jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
          vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << "comparing CPU with GPU function values ... " << endl;
      errsum = cmplx4_error2sum(dim,degp1,
                  funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
                  funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
                  funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
                  funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
                  "funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = cmplx4_error3sum(degp1,dim,dim,
                  jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
                  jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
                  jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
                  jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
                  "jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = cmplx4_error2sum(degp1,dim,
                  rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
                  rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
                  rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
                  rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
                  "rhs",vrblvl);
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
      if(vrblvl > 0)
         cout << "calling CPU_cmplx4_qrbs_solve ..." << endl;

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
         cout << "calling CPU_cmplx4_linear_residue ..." << endl;

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
      if(vrblvl > 0)
         cout << "calling GPU_cmplx4_bals_solve ..." << endl;

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
         cout << "calling GPU_cmplx4_linear_residue ..." << endl;

         double elapsedms;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_cmplx4_linear_residue
            (dim,degp1,szt,nbt,
             jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
             jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
             rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
             rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
             solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
             solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
             resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
             resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
             &elapsedms,&addcnt,&mulcnt,vrblvl);

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
      errsum = cmplx4_error2sum(dim,dim,
                  Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
                  Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
                  Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
                  Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = cmplx4_error2sum(dim,dim,
                  Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
                  Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
                  Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
                  Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = cmplx4_error2sum(degp1,dim,
                  urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
                  urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
                  urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
                  urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
                  "urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = cmplx4_error2sum(degp1,dim,
                  solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
                  solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
                  solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
                  solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
                  "sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = cmplx4_error2sum(dim,degp1,
                  inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
                  inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
                  inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
                  inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
                  "input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl4_complex_newton
 ( int szt, int nbt, int dim, int deg,
   int *nvr, int **idx, int **exp, int *nbrfac, int **expfac,
   double dpr, int nbsteps, int mode, int vrblvl )
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
   cmplx4_start_series_vector
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
   if(vrblvl > 0) cout << scientific << setprecision(16);

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step << " ***" << endl;

      cmplx4_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,dpr,
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
   if(vrblvl < 2)
   {
      cout << scientific << setprecision(16); // just in case vrblvl == 0
      cout << "The solution series : " << endl;
      for(int j=0; j<degp1; j++)
      {
         cout << "coefficient of degree " << j << " :" << endl;
         for(int i=0; i<dim; i++)
         {
            if((mode == 0) || (mode == 2))
              cout << "x_d[" << i << "][" << j << "] : "
                             << inputrehihi_d[i][j] << "  "
                             << inputrelohi_d[i][j] << endl << "  "
                             << inputrehilo_d[i][j] << "  "
                             << inputrelolo_d[i][j] << endl << "  "
                             << inputimhihi_d[i][j] << "  "
                             << inputimlohi_d[i][j] << endl << "  "
                             << inputimhilo_d[i][j] << "  "
                             << inputimlolo_d[i][j] << endl;
            if((mode == 1) || (mode == 2))
              cout << "x_h[" << i << "][" << j << "] : "
                             << inputrehihi_h[i][j] << "  "
                             << inputrelohi_h[i][j] << endl << "  "
                             << inputrehilo_h[i][j] << "  "
                             << inputrelolo_h[i][j] << endl << "  "
                             << inputimhihi_h[i][j] << "  "
                             << inputimlohi_h[i][j] << endl << "  "
                             << inputimhilo_h[i][j] << "  "
                             << inputimlolo_h[i][j] << endl;
         }
      }
   }
   return 0;
}
