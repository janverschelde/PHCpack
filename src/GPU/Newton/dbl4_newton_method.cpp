// The file dbl4_newton_method.cpp defines the functions with prototypes in
// the file dbl4_newton_method.h.

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <vector_types.h>
#include <time.h>
#ifdef winwalltime
#include "gettimeofday4win.h"
#else
#include <sys/time.h>
#endif
#include "unimodular_matrices.h"
#include "random_numbers.h"
#include "random4_vectors.h"
#include "quad_double_functions.h"
#include "dbl4_indexed_coefficients.h"
#include "dbl4_convolutions_host.h"
#include "dbl4_monomials_host.h"
#include "dbl4_polynomials_host.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"
#include "job_makers.h"
#include "dbl4_polynomials_kernels.h"
#include "dbl4_factorizations.h"
#include "dbl4_monomial_systems.h"
#include "dbl4_bals_host.h"
#include "dbl4_bals_kernels.h"
#include "dbl4_tail_kernels.h"
#include "dbl4_systems_host.h"
#include "dbl4_systems_kernels.h"
#include "dbl4_newton_testers.h"
#include "write_newton_times.h"

using namespace std;

int dbl4_errors_funjacrhs
 ( int dim, int deg,
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
   double **rhshilo_d, double **rhslolo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-50;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = dbl4_error2sum(dim,degp1,
                funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
                funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
                "funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = dbl4_error3sum(degp1,dim,dim,
                jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
                jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
                "jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   cout << scientific << setprecision(16);
   errsum = dbl4_error2sum(degp1,dim,
                rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
                rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,"rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl4_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-50;
   double errsum = 0.0;
   int fail = 0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   errsum = dbl4_error2sum(dim,dim,
                Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,
                Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = dbl4_error2sum(dim,dim,
                Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
                Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = dbl4_error2sum(degp1,dim,
                urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
                urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,"urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = dbl4_error2sum(degp1,dim,
                solhihi_h,sollohi_h,solhilo_h,sollolo_h,
                solhihi_d,sollohi_d,solhilo_d,sollolo_d,"sol",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = dbl4_error2sum(dim,degp1,
                inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
                inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
                "input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl4_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
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
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhshihi_h[i][j] = rhshihi_h[i][j];
            urhslohi_h[i][j] = rhslohi_h[i][j];
            urhshilo_h[i][j] = rhshilo_h[i][j];
            urhslolo_h[i][j] = rhslolo_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhihi_h[i][j] = 0.0; sollohi_h[i][j] = 0.0;
            solhilo_h[i][j] = 0.0; sollolo_h[i][j] = 0.0;
         }
 
      if(vrblvl > 0) cout << "calling CPU_dbl4_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_dbl4_qrbs_solve
         (dim,degp1,oldtail,
          jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
          urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
          solhihi_h,sollohi_h,solhilo_h,sollolo_h,
          Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
          workvechihi,workveclohi,workvechilo,workveclolo,
          zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl4_linear_residue ..." << endl;

         CPU_dbl4_linear_residue
            (dim,degp1,*tailidx_h-1,
             jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
             rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
             solhihi_h,sollohi_h,solhilo_h,sollolo_h,
             resvechihi,resveclohi,resvechilo,resveclolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhihi << endl;
      }
      dbl4_update_series
         (dim,degp1,*tailidx_h-1,
          inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
          solhihi_h,sollohi_h,solhilo_h,sollolo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhshihi_d[i][j] = rhshihi_d[i][j];
            urhslohi_d[i][j] = rhslohi_d[i][j];
            urhshilo_d[i][j] = rhshilo_d[i][j];
            urhslolo_d[i][j] = rhslolo_d[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhihi_d[i][j] = 0.0; sollohi_d[i][j] = 0.0;
            solhilo_d[i][j] = 0.0; sollolo_d[i][j] = 0.0;
         }
 
      if(vrblvl > 0) cout << "calling GPU_dbl4_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_dbl4_bals_solve
         (dim,degp1,szt,nbt,oldtail,
          jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
          Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
          urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
          solhihi_d,sollohi_d,solhilo_d,sollolo_d,zeroQ_d,noqr_d,
          upidx_d,bsidx_d,&newtail,totqrlapsedms,totqtblapsedms,
          totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl4_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_dbl4_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,
             jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
             rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
             solhihi_d,sollohi_d,solhilo_d,sollolo_d,
             resvechihi,resveclohi,resvechilo,resveclolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
             &elapsedms,&addcnt,&mulcnt,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhihi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;

         *totreslapsedms += elapsedms;
      }
      dbl4_update_series
         (dim,degp1,*tailidx_d-1,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
          solhihi_d,sollohi_d,solhilo_d,sollolo_d,vrblvl);
   }
   return 0;
}

int dbl4_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbhihi, double **mblohi, double **mbhilo, double **mblolo,
   double dpr,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double **acchihi, double **acclohi, double **acchilo, double **acclolo,
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
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms, 
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_dbl4_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl4_unit_series_vector
            (dim,deg,cffhihi[0],cfflohi[0],cffhilo[0],cfflolo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_dbl4_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhihi[0],cfflohi[0],cffhilo[0],cfflolo[0],
             acchihi[0],acclohi[0],acchilo[0],acclolo[0],
             inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
             outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,vrblvl);
      }
      else
         CPU_dbl4_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,
             cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
             inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
             funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
             jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl4_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl4_unit_series_vector
            (dim,deg,cffhihi[0],cfflohi[0],cffhilo[0],cfflolo[0]);
         // reset the coefficients

         GPU_dbl4_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhihi[0],cfflohi[0],cffhilo[0],cfflolo[0],
             acchihi[0],acclohi[0],acchilo[0],acclolo[0],
             inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
             outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
             totcnvlapsedms,vrblvl);
      }
      else
         GPU_dbl4_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,
             cffhihi,cfflohi,cffhilo,cfflolo,
             inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
             outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
             funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
             jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
             totcnvlapsedms,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      double errsum = 0.0;

      cout << "comparing CPU with GPU evaluations ... " << endl;

      errsum = dbl4_error3sum(dim,dim+1,degp1,
                  outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
                  outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
                  "output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl4_define_rhs
            (dim,degp1,mbhihi,mblohi,mbhilo,mblolo,
             funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
             rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalhihi_h[i][j][k] = 0.0; jacvallohi_h[i][j][k] = 0.0;
                  jacvalhilo_h[i][j][k] = 0.0; jacvallolo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl4_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbhihi,mblohi,mbhilo,mblolo,dpr,
             outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
             funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
             rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
             jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl4_define_rhs
            (dim,degp1,mbhihi,mblohi,mbhilo,mblolo,
             funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
             rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalhihi_d[i][j][k] = 0.0; jacvallohi_d[i][j][k] = 0.0;
                  jacvalhilo_d[i][j][k] = 0.0; jacvallolo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl4_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbhihi,mblohi,mbhilo,mblolo,dpr,
             outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
             funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
             rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
             jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl4_errors_funjacrhs(dim,deg,
                     funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
                     funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
                     jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
                     jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
                     rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
                     rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,vrblvl);
   }
   dbl4_update_newton_qrstep
       (szt,nbt,dim,deg,tailidx_h,tailidx_d,
        inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
        inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
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
        workvechihi,workveclohi,workvechilo,workveclolo,
        resvechihi,resveclohi,resvechilo,resveclolo,
        resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
        zeroQ_h,noqr_h,zeroQ_d,noqr_d,
        upidx_h,bsidx_h,upidx_d,bsidx_d,
        totqrlapsedms,totqtblapsedms,totbslapsedms,
        totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return dbl4_errors_inurhsQRsol(dim,deg,
                 inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
                 inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
                 Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,
                 Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
                 Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
                 Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
                 urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
                 urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
                 solhihi_h,sollohi_h,solhilo_h,sollolo_h,
                 solhihi_d,sollohi_d,solhilo_d,sollolo_d,vrblvl);
   }
   else
      return 0;
}

int dbl4_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **csthihi, double **cstlohi, double **csthilo, double **cstlolo,
   double ***cffhihi, double ***cfflohi, double ***cffhilo, double ***cfflolo,
   double dpr,
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
   double *workvechihi, double *workveclohi,
   double *workvechilo, double *workveclolo,
   double **resvechihi, double **resveclohi, 
   double **resvechilo, double **resveclolo, 
   double *resmaxhihi, double *resmaxlohi,
   double *resmaxhilo, double *resmaxlolo,
   bool *zeroQ_h, bool *noqr_h, bool *zeroQ_d, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_dbl4_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl4_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
            csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
            cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
            inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
            outputhihi_h[i],outputlohi_h[i],
            outputhilo_h[i],outputlolo_h[i],&lapsed,0);

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << lapsed << " seconds." << endl;

         timelapsed_h += lapsed;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_h << " seconds." << endl;
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl4_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ConvolutionJobs cnvjobs(dim);
         AdditionJobs addjobs(dim,nbr[i]);

         make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,false);

         GPU_dbl4_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
             cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
             inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
             outputhihi_d[i],outputlohi_d[i],outputhilo_d[i],outputlolo_d[i],
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0);

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << walltimes_d << " seconds." << endl;

         timelapsed_d += walltimes_d;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_d << " seconds." << endl;
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << scientific << setprecision(16)
           << "comparing CPU with GPU evaluations ... " << endl;
   
      double errsum = 0.0;

      errsum = dbl4_error3sum(dim,dim+1,degp1,
                  outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
                  outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
                  "output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(2);
      cout << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0)
      cout << "mapping the output to values of function and matrix series ..."
           << endl;

   if((mode == 1) || (mode == 2))
   {
      dbl4_map_evaldiff_output(dim,deg,
         outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
         funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
         jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshihi_h[i][j] = -funvalhihi_h[j][i];
            rhslohi_h[i][j] = -funvallohi_h[j][i];
            rhshilo_h[i][j] = -funvalhilo_h[j][i];
            rhslolo_h[i][j] = -funvallolo_h[j][i];
         }
   }
   if((mode == 0) || (mode == 2))
   {
      dbl4_map_evaldiff_output(dim,deg,
         outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
         funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
         jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshihi_d[i][j] = -funvalhihi_d[j][i];
            rhslohi_d[i][j] = -funvallohi_d[j][i];
            rhshilo_d[i][j] = -funvalhilo_d[j][i];
            rhslolo_d[i][j] = -funvallolo_d[j][i];
         }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl4_errors_funjacrhs(dim,deg,
                     funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
                     funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
                     jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
                     jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
                     rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
                     rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,vrblvl);
   }
   dbl4_update_newton_qrstep
       (szt,nbt,dim,deg,tailidx_h,tailidx_d,
        inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
        inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
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
        workvechihi,workveclohi,workvechilo,workveclolo,
        resvechihi,resveclohi,resvechilo,resveclolo,
        resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
        zeroQ_h,noqr_h,zeroQ_d,noqr_d,
        upidx_h,bsidx_h,upidx_d,bsidx_d,
        totqrlapsedms,totqtblapsedms,totbslapsedms,
        totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return dbl4_errors_inurhsQRsol(dim,deg,
                 inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
                 inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
                 Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,
                 Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
                 Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,
                 Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
                 urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
                 urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
                 solhihi_h,sollohi_h,solhilo_h,sollolo_h,
                 solhihi_d,sollohi_d,solhilo_d,sollolo_d,vrblvl);
   }
   else
      return 0;
}

int dbl4_allocate_inoutfunjac
 ( int dim, int deg, int mode,
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
   double ***jacvalhilo_d, double ***jacvallolo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputhihi_h[i] = new double[degp1];
         inputlohi_h[i] = new double[degp1];
         inputhilo_h[i] = new double[degp1];
         inputlolo_h[i] = new double[degp1];

         outputhihi_h[i] = new double*[dim+1];
         outputlohi_h[i] = new double*[dim+1];
         outputhilo_h[i] = new double*[dim+1];
         outputlolo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhihi_h[i][j] = new double[degp1];
            outputlohi_h[i][j] = new double[degp1];
            outputhilo_h[i][j] = new double[degp1];
            outputlolo_h[i][j] = new double[degp1];
         }
         funvalhihi_h[i] = new double[degp1];
         funvallohi_h[i] = new double[degp1];
         funvalhilo_h[i] = new double[degp1];
         funvallolo_h[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhihi_h[i] = new double*[dim];
         jacvallohi_h[i] = new double*[dim];
         jacvalhilo_h[i] = new double*[dim];
         jacvallolo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhihi_h[i][j] = new double[dim];
            jacvallohi_h[i][j] = new double[dim];
            jacvalhilo_h[i][j] = new double[dim];
            jacvallolo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputhihi_d[i] = new double[degp1];
         inputlohi_d[i] = new double[degp1];
         inputhilo_d[i] = new double[degp1];
         inputlolo_d[i] = new double[degp1];

         outputhihi_d[i] = new double*[dim+1];
         outputlohi_d[i] = new double*[dim+1];
         outputhilo_d[i] = new double*[dim+1];
         outputlolo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhihi_d[i][j] = new double[degp1];
            outputlohi_d[i][j] = new double[degp1];
            outputhilo_d[i][j] = new double[degp1];
            outputlolo_d[i][j] = new double[degp1];
         }
         funvalhihi_d[i] = new double[degp1];
         funvallohi_d[i] = new double[degp1];
         funvalhilo_d[i] = new double[degp1];
         funvallolo_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhihi_d[i] = new double*[dim];
         jacvallohi_d[i] = new double*[dim];
         jacvalhilo_d[i] = new double*[dim];
         jacvallolo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhihi_d[i][j] = new double[dim];
            jacvallohi_d[i][j] = new double[dim];
            jacvalhilo_d[i][j] = new double[dim];
            jacvallolo_d[i][j] = new double[dim];
         }
      }
   }
   return 0;
}

int dbl4_allocate_rhsqrsol
 ( int dim, int deg, int mode,
   double **rhshihi_h, double **rhslohi_h,
   double **rhshilo_h, double **rhslolo_h,
   double **rhshihi_d, double **rhslohi_d,
   double **rhshilo_d, double **rhslolo_d,
   double **urhshihi_h, double **urhslohi_h,
   double **urhshilo_h, double **urhslolo_h,
   double **urhshihi_d, double **urhslohi_d,
   double **urhshilo_d, double **urhslolo_d,
   double **Qhihi_h, double **Qlohi_h, double **Qhilo_h, double **Qlolo_h,
   double **Qhihi_d, double **Qlohi_d, double **Qhilo_d, double **Qlolo_d,
   double **Rhihi_h, double **Rlohi_h, double **Rhilo_h, double **Rlolo_h,
   double **Rhihi_d, double **Rlohi_d, double **Rhilo_d, double **Rlolo_d,
   double **solhihi_h, double **sollohi_h,
   double **solhilo_h, double **sollolo_h,
   double **solhihi_d, double **sollohi_d,
   double **solhilo_d, double **sollolo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++)
      {
         rhshihi_h[i] = new double[dim];
         rhslohi_h[i] = new double[dim];
         rhshilo_h[i] = new double[dim];
         rhslolo_h[i] = new double[dim];
         urhshihi_h[i] = new double[dim];
         urhslohi_h[i] = new double[dim];
         urhshilo_h[i] = new double[dim];
         urhslolo_h[i] = new double[dim];
         solhihi_h[i] = new double[dim];
         sollohi_h[i] = new double[dim];
         solhilo_h[i] = new double[dim];
         sollolo_h[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qhihi_h[i] = new double[dim];
         Qlohi_h[i] = new double[dim];
         Qhilo_h[i] = new double[dim];
         Qlolo_h[i] = new double[dim];
         Rhihi_h[i] = new double[dim];
         Rlohi_h[i] = new double[dim];
         Rhilo_h[i] = new double[dim];
         Rlolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++)
      {
         rhshihi_d[i] = new double[dim];
         rhslohi_d[i] = new double[dim];
         rhshilo_d[i] = new double[dim];
         rhslolo_d[i] = new double[dim];
         urhshihi_d[i] = new double[dim];
         urhslohi_d[i] = new double[dim];
         urhshilo_d[i] = new double[dim];
         urhslolo_d[i] = new double[dim];
         solhihi_d[i] = new double[dim];
         sollohi_d[i] = new double[dim];
         solhilo_d[i] = new double[dim];
         sollolo_d[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qhihi_d[i] = new double[dim];
         Qlohi_d[i] = new double[dim];
         Qhilo_d[i] = new double[dim];
         Qlolo_d[i] = new double[dim];
         Rhihi_d[i] = new double[dim];
         Rlohi_d[i] = new double[dim];
         Rhilo_d[i] = new double[dim];
         Rlolo_d[i] = new double[dim];
      }
   }
   return 0;
}

void dbl4_start_setup
 ( int dim, int deg,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d, int mode, int vrblvl )
{
   double *start0hihi = new double[dim];
   double *start0lohi = new double[dim];
   double *start0hilo = new double[dim];
   double *start0lolo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0hihi[i] = testsolhihi[i][0];
      start0lohi[i] = testsollohi[i][0];
      start0hilo[i] = testsolhilo[i][0];
      start0lolo[i] = testsollolo[i][0];
   }
   if((mode == 1) || (mode == 2))
      real4_start_series_vector
         (dim,deg,start0hihi,start0lohi,start0hilo,start0lolo,
          inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h);
   else
      real4_start_series_vector
         (dim,deg,start0hihi,start0lohi,start0hilo,start0lolo,
          inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
         {
            inputhihi_d[i][j] = inputhihi_h[i][j];
            inputlohi_d[i][j] = inputlohi_h[i][j];
            inputhilo_d[i][j] = inputhilo_h[i][j];
            inputlolo_d[i][j] = inputlolo_h[i][j];
         }
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;

      if((mode == 1) || (mode == 2))
      {
         for(int i=0; i<dim; i++)
         {
            cout << i << " : " << inputhihi_h[i][0] << "  "
                               << inputlohi_h[i][0] << endl;
            cout << "     " << inputhilo_h[i][0] << "  "
                            << inputlolo_h[i][0] << endl;
         }
      }
      else
      {
         for(int i=0; i<dim; i++)
         {
            cout << i << " : " << inputhihi_d[i][0] << "  "
                               << inputlohi_d[i][0] << endl;
            cout << "     " << inputhilo_d[i][0] << "  "
                            << inputlolo_d[i][0] << endl;
         }
      }
   }
}

void dbl4_column_setup
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **mbrhshihi, double **mbrhslohi,
   double **mbrhshilo, double **mbrhslolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolhihi[i] = new double[degp1];
      testsollohi[i] = new double[degp1];
      testsolhilo[i] = new double[degp1];
      testsollolo[i] = new double[degp1];
   }
   make_real4_exponentials(dim,deg,
      testsolhihi,testsollohi,testsolhilo,testsollolo);

   // compute the right hand sides via evaluation

   for(int i=0; i<dim; i++)
   {
      mbrhshihi[i] = new double[degp1];
      mbrhslohi[i] = new double[degp1];
      mbrhshilo[i] = new double[degp1];
      mbrhslolo[i] = new double[degp1];

      mbrhshihi[i][0] = 1.0;     // initialize product to one
      mbrhslohi[i][0] = 0.0;
      mbrhshilo[i][0] = 0.0;
      mbrhslolo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhshihi[i][k] = 0.0;
         mbrhslohi[i][k] = 0.0;
         mbrhshilo[i][k] = 0.0;
         mbrhslolo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_real4_monomials
         (dim,deg,rowsA,testsolhihi,testsollohi,testsolhilo,testsollolo,
          mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo);
   else
      evaluate_real4_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffhihi,cfflohi,cffhilo,cfflolo,
          testsolhihi,testsollohi,testsolhilo,testsollolo,
          mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);

      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshihi[i][j] << "  " << mbrhslohi[i][j] << endl
                 << "  "
                 << mbrhshilo[i][j] << "  " << mbrhslolo[i][j] << endl;
   }
   dbl4_start_setup(dim,deg,
      testsolhihi,testsollohi,testsolhilo,testsollolo,
      inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
      inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,mode,vrblvl);
}

void dbl4_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihi, double **cstlohi,
   double **csthilo, double **cstlolo,
   double ***cffhihi, double ***cfflohi,
   double ***cffhilo, double ***cfflolo,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d,
   double ***outputhihi_h, double ***outputlohi_h,
   double ***outputhilo_h, double ***outputlolo_h,
   double ***outputhihi_d, double ***outputlohi_d,
   double ***outputhilo_d, double ***outputlolo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++) 
   {
      testsolhihi[i] = new double[degp1];
      testsollohi[i] = new double[degp1];
      testsolhilo[i] = new double[degp1];
      testsollolo[i] = new double[degp1];
   }
   make_real4_exponentials(dim,deg,
      testsolhihi,testsollohi,testsolhilo,testsollolo);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl4_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
            csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
            cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
            testsolhihi,testsollohi,testsolhilo,testsollolo,
            outputhihi_h[i],outputlohi_h[i],
            outputhilo_h[i],outputlolo_h[i],&lapsed,0);

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << lapsed << " seconds." << endl;

         timelapsed_h += lapsed;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_h << " seconds." << endl;

      for(int i=0; i<dim; i++) // adjust constant coefficients
      {
         for(int j=0; j<=deg; j++) // cst[i][j] -= output_h[i][dim][j];
            qdf_dec(&csthihi[i][j],&cstlohi[i][j],
                    &csthilo[i][j],&cstlolo[i][j],
                    outputhihi_h[i][dim][j],outputlohi_h[i][dim][j],
                    outputhilo_h[i][dim][j],outputlolo_h[i][dim][j]);
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_dbl4_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
               csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
               cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
               testsolhihi,testsollohi,testsolhilo,testsollolo,
               outputhihi_h[i],outputlohi_h[i],
               outputhilo_h[i],outputlolo_h[i],&lapsed,0);

            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << lapsed << " seconds." << endl;

            timelapsed_h += lapsed;
         }
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_h << " seconds." << endl;

         double errsum = 0.0;

         for(int i=0; i<dim; i++)
            for(int j=0; j<=deg; j++)
               errsum = errsum + outputhihi_h[i][dim][j]
                               + outputlohi_h[i][dim][j]
                               + outputhilo_h[i][dim][j]
                               + outputlolo_h[i][dim][j];

         cout << scientific << setprecision(2)
              << "residual of test solution : " << errsum << endl;
      }
   }
   else // GPU is faster
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the device ..." << endl;

      const bool vrb = false; // no output (vrblvl > 1);
      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ConvolutionJobs cnvjobs(dim);
         AdditionJobs addjobs(dim,nbr[i]);

         make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrb);

         GPU_dbl4_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
             cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
             testsolhihi,testsollohi,testsolhilo,testsollolo,
             outputhihi_d[i],outputlohi_d[i],outputhilo_d[i],outputlolo_d[i],
             cnvjobs,addjobs,&cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0);

         if(vrblvl > 0)
            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << walltimes_d << " seconds." << endl;

         timelapsed_d += walltimes_d;
      }
      if(vrblvl > 0)
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_d << " seconds." << endl;

      for(int i=0; i<dim; i++) // adjust constant coefficients
      {
         for(int j=0; j<=deg; j++) // cst[i][j] -= output_d[i][dim][j];
            qdf_dec(&csthihi[i][j],&cstlohi[i][j],
                    &csthilo[i][j],&cstlolo[i][j],
                    outputhihi_d[i][dim][j],outputlohi_d[i][dim][j],
                    outputhilo_d[i][dim][j],outputlolo_d[i][dim][j]);
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double cnvlapms,addlapms,timelapms_d,walltimes_d;
   
            ConvolutionJobs cnvjobs(dim);
            AdditionJobs addjobs(dim,nbr[i]);

            make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,vrb);

            GPU_dbl4_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],
                csthihi[i],cstlohi[i],csthilo[i],cstlolo[i],
                cffhihi[i],cfflohi[i],cffhilo[i],cfflolo[i],
                testsolhihi,testsollohi,testsolhilo,testsollolo,
                outputhihi_d[i],outputlohi_d[i],
                outputhilo_d[i],outputlolo_d[i],cnvjobs,addjobs,
                &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0);

            cout << fixed << setprecision(3)
                 << "Evaluated and differentiated polynomial " << i
                 << " in " << walltimes_d << " seconds." << endl;

            timelapsed_d += walltimes_d;
         }
         cout << fixed << setprecision(3)
              << "Evaluated and differentiated system in "
              << timelapsed_d << " seconds." << endl;

         double errsum = 0.0;

         for(int i=0; i<dim; i++)
            for(int j=0; j<=deg; j++)
               errsum = errsum + outputhihi_d[i][dim][j]
                               + outputlohi_d[i][dim][j]
                               + outputhilo_d[i][dim][j]
                               + outputlolo_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "residual of test solution : " << errsum << endl;
      }
   }
   dbl4_start_setup(dim,deg,
       testsolhihi,testsollohi,testsolhilo,testsollolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,mode,vrblvl);
}

int dbl4_error_testsol
 ( int dim, int deg, int mode,
   double **testsolhihi, double **testsollohi,
   double **testsolhilo, double **testsollolo,
   double **inputhihi_h, double **inputlohi_h,
   double **inputhilo_h, double **inputlolo_h,
   double **inputhihi_d, double **inputlohi_d,
   double **inputhilo_d, double **inputlolo_d )
{
   double errsum = 0.0;

   cout << scientific << setprecision(16); // just in case vrblvl == 0
   cout << "The solution series : " << endl;

   for(int j=0; j<=deg; j++)
   {
      cout << "coefficient of degree " << j << " :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "sol[" << i << "][" << j << "] : "
                        << testsolhihi[i][j] << "  "
                        << testsollohi[i][j] << endl << "  "
                        << testsolhilo[i][j] << "  "
                        << testsollolo[i][j] << endl;

         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhihi_d[i][j] << "  "
                           << inputlohi_d[i][j] << endl << "  "
                           << inputhilo_d[i][j] << "  "
                           << inputlolo_d[i][j] << endl;

            errsum += abs(testsolhihi[i][j] - inputhihi_d[i][j])
                    + abs(testsollohi[i][j] - inputlohi_d[i][j])
                    + abs(testsolhilo[i][j] - inputhilo_d[i][j])
                    + abs(testsollolo[i][j] - inputlolo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhihi_h[i][j] << "  "
                           << inputlohi_h[i][j] << endl << "  "
                           << inputhilo_h[i][j] << "  "
                           << inputlolo_h[i][j] << endl;

            errsum += abs(testsolhihi[i][j] - inputhihi_h[i][j])
                    + abs(testsollohi[i][j] - inputlohi_h[i][j])
                    + abs(testsolhilo[i][j] - inputhilo_h[i][j])
                    + abs(testsollolo[i][j] - inputlolo_h[i][j]);
         }
      }
   }
   cout << scientific << setprecision(2)
        << "error : " << errsum << endl;

   return (errsum > 1.0e-50);
}

int test_dbl4_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **acchihi = new double*[dim+1]; // accumulated power series
   double **acclohi = new double*[dim+1]; // in one column
   double **acchilo = new double*[dim+1];
   double **acclolo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchihi[i] = new double[degp1];
      acclohi[i] = new double[degp1];
      acchilo[i] = new double[degp1];
      acclolo[i] = new double[degp1];
   }
   double ***cffhihi = new double**[nbrcol]; // coefficients of monomials
   double ***cfflohi = new double**[nbrcol];
   double ***cffhilo = new double**[nbrcol];
   double ***cfflolo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffhihi[i] = new double*[dim];
      cfflohi[i] = new double*[dim];
      cffhilo[i] = new double*[dim];
      cfflolo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffhihi[i][j] = new double[degp1];
         cfflohi[i][j] = new double[degp1];
         cffhilo[i][j] = new double[degp1];
         cfflolo[i][j] = new double[degp1];
      }
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_real4_coefficients(nbrcol,dim,cffhihi,cfflohi,cffhilo,cfflolo);

   double **inputhihi_h;
   double **inputlohi_h;
   double **inputhilo_h;
   double **inputlolo_h;
   double **inputhihi_d;
   double **inputlohi_d;
   double **inputhilo_d;
   double **inputlolo_d;
   double ***outputhihi_h;
   double ***outputlohi_h;
   double ***outputhilo_h;
   double ***outputlolo_h;
   double ***outputhihi_d;
   double ***outputlohi_d;
   double ***outputhilo_d;
   double ***outputlolo_d;
   double **funvalhihi_h;
   double **funvallohi_h;
   double **funvalhilo_h;
   double **funvallolo_h;
   double **funvalhihi_d;
   double **funvallohi_d;
   double **funvalhilo_d;
   double **funvallolo_d;
   double ***jacvalhihi_h;
   double ***jacvallohi_h;
   double ***jacvalhilo_h;
   double ***jacvallolo_h;
   double ***jacvalhihi_d;
   double ***jacvallohi_d;
   double ***jacvalhilo_d;
   double ***jacvallolo_d;

   if((mode == 1) || (mode == 2))
   {
      inputhihi_h = new double*[dim];
      inputlohi_h = new double*[dim];
      inputhilo_h = new double*[dim];
      inputlolo_h = new double*[dim];
      outputhihi_h = new double**[dim];
      outputlohi_h = new double**[dim];
      outputhilo_h = new double**[dim];
      outputlolo_h = new double**[dim];
      funvalhihi_h = new double*[dim];
      funvallohi_h = new double*[dim];
      funvalhilo_h = new double*[dim];
      funvallolo_h = new double*[dim];
      jacvalhihi_h = new double**[degp1];
      jacvallohi_h = new double**[degp1];
      jacvalhilo_h = new double**[degp1];
      jacvallolo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputhihi_d = new double*[dim];
      inputlohi_d = new double*[dim];
      inputhilo_d = new double*[dim];
      inputlolo_d = new double*[dim];
      outputhihi_d = new double**[dim];
      outputlohi_d = new double**[dim];
      outputhilo_d = new double**[dim];
      outputlolo_d = new double**[dim];
      funvalhihi_d = new double*[dim];
      funvallohi_d = new double*[dim];
      funvalhilo_d = new double*[dim];
      funvallolo_d = new double*[dim];
      jacvalhihi_d = new double**[degp1];
      jacvallohi_d = new double**[degp1];
      jacvalhilo_d = new double**[degp1];
      jacvallolo_d = new double**[degp1];
   }
   dbl4_allocate_inoutfunjac(dim,deg,mode,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
       outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
       outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
       funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
       funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
       jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
       jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhshihi_h;
   double **rhslohi_h;
   double **rhshilo_h;
   double **rhslolo_h;
   double **rhshihi_d;
   double **rhslohi_d;
   double **rhshilo_d;
   double **rhslolo_d;
   double **urhshihi_h;
   double **urhslohi_h;
   double **urhshilo_h;
   double **urhslolo_h;
   double **urhshihi_d;
   double **urhslohi_d;
   double **urhshilo_d;
   double **urhslolo_d;
   double **solhihi_h;
   double **sollohi_h;
   double **solhilo_h;
   double **sollolo_h;
   double **solhihi_d;
   double **sollohi_d;
   double **solhilo_d;
   double **sollolo_d;
   double **Qhihi_h;
   double **Qlohi_h;
   double **Qhilo_h;
   double **Qlolo_h;
   double **Qhihi_d;
   double **Qlohi_d;
   double **Qhilo_d;
   double **Qlolo_d;
   double **Rhihi_h;
   double **Rlohi_h;
   double **Rhilo_h;
   double **Rlolo_h;
   double **Rhihi_d;
   double **Rlohi_d;
   double **Rhilo_d;
   double **Rlolo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshihi_h = new double*[degp1];
      rhslohi_h = new double*[degp1];
      rhshilo_h = new double*[degp1];
      rhslolo_h = new double*[degp1];
      urhshihi_h = new double*[degp1];
      urhslohi_h = new double*[degp1];
      urhshilo_h = new double*[degp1];
      urhslolo_h = new double*[degp1];
      solhihi_h = new double*[degp1];
      sollohi_h = new double*[degp1];
      solhilo_h = new double*[degp1];
      sollolo_h = new double*[degp1];
      Qhihi_h = new double*[dim];
      Qlohi_h = new double*[dim];
      Qhilo_h = new double*[dim];
      Qlolo_h = new double*[dim];
      Rhihi_h = new double*[dim];
      Rlohi_h = new double*[dim];
      Rhilo_h = new double*[dim];
      Rlolo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhshihi_d = new double*[degp1];
      rhslohi_d = new double*[degp1];
      rhshilo_d = new double*[degp1];
      rhslolo_d = new double*[degp1];
      urhshihi_d = new double*[degp1];
      urhslohi_d = new double*[degp1];
      urhshilo_d = new double*[degp1];
      urhslolo_d = new double*[degp1];
      solhihi_d = new double*[degp1];
      sollohi_d = new double*[degp1];
      solhilo_d = new double*[degp1];
      sollolo_d = new double*[degp1];
      Qhihi_d = new double*[dim];
      Qlohi_d = new double*[dim];
      Qhilo_d = new double*[dim];
      Qlolo_d = new double*[dim];
      Rhihi_d = new double*[dim];
      Rlohi_d = new double*[dim];
      Rhilo_d = new double*[dim];
      Rlolo_d = new double*[dim];
   }
   dbl4_allocate_rhsqrsol(dim,deg,mode,
       rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
       rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
       urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
       urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
       Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
       Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
       solhihi_h,sollohi_h,solhilo_h,sollolo_h,
       solhihi_d,sollohi_d,solhilo_d,sollolo_d);

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
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolhihi = new double*[dim];
   double **testsollohi = new double*[dim];
   double **testsolhilo = new double*[dim];
   double **testsollolo = new double*[dim];
   double **mbrhshihi = new double*[dim];
   double **mbrhslohi = new double*[dim];
   double **mbrhshilo = new double*[dim];
   double **mbrhslolo = new double*[dim];

   dbl4_column_setup
      (dim,deg,nbrcol,nvr,idx,rowsA,cffhihi,cfflohi,cffhilo,cfflolo,
       testsolhihi,testsollohi,testsolhilo,testsollolo,
       mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,mode,vrblvl);

   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = true;
   bool zeroQ_d = true;
   bool noqr_h = false;
   bool noqr_d = false;
   int tailidx_h = 1;
   int tailidx_d = 1;
   int wrkdeg = 0; // working degree of precision
   int stepcnt = 0;

   double totcnvlapsedms = 0.0;
   double totqrlapsedms = 0.0;
   double totqtblapsedms = 0.0;
   double totbslapsedms = 0.0;
   double totupdlapsedms = 0.0;
   double totreslapsedms = 0.0;

   struct timeval begintime,endtime; // wall clock time of computations
   gettimeofday(&begintime,0);

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      dbl4_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,
          &tailidx_h,&tailidx_d,nvr,idx,exp,nbrfac,expfac,
          mbrhshihi,mbrhslohi,mbrhshilo,mbrhslolo,dpr,
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
          workvechihi,workveclohi,workvechilo,workveclolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,
          &zeroQ_h,&noqr_h,&zeroQ_d,&noqr_d,
          &upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
          &totcnvlapsedms,&totqrlapsedms,&totqtblapsedms,&totbslapsedms,
          &totupdlapsedms,&totreslapsedms,vrblvl,mode);

      stepcnt = stepcnt + 1;

      if(vrblvl > 0)
         cout << "up_h : " << upidx_h << "  bs_h : " << bsidx_h
              << "  tail_h : " << tailidx_h
              << "  up_d : " << upidx_d << "  bs_d : " << bsidx_d
              << "  tail_d : " << tailidx_d
              << "  wdeg : " << wrkdeg << endl;

      if((mode == 1) || (mode == 2)) if(tailidx_h >= deg) break;
      if((mode == 0) || (mode == 2)) if(tailidx_d >= deg) break;

      wrkdeg = wrkdeg + 1 + wrkdeg/2;
      if(wrkdeg > deg) wrkdeg = deg;
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   double walltimesec = seconds + microseconds*1.0e-6;

   dbl4_error_testsol
      (dim,deg,mode,testsolhihi,testsollohi,testsolhilo,testsollolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_dbl4_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;
   double **csthihi = new double*[dim];
   double **cstlohi = new double*[dim];
   double **csthilo = new double*[dim];
   double **cstlolo = new double*[dim];
   double ***cffhihi = new double**[dim];
   double ***cfflohi = new double**[dim];
   double ***cffhilo = new double**[dim];
   double ***cfflolo = new double**[dim];

   dbl4_make_coefficients(dim,deg,nbr,nvr,idx,
       csthihi,cstlohi,csthilo,cstlolo,
       cffhihi,cfflohi,cffhilo,cfflolo,vrblvl);
   // must randomize the leading coefficients,
   // because for exponenials all leading coefficients are one
   for(int i=0; i<dim; i++)
      for(int j=0; j<nbr[i]; j++)
         random_quad_double(&cffhihi[i][j][0],&cfflohi[i][j][0],
                            &cffhilo[i][j][0],&cfflolo[i][j][0]);

   double **inputhihi_h;
   double **inputlohi_h;
   double **inputhilo_h;
   double **inputlolo_h;
   double **inputhihi_d;
   double **inputlohi_d;
   double **inputhilo_d;
   double **inputlolo_d;
   double ***outputhihi_h;
   double ***outputlohi_h;
   double ***outputhilo_h;
   double ***outputlolo_h;
   double ***outputhihi_d;
   double ***outputlohi_d;
   double ***outputhilo_d;
   double ***outputlolo_d;
   double **funvalhihi_h;
   double **funvallohi_h;
   double **funvalhilo_h;
   double **funvallolo_h;
   double **funvalhihi_d;
   double **funvallohi_d;
   double **funvalhilo_d;
   double **funvallolo_d;
   double ***jacvalhihi_h;
   double ***jacvallohi_h;
   double ***jacvalhilo_h;
   double ***jacvallolo_h;
   double ***jacvalhihi_d;
   double ***jacvallohi_d;
   double ***jacvalhilo_d;
   double ***jacvallolo_d;

   if((mode == 1) || (mode == 2))
   {
      inputhihi_h = new double*[dim];
      inputlohi_h = new double*[dim];
      inputhilo_h = new double*[dim];
      inputlolo_h = new double*[dim];
      outputhihi_h = new double**[dim];
      outputlohi_h = new double**[dim];
      outputhilo_h = new double**[dim];
      outputlolo_h = new double**[dim];
      funvalhihi_h = new double*[dim];
      funvallohi_h = new double*[dim];
      funvalhilo_h = new double*[dim];
      funvallolo_h = new double*[dim];
      jacvalhihi_h = new double**[degp1];
      jacvallohi_h = new double**[degp1];
      jacvalhilo_h = new double**[degp1];
      jacvallolo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputhihi_d = new double*[dim];
      inputlohi_d = new double*[dim];
      inputhilo_d = new double*[dim];
      inputlolo_d = new double*[dim];
      outputhihi_d = new double**[dim];
      outputlohi_d = new double**[dim];
      outputhilo_d = new double**[dim];
      outputlolo_d = new double**[dim];
      funvalhihi_d = new double*[dim];
      funvallohi_d = new double*[dim];
      funvalhilo_d = new double*[dim];
      funvallolo_d = new double*[dim];
      jacvalhihi_d = new double**[degp1];
      jacvallohi_d = new double**[degp1];
      jacvalhilo_d = new double**[degp1];
      jacvallolo_d = new double**[degp1];
   }
   dbl4_allocate_inoutfunjac(dim,deg,mode,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
       outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
       outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
       funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
       funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
       jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
       jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhshihi_h;
   double **rhslohi_h;
   double **rhshilo_h;
   double **rhslolo_h;
   double **rhshihi_d;
   double **rhslohi_d;
   double **rhshilo_d;
   double **rhslolo_d;
   double **urhshihi_h;
   double **urhslohi_h;
   double **urhshilo_h;
   double **urhslolo_h;
   double **urhshihi_d;
   double **urhslohi_d;
   double **urhshilo_d;
   double **urhslolo_d;
   double **solhihi_h;
   double **sollohi_h;
   double **solhilo_h;
   double **sollolo_h;
   double **solhihi_d;
   double **sollohi_d;
   double **solhilo_d;
   double **sollolo_d;
   double **Qhihi_h;
   double **Qlohi_h;
   double **Qhilo_h;
   double **Qlolo_h;
   double **Qhihi_d;
   double **Qlohi_d;
   double **Qhilo_d;
   double **Qlolo_d;
   double **Rhihi_h;
   double **Rlohi_h;
   double **Rhilo_h;
   double **Rlolo_h;
   double **Rhihi_d;
   double **Rlohi_d;
   double **Rhilo_d;
   double **Rlolo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshihi_h = new double*[degp1];
      rhslohi_h = new double*[degp1];
      rhshilo_h = new double*[degp1];
      rhslolo_h = new double*[degp1];
      urhshihi_h = new double*[degp1];
      urhslohi_h = new double*[degp1];
      urhshilo_h = new double*[degp1];
      urhslolo_h = new double*[degp1];
      solhihi_h = new double*[degp1];
      sollohi_h = new double*[degp1];
      solhilo_h = new double*[degp1];
      sollolo_h = new double*[degp1];
      Qhihi_h = new double*[dim];
      Qlohi_h = new double*[dim];
      Qhilo_h = new double*[dim];
      Qlolo_h = new double*[dim];
      Rhihi_h = new double*[dim];
      Rlohi_h = new double*[dim];
      Rhilo_h = new double*[dim];
      Rlolo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhshihi_d = new double*[degp1];
      rhslohi_d = new double*[degp1];
      rhshilo_d = new double*[degp1];
      rhslolo_d = new double*[degp1];
      urhshihi_d = new double*[degp1];
      urhslohi_d = new double*[degp1];
      urhshilo_d = new double*[degp1];
      urhslolo_d = new double*[degp1];
      solhihi_d = new double*[degp1];
      sollohi_d = new double*[degp1];
      solhilo_d = new double*[degp1];
      sollolo_d = new double*[degp1];
      Qhihi_d = new double*[dim];
      Qlohi_d = new double*[dim];
      Qhilo_d = new double*[dim];
      Qlolo_d = new double*[dim];
      Rhihi_d = new double*[dim];
      Rlohi_d = new double*[dim];
      Rhilo_d = new double*[dim];
      Rlolo_d = new double*[dim];
   }
   dbl4_allocate_rhsqrsol(dim,deg,mode,
       rhshihi_h,rhslohi_h,rhshilo_h,rhslolo_h,
       rhshihi_d,rhslohi_d,rhshilo_d,rhslolo_d,
       urhshihi_h,urhslohi_h,urhshilo_h,urhslolo_h,
       urhshihi_d,urhslohi_d,urhshilo_d,urhslolo_d,
       Qhihi_h,Qlohi_h,Qhilo_h,Qlolo_h,Qhihi_d,Qlohi_d,Qhilo_d,Qlolo_d,
       Rhihi_h,Rlohi_h,Rhilo_h,Rlolo_h,Rhihi_d,Rlohi_d,Rhilo_d,Rlolo_d,
       solhihi_h,sollohi_h,solhilo_h,sollolo_h,
       solhihi_d,sollohi_d,solhilo_d,sollolo_d);
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolhihi = new double*[dim];
   double **testsollohi = new double*[dim];
   double **testsolhilo = new double*[dim];
   double **testsollolo = new double*[dim];

   dbl4_row_setup(dim,deg,nbr,nvr,idx,
      csthihi,cstlohi,csthilo,cstlolo,cffhihi,cfflohi,cffhilo,cfflolo,
      testsolhihi,testsollohi,testsolhilo,testsollolo,
      inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
      inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
      outputhihi_h,outputlohi_h,outputhilo_h,outputlolo_h,
      outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,mode,vrblvl);

   // allocating extra work space

   double *workvechihi = new double[dim];
   double *workveclohi = new double[dim];
   double *workvechilo = new double[dim];
   double *workveclolo = new double[dim];
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

   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = true;
   bool zeroQ_d = true;
   bool noqr_h = false;
   bool noqr_d = false;
   int tailidx_h = 1;
   int tailidx_d = 1;
   int wrkdeg = 0; // working degree of precision
   int stepcnt = 0;

   double totcnvlapsedms = 0.0;
   double totqrlapsedms = 0.0;
   double totqtblapsedms = 0.0;
   double totbslapsedms = 0.0;
   double totupdlapsedms = 0.0;
   double totreslapsedms = 0.0;

   struct timeval begintime,endtime; // wall clock time of computations
   gettimeofday(&begintime,0);

   for(int step=0; step<nbsteps; step++)
   {
      if(vrblvl > 0)
         cout << "*** running Newton step " << step
              << " at degree " << wrkdeg << " ***" << endl;

      dbl4_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,nbr,nvr,idx,
          csthihi,cstlohi,csthilo,cstlolo,
          cffhihi,cfflohi,cffhilo,cfflolo,dpr,
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
          workvechihi,workveclohi,workvechilo,workveclolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,
          &zeroQ_h,&noqr_h,&zeroQ_d,&noqr_d,
          &upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
          &totcnvlapsedms,&totqrlapsedms,&totqtblapsedms,&totbslapsedms,
          &totupdlapsedms,&totreslapsedms,vrblvl,mode);

      stepcnt = stepcnt + 1;

      if(vrblvl > 0)
         cout << "up_h : " << upidx_h << "  bs_h : " << bsidx_h
              << "  tail_h : " << tailidx_h
              << "  up_d : " << upidx_d << "  bs_d : " << bsidx_d
              << "  tail_d : " << tailidx_d
              << "  wdeg : " << wrkdeg << endl;

      if((mode == 1) || (mode == 2)) if(tailidx_h >= deg) break;
      if((mode == 0) || (mode == 2)) if(tailidx_d >= deg) break;

      wrkdeg = wrkdeg + 1 + wrkdeg/2;
      if(wrkdeg > deg) wrkdeg = deg;
   }
   gettimeofday(&endtime,0);
   long seconds = endtime.tv_sec - begintime.tv_sec;
   long microseconds = endtime.tv_usec - begintime.tv_usec;
   double walltimesec = seconds + microseconds*1.0e-6;

   dbl4_error_testsol
      (dim,deg,mode,testsolhihi,testsollohi,testsolhilo,testsollolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
