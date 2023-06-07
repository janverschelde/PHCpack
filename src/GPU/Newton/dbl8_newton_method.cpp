// The file dbl8_newton_method.cpp defines the functions with prototypes in
// the file dbl8_newton_method.h.

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
#include "random8_vectors.h"
#include "octo_double_functions.h"
#include "dbl8_indexed_coefficients.h"
#include "dbl8_polynomials_host.h"
#include "convolution_jobs.h"
#include "addition_jobs.h"
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"
#include "complexadd_jobs.h"
#include "job_makers.h"
#include "dbl8_polynomials_kernels.h"
#include "dbl8_convolutions_host.h"
#include "dbl8_monomials_host.h"
#include "dbl8_factorizations.h"
#include "dbl8_monomial_systems.h"
#include "dbl8_bals_host.h"
#include "dbl8_bals_kernels.h"
#include "dbl8_tail_kernels.h"
#include "dbl8_systems_host.h"
#include "dbl8_systems_kernels.h"
#include "dbl8_newton_testers.h"
#include "write_newton_times.h"

using namespace std;

int dbl8_errors_funjacrhs
 ( int dim, int deg,
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
   double **rhshilolo_d, double **rhslololo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-100;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = dbl8_error2sum(dim,degp1,
               funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
               funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
               funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
               funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
               "funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = dbl8_error3sum(degp1,dim,dim,
               jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
               jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
               jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
               jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
               "jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   errsum = dbl8_error2sum(degp1,dim,
               rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
               rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
               rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
               rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
               "rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl8_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
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
   double **solhilolo_d, double **sollololo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-100;
   double errsum = 0.0;
   int fail = 0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   errsum = dbl8_error2sum(dim,dim,
               Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
               Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
               Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
               Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = dbl8_error2sum(dim,dim,
               Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
               Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
               Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
               Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = dbl8_error2sum(degp1,dim,
               urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
               urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
               urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
               urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
               "urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = dbl8_error2sum(degp1,dim,
               solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
               solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
               solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
               solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
               "sol",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = dbl8_error2sum(dim,degp1,
               inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
               inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
               inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
               inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
               "input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int dbl8_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d,
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
   double *resmaxhilolo, double *resmaxlololo,
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
            urhshihihi_h[i][j] = rhshihihi_h[i][j];
            urhslohihi_h[i][j] = rhslohihi_h[i][j];
            urhshilohi_h[i][j] = rhshilohi_h[i][j];
            urhslolohi_h[i][j] = rhslolohi_h[i][j];
            urhshihilo_h[i][j] = rhshihilo_h[i][j];
            urhslohilo_h[i][j] = rhslohilo_h[i][j];
            urhshilolo_h[i][j] = rhshilolo_h[i][j];
            urhslololo_h[i][j] = rhslololo_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhihihi_h[i][j] = 0.0; sollohihi_h[i][j] = 0.0;
            solhilohi_h[i][j] = 0.0; sollolohi_h[i][j] = 0.0;
            solhihilo_h[i][j] = 0.0; sollohilo_h[i][j] = 0.0;
            solhilolo_h[i][j] = 0.0; sollololo_h[i][j] = 0.0;
         }
 
      if(vrblvl > 0) cout << "calling CPU_dbl8_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_dbl8_qrbs_solve
         (dim,degp1,oldtail,
          jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
          jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
          urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
          urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
          Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
          Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
          Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
          Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl8_linear_residue ..." << endl;
         CPU_dbl8_linear_residue
            (dim,degp1,*tailidx_h-1,
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

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhihihi << endl;
      }
      dbl8_update_series
         (dim,degp1,*tailidx_h-1,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
          solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
          solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhshihihi_d[i][j] = rhshihihi_d[i][j];
            urhslohihi_d[i][j] = rhslohihi_d[i][j];
            urhshilohi_d[i][j] = rhshilohi_d[i][j];
            urhslolohi_d[i][j] = rhslolohi_d[i][j];
            urhshihilo_d[i][j] = rhshihilo_d[i][j];
            urhslohilo_d[i][j] = rhslohilo_d[i][j];
            urhshilolo_d[i][j] = rhshilolo_d[i][j];
            urhslololo_d[i][j] = rhslololo_d[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhihihi_d[i][j] = 0.0; sollohihi_d[i][j] = 0.0;
            solhilohi_d[i][j] = 0.0; sollolohi_d[i][j] = 0.0;
            solhihilo_d[i][j] = 0.0; sollohilo_d[i][j] = 0.0;
            solhilolo_d[i][j] = 0.0; sollololo_d[i][j] = 0.0;
         }
 
      if(vrblvl > 0) cout << "calling GPU_dbl8_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_dbl8_bals_solve
         (dim,degp1,szt,nbt,oldtail,
          jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
          jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
          Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
          Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
          Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
          Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
          urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
          urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,zeroQ_d,noqr_d,
          upidx_d,bsidx_d,&newtail,totqrlapsedms,totqtblapsedms,
          totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl8_linear_residue ..." << endl;
         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_dbl8_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,
             jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
             jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
             rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
             rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
             solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
             solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,
             resvechihihi,resveclohihi,resvechilohi,resveclolohi,
             resvechihilo,resveclohilo,resvechilolo,resveclololo,
             resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
             resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,
             &elapsedms,&addcnt,&mulcnt,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhihihi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;

         *totreslapsedms += elapsedms;
      }
      dbl8_update_series
         (dim,degp1,*tailidx_d-1,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
          solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
          solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);
   }
   return 0;
}

int dbl8_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbhihihi, double **mblohihi, double **mbhilohi, double **mblolohi,
   double **mbhihilo, double **mblohilo, double **mbhilolo, double **mblololo,
   double dpr,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **acchihihi, double **acclohihi,
   double **acchilohi, double **acclolohi,
   double **acchihilo, double **acclohilo,
   double **acchilolo, double **acclololo,
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
   double *resmaxhilolo, double *resmaxlololo,
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
         cout << "calling CPU_dbl8_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl8_unit_series_vector
            (dim,deg,cffhihihi[0],cfflohihi[0],cffhilohi[0],cfflolohi[0],
                     cffhihilo[0],cfflohilo[0],cffhilolo[0],cfflololo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_dbl8_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhihihi[0],cfflohihi[0],cffhilohi[0],cfflolohi[0],
             cffhihilo[0],cfflohilo[0],cffhilolo[0],cfflololo[0],
             acchihihi[0],acclohihi[0],acchilohi[0],acclolohi[0],
             acchihilo[0],acclohilo[0],acchilolo[0],acclololo[0],
             inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
             inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
             outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
             outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
             vrblvl);
      }
      else
         CPU_dbl8_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,
             cffhihihi,cfflohihi,cffhilohi,cfflolohi,
             cffhihilo,cfflohilo,cffhilolo,cfflololo,
             acchihihi,acclohihi,acchilohi,acclolohi,
             acchihilo,acclohilo,acchilolo,acclololo,
             inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
             inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
             funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
             funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
             jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
             jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
             vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl8_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl8_unit_series_vector
            (dim,deg,cffhihihi[0],cfflohihi[0],cffhilohi[0],cfflolohi[0],
                     cffhihilo[0],cfflohilo[0],cffhilolo[0],cfflololo[0]);
         // reset the coefficients

         GPU_dbl8_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhihihi[0],cfflohihi[0],cffhilohi[0],cfflolohi[0],
             cffhihilo[0],cfflohilo[0],cffhilolo[0],cfflololo[0],
             acchihihi[0],acclohihi[0],acchilohi[0],acclolohi[0],
             acchihilo[0],acclohilo[0],acchilolo[0],acclololo[0],
             inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
             inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
             outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
             outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
             totcnvlapsedms,vrblvl);
      }
      else
         GPU_dbl8_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,
             cffhihihi,cfflohihi,cffhilohi,cfflolohi,
             cffhihilo,cfflohilo,cffhilolo,cfflololo,
             inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
             inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
             outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
             outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
             funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
             funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
             jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
             jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
             totcnvlapsedms,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;

      double errsum = 0.0;

      errsum = dbl8_error3sum(dim,dim+1,degp1,
                  outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
                  outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
                  outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
                  outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
                  "output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(3)
           << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;

   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl8_define_rhs
            (dim,degp1,
             mbhihihi,mblohihi,mbhilohi,mblolohi,
             mbhihilo,mblohilo,mbhilolo,mblololo,
             funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
             funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
             rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
             rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalhihihi_h[i][j][k] = 0.0;
                  jacvallohihi_h[i][j][k] = 0.0;
                  jacvalhilohi_h[i][j][k] = 0.0;
                  jacvallolohi_h[i][j][k] = 0.0;
                  jacvalhihilo_h[i][j][k] = 0.0;
                  jacvallohilo_h[i][j][k] = 0.0;
                  jacvalhilolo_h[i][j][k] = 0.0;
                  jacvallololo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl8_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],
             mbhihihi,mblohihi,mbhilohi,mblolohi,
             mbhihilo,mblohilo,mbhilolo,mblololo,dpr,
             outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
             outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
             funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
             funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
             rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
             rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
             jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
             jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
             vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl8_define_rhs
            (dim,degp1,
             mbhihihi,mblohihi,mbhilohi,mblolohi,
             mbhihilo,mblohilo,mbhilolo,mblololo,
             funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
             funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
             rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
             rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalhihihi_d[i][j][k] = 0.0;
                  jacvallohihi_d[i][j][k] = 0.0;
                  jacvalhilohi_d[i][j][k] = 0.0;
                  jacvallolohi_d[i][j][k] = 0.0;
                  jacvalhihilo_d[i][j][k] = 0.0;
                  jacvallohilo_d[i][j][k] = 0.0;
                  jacvalhilolo_d[i][j][k] = 0.0;
                  jacvallololo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl8_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],
             mbhihihi,mblohihi,mbhilohi,mblolohi,
             mbhihilo,mblohilo,mbhilolo,mblololo,dpr,
             outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
             outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
             funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
             funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
             rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
             rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
             jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
             jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
             vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl8_errors_funjacrhs(dim,deg,
                    funvalhihihi_h,funvallohihi_h,
                    funvalhilohi_h,funvallolohi_h,
                    funvalhihilo_h,funvallohilo_h,
                    funvalhilolo_h,funvallololo_h,
                    funvalhihihi_d,funvallohihi_d,
                    funvalhilohi_d,funvallolohi_d,
                    funvalhihilo_d,funvallohilo_d,
                    funvalhilolo_d,funvallololo_d,
                    jacvalhihihi_h,jacvallohihi_h,
                    jacvalhilohi_h,jacvallolohi_h,
                    jacvalhihilo_h,jacvallohilo_h,
                    jacvalhilolo_h,jacvallololo_h,
                    jacvalhihihi_d,jacvallohihi_d,
                    jacvalhilohi_d,jacvallolohi_d,
                    jacvalhihilo_d,jacvallohilo_d,
                    jacvalhilolo_d,jacvallololo_d,
                    rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
                    rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
                    rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
                    rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,vrblvl);
   }
   dbl8_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
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
       workvechihihi,workveclohihi,workvechilohi,workveclolohi,
       workvechihilo,workveclohilo,workvechilolo,workveclololo,
       resvechihihi,resveclohihi,resvechilohi,resveclolohi,
       resvechihilo,resveclohilo,resvechilolo,resveclololo,
       resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
       resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return dbl8_errors_inurhsQRsol(dim,deg,
                 inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
                 inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
                 inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
                 inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
                 Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
                 Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
                 Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
                 Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
                 Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
                 Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
                 Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
                 Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
                 urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
                 urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
                 urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
                 urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
                 solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
                 solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
                 solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
                 solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);
   }
   else
      return 0;
}

int dbl8_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx, 
   double **csthihihi, double **cstlohihi,
   double **csthilohi, double **cstlolohi,
   double **csthihilo, double **cstlohilo,
   double **csthilolo, double **cstlololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo, double dpr,
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
   double *resmaxhilolo, double *resmaxlololo,
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
         cout << "calling CPU_dbl8_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl8_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
            csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
            csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
            cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
            cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
            inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
            inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
            outputhihihi_h[i],outputlohihi_h[i],
            outputhilohi_h[i],outputlolohi_h[i],
            outputhihilo_h[i],outputlohilo_h[i],
            outputhilolo_h[i],outputlololo_h[i],&lapsed,0);

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
         cout << "calling GPU_dbl8_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ConvolutionJobs cnvjobs(dim);
         AdditionJobs addjobs(dim,nbr[i]);

         make_all_jobs(dim,nbr[i],nvr[i],idx[i],&cnvjobs,&addjobs,false);

         GPU_dbl8_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
             csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
             cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
             cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
             inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
             inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
             outputhihihi_d[i],outputlohihi_d[i],
             outputhilohi_d[i],outputlolohi_d[i],
             outputhihilo_d[i],outputlohilo_d[i],
             outputhilolo_d[i],outputlololo_d[i],
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

      errsum = dbl8_error3sum(dim,dim+1,degp1,
                  outputhihihi_h,outputlohihi_h,
                  outputhilohi_h,outputlolohi_h,
                  outputhihilo_h,outputlohilo_h,
                  outputhilolo_h,outputlololo_h,
                  outputhihihi_d,outputlohihi_d,
                  outputhilohi_d,outputlolohi_d,
                  outputhihilo_d,outputlohilo_d,
                  outputhilolo_d,outputlololo_d,
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
      dbl8_map_evaldiff_output(dim,deg,
         outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
         outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
         funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
         funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
         jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
         jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
         vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshihihi_h[i][j] = -funvalhihihi_h[j][i];
            rhslohihi_h[i][j] = -funvallohihi_h[j][i];
            rhshilohi_h[i][j] = -funvalhilohi_h[j][i];
            rhslolohi_h[i][j] = -funvallolohi_h[j][i];
            rhshihilo_h[i][j] = -funvalhihilo_h[j][i];
            rhslohilo_h[i][j] = -funvallohilo_h[j][i];
            rhshilolo_h[i][j] = -funvalhilolo_h[j][i];
            rhslololo_h[i][j] = -funvallololo_h[j][i];
         }
   }
   if((mode == 0) || (mode == 2))
   {
      dbl8_map_evaldiff_output(dim,deg,
         outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
         outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
         funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
         funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
         jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
         jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
         vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhshihihi_d[i][j] = -funvalhihihi_d[j][i];
            rhslohihi_d[i][j] = -funvallohihi_d[j][i];
            rhshilohi_d[i][j] = -funvalhilohi_d[j][i];
            rhslolohi_d[i][j] = -funvallolohi_d[j][i];
            rhshihilo_d[i][j] = -funvalhihilo_d[j][i];
            rhslohilo_d[i][j] = -funvallohilo_d[j][i];
            rhshilolo_d[i][j] = -funvalhilolo_d[j][i];
            rhslololo_d[i][j] = -funvallololo_d[j][i];
         }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = dbl8_errors_funjacrhs(dim,deg,
                    funvalhihihi_h,funvallohihi_h,
                    funvalhilohi_h,funvallolohi_h,
                    funvalhihilo_h,funvallohilo_h,
                    funvalhilolo_h,funvallololo_h,
                    funvalhihihi_d,funvallohihi_d,
                    funvalhilohi_d,funvallolohi_d,
                    funvalhihilo_d,funvallohilo_d,
                    funvalhilolo_d,funvallololo_d,
                    jacvalhihihi_h,jacvallohihi_h,
                    jacvalhilohi_h,jacvallolohi_h,
                    jacvalhihilo_h,jacvallohilo_h,
                    jacvalhilolo_h,jacvallololo_h,
                    jacvalhihihi_d,jacvallohihi_d,
                    jacvalhilohi_d,jacvallolohi_d,
                    jacvalhihilo_d,jacvallohilo_d,
                    jacvalhilolo_d,jacvallololo_d,
                    rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
                    rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
                    rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
                    rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,vrblvl);
   }
   dbl8_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
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
       workvechihihi,workveclohihi,workvechilohi,workveclolohi,
       workvechihilo,workveclohilo,workvechilolo,workveclololo,
       resvechihihi,resveclohihi,resvechilohi,resveclolohi,
       resvechihilo,resveclohilo,resvechilolo,resveclololo,
       resmaxhihihi,resmaxlohihi,resmaxhilohi,resmaxlolohi,
       resmaxhihilo,resmaxlohilo,resmaxhilolo,resmaxlololo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return dbl8_errors_inurhsQRsol(dim,deg,
                 inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
                 inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
                 inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
                 inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
                 Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
                 Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
                 Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
                 Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
                 Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
                 Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
                 Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
                 Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
                 urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
                 urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
                 urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
                 urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
                 solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
                 solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
                 solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
                 solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d,vrblvl);
   }
   else
      return 0;
}

int dbl8_allocate_inoutfunjac
 ( int dim, int deg, int mode,
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
   double ***jacvalhilolo_d, double ***jacvallololo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
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

         outputhihihi_h[i] = new double*[dim+1];
         outputlohihi_h[i] = new double*[dim+1];
         outputhilohi_h[i] = new double*[dim+1];
         outputlolohi_h[i] = new double*[dim+1];
         outputhihilo_h[i] = new double*[dim+1];
         outputlohilo_h[i] = new double*[dim+1];
         outputhilolo_h[i] = new double*[dim+1];
         outputlololo_h[i] = new double*[dim+1];
 
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
         }
         funvalhihihi_h[i] = new double[degp1];
         funvallohihi_h[i] = new double[degp1];
         funvalhilohi_h[i] = new double[degp1];
         funvallolohi_h[i] = new double[degp1];
         funvalhihilo_h[i] = new double[degp1];
         funvallohilo_h[i] = new double[degp1];
         funvalhilolo_h[i] = new double[degp1];
         funvallololo_h[i] = new double[degp1];
      }
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
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputhihihi_d[i] = new double[degp1];
         inputlohihi_d[i] = new double[degp1];
         inputhilohi_d[i] = new double[degp1];
         inputlolohi_d[i] = new double[degp1];
         inputhihilo_d[i] = new double[degp1];
         inputlohilo_d[i] = new double[degp1];
         inputhilolo_d[i] = new double[degp1];
         inputlololo_d[i] = new double[degp1];

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
            outputhihihi_d[i][j] = new double[degp1];
            outputlohihi_d[i][j] = new double[degp1];
            outputhilohi_d[i][j] = new double[degp1];
            outputlolohi_d[i][j] = new double[degp1];
            outputhihilo_d[i][j] = new double[degp1];
            outputlohilo_d[i][j] = new double[degp1];
            outputhilolo_d[i][j] = new double[degp1];
            outputlololo_d[i][j] = new double[degp1];
         }
         funvalhihihi_d[i] = new double[degp1];
         funvallohihi_d[i] = new double[degp1];
         funvalhilohi_d[i] = new double[degp1];
         funvallolohi_d[i] = new double[degp1];
         funvalhihilo_d[i] = new double[degp1];
         funvallohilo_d[i] = new double[degp1];
         funvalhilolo_d[i] = new double[degp1];
         funvallololo_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
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
   }
   return 0;
}

int dbl8_allocate_rhsqrsol
 ( int dim, int deg, int mode,
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
   double **solhihihi_h, double **sollohihi_h,
   double **solhilohi_h, double **sollolohi_h,
   double **solhihilo_h, double **sollohilo_h,
   double **solhilolo_h, double **sollololo_h,
   double **solhihihi_d, double **sollohihi_d,
   double **solhilohi_d, double **sollolohi_d,
   double **solhihilo_d, double **sollohilo_d,
   double **solhilolo_d, double **sollololo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
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
         urhshihihi_h[i] = new double[dim];
         urhslohihi_h[i] = new double[dim];
         urhshilohi_h[i] = new double[dim];
         urhslolohi_h[i] = new double[dim];
         urhshihilo_h[i] = new double[dim];
         urhslohilo_h[i] = new double[dim];
         urhshilolo_h[i] = new double[dim];
         urhslololo_h[i] = new double[dim];
         solhihihi_h[i] = new double[dim];
         sollohihi_h[i] = new double[dim];
         solhilohi_h[i] = new double[dim];
         sollolohi_h[i] = new double[dim];
         solhihilo_h[i] = new double[dim];
         sollohilo_h[i] = new double[dim];
         solhilolo_h[i] = new double[dim];
         sollololo_h[i] = new double[dim];
      }
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
         Rhihihi_h[i] = new double[dim];
         Rlohihi_h[i] = new double[dim];
         Rhilohi_h[i] = new double[dim];
         Rlolohi_h[i] = new double[dim];
         Rhihilo_h[i] = new double[dim];
         Rlohilo_h[i] = new double[dim];
         Rhilolo_h[i] = new double[dim];
         Rlololo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++)
      {
         rhshihihi_d[i] = new double[dim];
         rhslohihi_d[i] = new double[dim];
         rhshilohi_d[i] = new double[dim];
         rhslolohi_d[i] = new double[dim];
         rhshihilo_d[i] = new double[dim];
         rhslohilo_d[i] = new double[dim];
         rhshilolo_d[i] = new double[dim];
         rhslololo_d[i] = new double[dim];
         urhshihihi_d[i] = new double[dim];
         urhslohihi_d[i] = new double[dim];
         urhshilohi_d[i] = new double[dim];
         urhslolohi_d[i] = new double[dim];
         urhshihilo_d[i] = new double[dim];
         urhslohilo_d[i] = new double[dim];
         urhshilolo_d[i] = new double[dim];
         urhslololo_d[i] = new double[dim];
         solhihihi_d[i] = new double[dim];
         sollohihi_d[i] = new double[dim];
         solhilohi_d[i] = new double[dim];
         sollolohi_d[i] = new double[dim];
         solhihilo_d[i] = new double[dim];
         sollohilo_d[i] = new double[dim];
         solhilolo_d[i] = new double[dim];
         sollololo_d[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qhihihi_d[i] = new double[dim];
         Qlohihi_d[i] = new double[dim];
         Qhilohi_d[i] = new double[dim];
         Qlolohi_d[i] = new double[dim];
         Qhihilo_d[i] = new double[dim];
         Qlohilo_d[i] = new double[dim];
         Qhilolo_d[i] = new double[dim];
         Qlololo_d[i] = new double[dim];
         Rhihihi_d[i] = new double[dim];
         Rlohihi_d[i] = new double[dim];
         Rhilohi_d[i] = new double[dim];
         Rlolohi_d[i] = new double[dim];
         Rhihilo_d[i] = new double[dim];
         Rlohilo_d[i] = new double[dim];
         Rhilolo_d[i] = new double[dim];
         Rlololo_d[i] = new double[dim];
      }
   }
   return 0;
}

void dbl8_start_setup
 ( int dim, int deg,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d, int mode, int vrblvl )
{
   double *start0hihihi = new double[dim];
   double *start0lohihi = new double[dim];
   double *start0hilohi = new double[dim];
   double *start0lolohi = new double[dim];
   double *start0hihilo = new double[dim];
   double *start0lohilo = new double[dim];
   double *start0hilolo = new double[dim];
   double *start0lololo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0hihihi[i] = testsolhihihi[i][0];
      start0lohihi[i] = testsollohihi[i][0];
      start0hilohi[i] = testsolhilohi[i][0];
      start0lolohi[i] = testsollolohi[i][0];
      start0hihilo[i] = testsolhihilo[i][0];
      start0lohilo[i] = testsollohilo[i][0];
      start0hilolo[i] = testsolhilolo[i][0];
      start0lololo[i] = testsollololo[i][0];
   }
   if((mode == 1) || (mode == 2))
      real8_start_series_vector
         (dim,deg,start0hihihi,start0lohihi,start0hilohi,start0lolohi,
                  start0hihilo,start0lohilo,start0hilolo,start0lololo,
          inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
          inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h);
   else
      real8_start_series_vector
         (dim,deg,start0hihihi,start0lohihi,start0hilohi,start0lolohi,
                  start0hihilo,start0lohilo,start0hilolo,start0lololo,
          inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
          inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
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
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;

      if((mode == 1) || (mode == 2))
      {
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
      else
      {
         for(int i=0; i<dim; i++)
         {
            cout << i << " : " << inputhihihi_d[i][0] << "  "
                               << inputlohihi_d[i][0] << endl;
            cout << "     " << inputhilohi_d[i][0] << "  "
                            << inputlolohi_d[i][0] << endl;
            cout << "     " << inputhihilo_d[i][0] << "  "
                            << inputlohilo_d[i][0] << endl;
            cout << "     " << inputhilolo_d[i][0] << "  "
                            << inputlololo_d[i][0] << endl;
         }
      }
   }
}

void dbl8_column_setup
 ( int dim, int deg, int nbrcol, int **nvr, int ***idx, int **rowsA,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
   double **mbrhshihihi, double **mbrhslohihi,
   double **mbrhshilohi, double **mbrhslolohi,
   double **mbrhshihilo, double **mbrhslohilo,
   double **mbrhshilolo, double **mbrhslololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolhihihi[i] = new double[degp1];
      testsollohihi[i] = new double[degp1];
      testsolhilohi[i] = new double[degp1];
      testsollolohi[i] = new double[degp1];
      testsolhihilo[i] = new double[degp1];
      testsollohilo[i] = new double[degp1];
      testsolhilolo[i] = new double[degp1];
      testsollololo[i] = new double[degp1];
   }
   make_real8_exponentials
      (dim,deg,testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
               testsolhihilo,testsollohilo,testsolhilolo,testsollololo);

   // compute the right hand sides via evaluation

   for(int i=0; i<dim; i++)
   {
      mbrhshihihi[i] = new double[degp1];
      mbrhslohihi[i] = new double[degp1];
      mbrhshilohi[i] = new double[degp1];
      mbrhslolohi[i] = new double[degp1];
      mbrhshihilo[i] = new double[degp1];
      mbrhslohilo[i] = new double[degp1];
      mbrhshilolo[i] = new double[degp1];
      mbrhslololo[i] = new double[degp1];

      mbrhshihihi[i][0] = 1.0;     // initialize product to one
      mbrhslohihi[i][0] = 0.0;
      mbrhshilohi[i][0] = 0.0;
      mbrhslolohi[i][0] = 0.0;
      mbrhshihilo[i][0] = 0.0; 
      mbrhslohilo[i][0] = 0.0;
      mbrhshilolo[i][0] = 0.0;
      mbrhslololo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhshihihi[i][k] = 0.0; mbrhslohihi[i][k] = 0.0;
         mbrhshilohi[i][k] = 0.0; mbrhslolohi[i][k] = 0.0;
         mbrhshihilo[i][k] = 0.0; mbrhslohilo[i][k] = 0.0;
         mbrhshilolo[i][k] = 0.0; mbrhslololo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_real8_monomials
         (dim,deg,rowsA,
          testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
          testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
          mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
          mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo);
   else
      evaluate_real8_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
          testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
          mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
          mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshihihi[i][j] << "  " << mbrhslohihi[i][j] << endl
                 << "  "
                 << mbrhshilohi[i][j] << "  " << mbrhslolohi[i][j] << endl
                 << "  "
                 << mbrhshihilo[i][j] << "  " << mbrhslohilo[i][j] << endl
                 << "  "
                 << mbrhshilolo[i][j] << "  " << mbrhslololo[i][j] << endl;
   }
   dbl8_start_setup(dim,deg,
      testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
      testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
      inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
      inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
      inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
      inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,mode,vrblvl);
}

void dbl8_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **csthihihi, double **cstlohihi,
   double **csthilohi, double **cstlolohi,
   double **csthihilo, double **cstlohilo,
   double **csthilolo, double **cstlololo,
   double ***cffhihihi, double ***cfflohihi,
   double ***cffhilohi, double ***cfflolohi,
   double ***cffhihilo, double ***cfflohilo,
   double ***cffhilolo, double ***cfflololo,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
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
   double ***outputhilolo_d, double ***outputlololo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++) 
   {
      testsolhihihi[i] = new double[degp1];
      testsollohihi[i] = new double[degp1];
      testsolhilohi[i] = new double[degp1];
      testsollolohi[i] = new double[degp1];
      testsolhihilo[i] = new double[degp1];
      testsollohilo[i] = new double[degp1];
      testsolhilolo[i] = new double[degp1];
      testsollololo[i] = new double[degp1];
   }
   make_real8_exponentials(dim,deg,
      testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
      testsolhihilo,testsollohilo,testsolhilolo,testsollololo);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_dbl8_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
            csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
            csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
            cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
            cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
            testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
            testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
            outputhihihi_h[i],outputlohihi_h[i],
            outputhilohi_h[i],outputlolohi_h[i],
            outputhihilo_h[i],outputlohilo_h[i],
            outputhilolo_h[i],outputlololo_h[i],&lapsed,0);

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
            odf_dec(&csthihihi[i][j],&cstlohihi[i][j],
                    &csthilohi[i][j],&cstlolohi[i][j],
                    &csthihilo[i][j],&cstlohilo[i][j],
                    &csthilolo[i][j],&cstlololo[i][j],
                    outputhihihi_h[i][dim][j],outputlohihi_h[i][dim][j],
                    outputhilohi_h[i][dim][j],outputlolohi_h[i][dim][j],
                    outputhihilo_h[i][dim][j],outputlohilo_h[i][dim][j],
                    outputhilolo_h[i][dim][j],outputlololo_h[i][dim][j]);
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_dbl8_poly_evaldiff(dim,nbr[i],deg,nvr[i],idx[i],
               csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
               csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
               cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
               cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
               testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
               testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
               outputhihihi_h[i],outputlohihi_h[i],
               outputhilohi_h[i],outputlolohi_h[i],
               outputhihilo_h[i],outputlohilo_h[i],
               outputhilolo_h[i],outputlololo_h[i],&lapsed,0);

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
               errsum = errsum + outputhihihi_h[i][dim][j]
                               + outputlohihi_h[i][dim][j]
                               + outputhilohi_h[i][dim][j]
                               + outputlolohi_h[i][dim][j]
                               + outputhihilo_h[i][dim][j]
                               + outputlohilo_h[i][dim][j]
                               + outputhilolo_h[i][dim][j]
                               + outputlololo_h[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
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

         GPU_dbl8_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
             csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
             cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
             cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
             testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
             testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
             outputhihihi_d[i],outputlohihi_d[i],
             outputhilohi_d[i],outputlolohi_d[i],
             outputhihilo_d[i],outputlohilo_d[i],
             outputhilolo_d[i],outputlololo_d[i],
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
            odf_dec(&csthihihi[i][j],&cstlohihi[i][j],
                    &csthilohi[i][j],&cstlolohi[i][j],
                    &csthihilo[i][j],&cstlohilo[i][j],
                    &csthilolo[i][j],&cstlololo[i][j],
                    outputhihihi_d[i][dim][j],outputlohihi_d[i][dim][j],
                    outputhilohi_d[i][dim][j],outputlolohi_d[i][dim][j],
                    outputhihilo_d[i][dim][j],outputlohilo_d[i][dim][j],
                    outputhilolo_d[i][dim][j],outputlololo_d[i][dim][j]);
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

            GPU_dbl8_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],
                csthihihi[i],cstlohihi[i],csthilohi[i],cstlolohi[i],
                csthihilo[i],cstlohilo[i],csthilolo[i],cstlololo[i],
                cffhihihi[i],cfflohihi[i],cffhilohi[i],cfflolohi[i],
                cffhihilo[i],cfflohilo[i],cffhilolo[i],cfflololo[i],
                testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
                testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
                outputhihihi_d[i],outputlohihi_d[i],
                outputhilohi_d[i],outputlolohi_d[i],
                outputhihilo_d[i],outputlohilo_d[i],
                outputhilolo_d[i],outputlololo_d[i],cnvjobs,addjobs,
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
               errsum = errsum + outputhihihi_d[i][dim][j]
                               + outputlohihi_d[i][dim][j]
                               + outputhilohi_d[i][dim][j]
                               + outputlolohi_d[i][dim][j]
                               + outputhihilo_d[i][dim][j]
                               + outputlohilo_d[i][dim][j]
                               + outputhilolo_d[i][dim][j]
                               + outputlololo_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
      }
   }
   dbl8_start_setup(dim,deg,
       testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
       testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,mode,vrblvl);
}

int dbl8_error_testsol
 ( int dim, int deg, int mode,
   double **testsolhihihi, double **testsollohihi,
   double **testsolhilohi, double **testsollolohi,
   double **testsolhihilo, double **testsollohilo,
   double **testsolhilolo, double **testsollololo,
   double **inputhihihi_h, double **inputlohihi_h,
   double **inputhilohi_h, double **inputlolohi_h,
   double **inputhihilo_h, double **inputlohilo_h,
   double **inputhilolo_h, double **inputlololo_h,
   double **inputhihihi_d, double **inputlohihi_d,
   double **inputhilohi_d, double **inputlolohi_d,
   double **inputhihilo_d, double **inputlohilo_d,
   double **inputhilolo_d, double **inputlololo_d )
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
                        << testsolhihihi[i][j] << "  "
                        << testsollohihi[i][j] << endl << "  "
                        << testsolhilohi[i][j] << "  "
                        << testsollolohi[i][j] << endl << "  "
                        << testsolhihilo[i][j] << "  "
                        << testsollohilo[i][j] << endl << "  "
                        << testsolhilolo[i][j] << "  "
                        << testsollololo[i][j] << endl;

         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhihihi_d[i][j] << "  "
                           << inputlohihi_d[i][j] << endl << "  "
                           << inputhilohi_d[i][j] << "  "
                           << inputlolohi_d[i][j] << endl << "  "
                           << inputhihilo_d[i][j] << "  "
                           << inputlohilo_d[i][j] << endl << "  "
                           << inputhilolo_d[i][j] << "  "
                           << inputlololo_d[i][j] << endl;

            errsum += abs(testsolhihihi[i][j] - inputhihihi_d[i][j])
                    + abs(testsollohihi[i][j] - inputlohihi_d[i][j])
                    + abs(testsolhilohi[i][j] - inputhilohi_d[i][j])
                    + abs(testsollolohi[i][j] - inputlolohi_d[i][j])
                    + abs(testsolhihilo[i][j] - inputhihilo_d[i][j])
                    + abs(testsollohilo[i][j] - inputlohilo_d[i][j])
                    + abs(testsolhilolo[i][j] - inputhilolo_d[i][j])
                    + abs(testsollololo[i][j] - inputlololo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhihihi_h[i][j] << "  "
                           << inputlohihi_h[i][j] << endl << "  "
                           << inputhilohi_h[i][j] << "  "
                           << inputlolohi_h[i][j] << endl << "  "
                           << inputhihilo_h[i][j] << "  "
                           << inputlohilo_h[i][j] << endl << "  "
                           << inputhilolo_h[i][j] << "  "
                           << inputlololo_h[i][j] << endl;

            errsum += abs(testsolhihihi[i][j] - inputhihihi_h[i][j])
                    + abs(testsollohihi[i][j] - inputlohihi_h[i][j])
                    + abs(testsolhilohi[i][j] - inputhilohi_h[i][j])
                    + abs(testsollolohi[i][j] - inputlolohi_h[i][j])
                    + abs(testsolhihilo[i][j] - inputhihilo_h[i][j])
                    + abs(testsollohilo[i][j] - inputlohilo_h[i][j])
                    + abs(testsolhilolo[i][j] - inputhilolo_h[i][j])
                    + abs(testsollololo[i][j] - inputlololo_h[i][j]);
         }
      }
   }
   cout << scientific << setprecision(2)
        << "error : " << errsum << endl;

   return (errsum > 1.0e-100);
}

int test_dbl8_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **acchihihi = new double*[dim+1]; // accumulated power series
   double **acclohihi = new double*[dim+1];
   double **acchilohi = new double*[dim+1];
   double **acclolohi = new double*[dim+1];
   double **acchihilo = new double*[dim+1];
   double **acclohilo = new double*[dim+1];
   double **acchilolo = new double*[dim+1];
   double **acclololo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchihihi[i] = new double[degp1];
      acclohihi[i] = new double[degp1];
      acchilohi[i] = new double[degp1];
      acclolohi[i] = new double[degp1];
      acchihilo[i] = new double[degp1];
      acclohilo[i] = new double[degp1];
      acchilolo[i] = new double[degp1];
      acclololo[i] = new double[degp1];
   }
   double ***cffhihihi = new double**[nbrcol]; // coefficients of monomials
   double ***cfflohihi = new double**[nbrcol];
   double ***cffhilohi = new double**[nbrcol];
   double ***cfflolohi = new double**[nbrcol];
   double ***cffhihilo = new double**[nbrcol];
   double ***cfflohilo = new double**[nbrcol];
   double ***cffhilolo = new double**[nbrcol];
   double ***cfflololo = new double**[nbrcol];

   for(int i=0; i<nbrcol; i++)
   {
      cffhihihi[i] = new double*[dim];
      cfflohihi[i] = new double*[dim];
      cffhilohi[i] = new double*[dim];
      cfflolohi[i] = new double*[dim];
      cffhihilo[i] = new double*[dim];
      cfflohilo[i] = new double*[dim];
      cffhilolo[i] = new double*[dim];
      cfflololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffhihihi[i][j] = new double[degp1];
         cfflohihi[i][j] = new double[degp1];
         cffhilohi[i][j] = new double[degp1];
         cfflolohi[i][j] = new double[degp1];
         cffhihilo[i][j] = new double[degp1];
         cfflohilo[i][j] = new double[degp1];
         cffhilolo[i][j] = new double[degp1];
         cfflololo[i][j] = new double[degp1];
      }
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_real8_coefficients
         (nbrcol,dim,cffhihihi,cfflohihi,cffhilohi,cfflolohi,
                     cffhihilo,cfflohilo,cffhilolo,cfflololo);

   double **inputhihihi_h;
   double **inputlohihi_h;
   double **inputhilohi_h;
   double **inputlolohi_h;
   double **inputhihilo_h;
   double **inputlohilo_h;
   double **inputhilolo_h;
   double **inputlololo_h;
   double **inputhihihi_d;
   double **inputlohihi_d;
   double **inputhilohi_d;
   double **inputlolohi_d;
   double **inputhihilo_d;
   double **inputlohilo_d;
   double **inputhilolo_d;
   double **inputlololo_d;
   double ***outputhihihi_h;
   double ***outputlohihi_h;
   double ***outputhilohi_h;
   double ***outputlolohi_h;
   double ***outputhihilo_h;
   double ***outputlohilo_h;
   double ***outputhilolo_h;
   double ***outputlololo_h;
   double ***outputhihihi_d;
   double ***outputlohihi_d;
   double ***outputhilohi_d;
   double ***outputlolohi_d;
   double ***outputhihilo_d;
   double ***outputlohilo_d;
   double ***outputhilolo_d;
   double ***outputlololo_d;
   double **funvalhihihi_h;
   double **funvallohihi_h;
   double **funvalhilohi_h;
   double **funvallolohi_h;
   double **funvalhihilo_h;
   double **funvallohilo_h;
   double **funvalhilolo_h;
   double **funvallololo_h;
   double **funvalhihihi_d;
   double **funvallohihi_d;
   double **funvalhilohi_d;
   double **funvallolohi_d;
   double **funvalhihilo_d;
   double **funvallohilo_d;
   double **funvalhilolo_d;
   double **funvallololo_d;
   double ***jacvalhihihi_h;
   double ***jacvallohihi_h;
   double ***jacvalhilohi_h;
   double ***jacvallolohi_h;
   double ***jacvalhihilo_h;
   double ***jacvallohilo_h;
   double ***jacvalhilolo_h;
   double ***jacvallololo_h;
   double ***jacvalhihihi_d;
   double ***jacvallohihi_d;
   double ***jacvalhilohi_d;
   double ***jacvallolohi_d;
   double ***jacvalhihilo_d;
   double ***jacvallohilo_d;
   double ***jacvalhilolo_d;
   double ***jacvallololo_d;

   if((mode == 1) || (mode == 2))
   {
      inputhihihi_h = new double*[dim];
      inputlohihi_h = new double*[dim];
      inputhilohi_h = new double*[dim];
      inputlolohi_h = new double*[dim];
      inputhihilo_h = new double*[dim];
      inputlohilo_h = new double*[dim];
      inputhilolo_h = new double*[dim];
      inputlololo_h = new double*[dim];
      outputhihihi_h = new double**[dim];
      outputlohihi_h = new double**[dim];
      outputhilohi_h = new double**[dim];
      outputlolohi_h = new double**[dim];
      outputhihilo_h = new double**[dim];
      outputlohilo_h = new double**[dim];
      outputhilolo_h = new double**[dim];
      outputlololo_h = new double**[dim];
      funvalhihihi_h = new double*[dim];
      funvallohihi_h = new double*[dim];
      funvalhilohi_h = new double*[dim];
      funvallolohi_h = new double*[dim];
      funvalhihilo_h = new double*[dim];
      funvallohilo_h = new double*[dim];
      funvalhilolo_h = new double*[dim];
      funvallololo_h = new double*[dim];
      jacvalhihihi_h = new double**[degp1];
      jacvallohihi_h = new double**[degp1];
      jacvalhilohi_h = new double**[degp1];
      jacvallolohi_h = new double**[degp1];
      jacvalhihilo_h = new double**[degp1];
      jacvallohilo_h = new double**[degp1];
      jacvalhilolo_h = new double**[degp1];
      jacvallololo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputhihihi_d = new double*[dim];
      inputlohihi_d = new double*[dim];
      inputhilohi_d = new double*[dim];
      inputlolohi_d = new double*[dim];
      inputhihilo_d = new double*[dim];
      inputlohilo_d = new double*[dim];
      inputhilolo_d = new double*[dim];
      inputlololo_d = new double*[dim];
      outputhihihi_d = new double**[dim];
      outputlohihi_d = new double**[dim];
      outputhilohi_d = new double**[dim];
      outputlolohi_d = new double**[dim];
      outputhihilo_d = new double**[dim];
      outputlohilo_d = new double**[dim];
      outputhilolo_d = new double**[dim];
      outputlololo_d = new double**[dim];
      funvalhihihi_d = new double*[dim];
      funvallohihi_d = new double*[dim];
      funvalhilohi_d = new double*[dim];
      funvallolohi_d = new double*[dim];
      funvalhihilo_d = new double*[dim];
      funvallohilo_d = new double*[dim];
      funvalhilolo_d = new double*[dim];
      funvallololo_d = new double*[dim];
      jacvalhihihi_d = new double**[degp1];
      jacvallohihi_d = new double**[degp1];
      jacvalhilohi_d = new double**[degp1];
      jacvallolohi_d = new double**[degp1];
      jacvalhihilo_d = new double**[degp1];
      jacvallohilo_d = new double**[degp1];
      jacvalhilolo_d = new double**[degp1];
      jacvallololo_d = new double**[degp1];
   }
   dbl8_allocate_inoutfunjac(dim,deg,mode,
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
       jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhshihihi_h;
   double **rhslohihi_h;
   double **rhshilohi_h;
   double **rhslolohi_h;
   double **rhshihilo_h;
   double **rhslohilo_h;
   double **rhshilolo_h;
   double **rhslololo_h;
   double **rhshihihi_d;
   double **rhslohihi_d;
   double **rhshilohi_d;
   double **rhslolohi_d;
   double **rhshihilo_d;
   double **rhslohilo_d;
   double **rhshilolo_d;
   double **rhslololo_d;
   double **urhshihihi_h;
   double **urhslohihi_h;
   double **urhshilohi_h;
   double **urhslolohi_h;
   double **urhshihilo_h;
   double **urhslohilo_h;
   double **urhshilolo_h;
   double **urhslololo_h;
   double **urhshihihi_d;
   double **urhslohihi_d;
   double **urhshilohi_d;
   double **urhslolohi_d;
   double **urhshihilo_d;
   double **urhslohilo_d;
   double **urhshilolo_d;
   double **urhslololo_d;
   double **solhihihi_h;
   double **sollohihi_h;
   double **solhilohi_h;
   double **sollolohi_h;
   double **solhihilo_h;
   double **sollohilo_h;
   double **solhilolo_h;
   double **sollololo_h;
   double **solhihihi_d;
   double **sollohihi_d;
   double **solhilohi_d;
   double **sollolohi_d;
   double **solhihilo_d;
   double **sollohilo_d;
   double **solhilolo_d;
   double **sollololo_d;
   double **Qhihihi_h;
   double **Qlohihi_h;
   double **Qhilohi_h;
   double **Qlolohi_h;
   double **Qhihilo_h;
   double **Qlohilo_h;
   double **Qhilolo_h;
   double **Qlololo_h;
   double **Qhihihi_d;
   double **Qlohihi_d;
   double **Qhilohi_d;
   double **Qlolohi_d;
   double **Qhihilo_d;
   double **Qlohilo_d;
   double **Qhilolo_d;
   double **Qlololo_d;
   double **Rhihihi_h;
   double **Rlohihi_h;
   double **Rhilohi_h;
   double **Rlolohi_h;
   double **Rhihilo_h;
   double **Rlohilo_h;
   double **Rhilolo_h;
   double **Rlololo_h;
   double **Rhihihi_d;
   double **Rlohihi_d;
   double **Rhilohi_d;
   double **Rlolohi_d;
   double **Rhihilo_d;
   double **Rlohilo_d;
   double **Rhilolo_d;
   double **Rlololo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshihihi_h = new double*[degp1];
      rhslohihi_h = new double*[degp1];
      rhshilohi_h = new double*[degp1];
      rhslolohi_h = new double*[degp1];
      rhshihilo_h = new double*[degp1];
      rhslohilo_h = new double*[degp1];
      rhshilolo_h = new double*[degp1];
      rhslololo_h = new double*[degp1];
      urhshihihi_h = new double*[degp1];
      urhslohihi_h = new double*[degp1];
      urhshilohi_h = new double*[degp1];
      urhslolohi_h = new double*[degp1];
      urhshihilo_h = new double*[degp1];
      urhslohilo_h = new double*[degp1];
      urhshilolo_h = new double*[degp1];
      urhslololo_h = new double*[degp1];
      solhihihi_h = new double*[degp1];
      sollohihi_h = new double*[degp1];
      solhilohi_h = new double*[degp1];
      sollolohi_h = new double*[degp1];
      solhihilo_h = new double*[degp1];
      sollohilo_h = new double*[degp1];
      solhilolo_h = new double*[degp1];
      sollololo_h = new double*[degp1];
      Qhihihi_h = new double*[dim];
      Qlohihi_h = new double*[dim];
      Qhilohi_h = new double*[dim];
      Qlolohi_h = new double*[dim];
      Qhihilo_h = new double*[dim];
      Qlohilo_h = new double*[dim];
      Qhilolo_h = new double*[dim];
      Qlololo_h = new double*[dim];
      Rhihihi_h = new double*[dim];
      Rlohihi_h = new double*[dim];
      Rhilohi_h = new double*[dim];
      Rlolohi_h = new double*[dim];
      Rhihilo_h = new double*[dim];
      Rlohilo_h = new double*[dim];
      Rhilolo_h = new double*[dim];
      Rlololo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhshihihi_d = new double*[degp1];
      rhslohihi_d = new double*[degp1];
      rhshilohi_d = new double*[degp1];
      rhslolohi_d = new double*[degp1];
      rhshihilo_d = new double*[degp1];
      rhslohilo_d = new double*[degp1];
      rhshilolo_d = new double*[degp1];
      rhslololo_d = new double*[degp1];
      urhshihihi_d = new double*[degp1];
      urhslohihi_d = new double*[degp1];
      urhshilohi_d = new double*[degp1];
      urhslolohi_d = new double*[degp1];
      urhshihilo_d = new double*[degp1];
      urhslohilo_d = new double*[degp1];
      urhshilolo_d = new double*[degp1];
      urhslololo_d = new double*[degp1];
      solhihihi_d = new double*[degp1];
      sollohihi_d = new double*[degp1];
      solhilohi_d = new double*[degp1];
      sollolohi_d = new double*[degp1];
      solhihilo_d = new double*[degp1];
      sollohilo_d = new double*[degp1];
      solhilolo_d = new double*[degp1];
      sollololo_d = new double*[degp1];
      Qhihihi_d = new double*[dim];
      Qlohihi_d = new double*[dim];
      Qhilohi_d = new double*[dim];
      Qlolohi_d = new double*[dim];
      Qhihilo_d = new double*[dim];
      Qlohilo_d = new double*[dim];
      Qhilolo_d = new double*[dim];
      Qlololo_d = new double*[dim];
      Rhihihi_d = new double*[dim];
      Rlohihi_d = new double*[dim];
      Rhilohi_d = new double*[dim];
      Rlolohi_d = new double*[dim];
      Rhihilo_d = new double*[dim];
      Rlohilo_d = new double*[dim];
      Rhilolo_d = new double*[dim];
      Rlololo_d = new double*[dim];
   }
   dbl8_allocate_rhsqrsol(dim,deg,mode,
      rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
      rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
      rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
      rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
      urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
      urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
      urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
      urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
      Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
      Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
      Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
      Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
      Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
      Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
      Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
      Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
      solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
      solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
      solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
      solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d);

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
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolhihihi = new double*[dim];
   double **testsollohihi = new double*[dim];
   double **testsolhilohi = new double*[dim];
   double **testsollolohi = new double*[dim];
   double **testsolhihilo = new double*[dim];
   double **testsollohilo = new double*[dim];
   double **testsolhilolo = new double*[dim];
   double **testsollololo = new double*[dim];
   double **mbrhshihihi = new double*[dim];
   double **mbrhslohihi = new double*[dim];
   double **mbrhshilohi = new double*[dim];
   double **mbrhslolohi = new double*[dim];
   double **mbrhshihilo = new double*[dim];
   double **mbrhslohilo = new double*[dim];
   double **mbrhshilolo = new double*[dim];
   double **mbrhslololo = new double*[dim];

   dbl8_column_setup
      (dim,deg,nbrcol,nvr,idx,rowsA,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
       testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
       mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
       mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,mode,vrblvl);

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

      dbl8_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,
          &tailidx_h,&tailidx_d,nvr,idx,exp,nbrfac,expfac,
          mbrhshihihi,mbrhslohihi,mbrhshilohi,mbrhslolohi,
          mbrhshihilo,mbrhslohilo,mbrhshilolo,mbrhslololo,dpr,
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
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
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

   dbl8_error_testsol
      (dim,deg,mode,
       testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
       testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_dbl8_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   double **csthihihi = new double*[dim];
   double **cstlohihi = new double*[dim];
   double **csthilohi = new double*[dim];
   double **cstlolohi = new double*[dim];
   double **csthihilo = new double*[dim];
   double **cstlohilo = new double*[dim];
   double **csthilolo = new double*[dim];
   double **cstlololo = new double*[dim];
   double ***cffhihihi = new double**[dim];
   double ***cfflohihi = new double**[dim];
   double ***cffhilohi = new double**[dim];
   double ***cfflolohi = new double**[dim];
   double ***cffhihilo = new double**[dim];
   double ***cfflohilo = new double**[dim];
   double ***cffhilolo = new double**[dim];
   double ***cfflololo = new double**[dim];

   dbl8_make_coefficients(dim,deg,nbr,nvr,idx,
       csthihihi,cstlohihi,csthilohi,cstlolohi,
       csthihilo,cstlohilo,csthilolo,cstlololo,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,vrblvl);
   // must randomize the leading coefficients,
   // because for exponenials all leading coefficients are one
   for(int i=0; i<dim; i++)
      for(int j=0; j<nbr[i]; j++)
         random_octo_double(&cffhihihi[i][j][0],&cfflohihi[i][j][0],
                            &cffhilohi[i][j][0],&cfflolohi[i][j][0],
                            &cffhihilo[i][j][0],&cfflohilo[i][j][0],
                            &cffhilolo[i][j][0],&cfflololo[i][j][0]);

   double **inputhihihi_h;
   double **inputlohihi_h;
   double **inputhilohi_h;
   double **inputlolohi_h;
   double **inputhihilo_h;
   double **inputlohilo_h;
   double **inputhilolo_h;
   double **inputlololo_h;
   double **inputhihihi_d;
   double **inputlohihi_d;
   double **inputhilohi_d;
   double **inputlolohi_d;
   double **inputhihilo_d;
   double **inputlohilo_d;
   double **inputhilolo_d;
   double **inputlololo_d;
   double ***outputhihihi_h;
   double ***outputlohihi_h;
   double ***outputhilohi_h;
   double ***outputlolohi_h;
   double ***outputhihilo_h;
   double ***outputlohilo_h;
   double ***outputhilolo_h;
   double ***outputlololo_h;
   double ***outputhihihi_d;
   double ***outputlohihi_d;
   double ***outputhilohi_d;
   double ***outputlolohi_d;
   double ***outputhihilo_d;
   double ***outputlohilo_d;
   double ***outputhilolo_d;
   double ***outputlololo_d;
   double **funvalhihihi_h;
   double **funvallohihi_h;
   double **funvalhilohi_h;
   double **funvallolohi_h;
   double **funvalhihilo_h;
   double **funvallohilo_h;
   double **funvalhilolo_h;
   double **funvallololo_h;
   double **funvalhihihi_d;
   double **funvallohihi_d;
   double **funvalhilohi_d;
   double **funvallolohi_d;
   double **funvalhihilo_d;
   double **funvallohilo_d;
   double **funvalhilolo_d;
   double **funvallololo_d;
   double ***jacvalhihihi_h;
   double ***jacvallohihi_h;
   double ***jacvalhilohi_h;
   double ***jacvallolohi_h;
   double ***jacvalhihilo_h;
   double ***jacvallohilo_h;
   double ***jacvalhilolo_h;
   double ***jacvallololo_h;
   double ***jacvalhihihi_d;
   double ***jacvallohihi_d;
   double ***jacvalhilohi_d;
   double ***jacvallolohi_d;
   double ***jacvalhihilo_d;
   double ***jacvallohilo_d;
   double ***jacvalhilolo_d;
   double ***jacvallololo_d;

   if((mode == 1) || (mode == 2))
   {
      inputhihihi_h = new double*[dim];
      inputlohihi_h = new double*[dim];
      inputhilohi_h = new double*[dim];
      inputlolohi_h = new double*[dim];
      inputhihilo_h = new double*[dim];
      inputlohilo_h = new double*[dim];
      inputhilolo_h = new double*[dim];
      inputlololo_h = new double*[dim];
      outputhihihi_h = new double**[dim];
      outputlohihi_h = new double**[dim];
      outputhilohi_h = new double**[dim];
      outputlolohi_h = new double**[dim];
      outputhihilo_h = new double**[dim];
      outputlohilo_h = new double**[dim];
      outputhilolo_h = new double**[dim];
      outputlololo_h = new double**[dim];
      funvalhihihi_h = new double*[dim];
      funvallohihi_h = new double*[dim];
      funvalhilohi_h = new double*[dim];
      funvallolohi_h = new double*[dim];
      funvalhihilo_h = new double*[dim];
      funvallohilo_h = new double*[dim];
      funvalhilolo_h = new double*[dim];
      funvallololo_h = new double*[dim];
      jacvalhihihi_h = new double**[degp1];
      jacvallohihi_h = new double**[degp1];
      jacvalhilohi_h = new double**[degp1];
      jacvallolohi_h = new double**[degp1];
      jacvalhihilo_h = new double**[degp1];
      jacvallohilo_h = new double**[degp1];
      jacvalhilolo_h = new double**[degp1];
      jacvallololo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputhihihi_d = new double*[dim];
      inputlohihi_d = new double*[dim];
      inputhilohi_d = new double*[dim];
      inputlolohi_d = new double*[dim];
      inputhihilo_d = new double*[dim];
      inputlohilo_d = new double*[dim];
      inputhilolo_d = new double*[dim];
      inputlololo_d = new double*[dim];
      outputhihihi_d = new double**[dim];
      outputlohihi_d = new double**[dim];
      outputhilohi_d = new double**[dim];
      outputlolohi_d = new double**[dim];
      outputhihilo_d = new double**[dim];
      outputlohilo_d = new double**[dim];
      outputhilolo_d = new double**[dim];
      outputlololo_d = new double**[dim];
      funvalhihihi_d = new double*[dim];
      funvallohihi_d = new double*[dim];
      funvalhilohi_d = new double*[dim];
      funvallolohi_d = new double*[dim];
      funvalhihilo_d = new double*[dim];
      funvallohilo_d = new double*[dim];
      funvalhilolo_d = new double*[dim];
      funvallololo_d = new double*[dim];
      jacvalhihihi_d = new double**[degp1];
      jacvallohihi_d = new double**[degp1];
      jacvalhilohi_d = new double**[degp1];
      jacvallolohi_d = new double**[degp1];
      jacvalhihilo_d = new double**[degp1];
      jacvallohilo_d = new double**[degp1];
      jacvalhilolo_d = new double**[degp1];
      jacvallololo_d = new double**[degp1];
   }
   dbl8_allocate_inoutfunjac(dim,deg,mode,
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
       jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhshihihi_h;
   double **rhslohihi_h;
   double **rhshilohi_h;
   double **rhslolohi_h;
   double **rhshihilo_h;
   double **rhslohilo_h;
   double **rhshilolo_h;
   double **rhslololo_h;
   double **rhshihihi_d;
   double **rhslohihi_d;
   double **rhshilohi_d;
   double **rhslolohi_d;
   double **rhshihilo_d;
   double **rhslohilo_d;
   double **rhshilolo_d;
   double **rhslololo_d;
   double **urhshihihi_h;
   double **urhslohihi_h;
   double **urhshilohi_h;
   double **urhslolohi_h;
   double **urhshihilo_h;
   double **urhslohilo_h;
   double **urhshilolo_h;
   double **urhslololo_h;
   double **urhshihihi_d;
   double **urhslohihi_d;
   double **urhshilohi_d;
   double **urhslolohi_d;
   double **urhshihilo_d;
   double **urhslohilo_d;
   double **urhshilolo_d;
   double **urhslololo_d;
   double **solhihihi_h;
   double **sollohihi_h;
   double **solhilohi_h;
   double **sollolohi_h;
   double **solhihilo_h;
   double **sollohilo_h;
   double **solhilolo_h;
   double **sollololo_h;
   double **solhihihi_d;
   double **sollohihi_d;
   double **solhilohi_d;
   double **sollolohi_d;
   double **solhihilo_d;
   double **sollohilo_d;
   double **solhilolo_d;
   double **sollololo_d;
   double **Qhihihi_h;
   double **Qlohihi_h;
   double **Qhilohi_h;
   double **Qlolohi_h;
   double **Qhihilo_h;
   double **Qlohilo_h;
   double **Qhilolo_h;
   double **Qlololo_h;
   double **Qhihihi_d;
   double **Qlohihi_d;
   double **Qhilohi_d;
   double **Qlolohi_d;
   double **Qhihilo_d;
   double **Qlohilo_d;
   double **Qhilolo_d;
   double **Qlololo_d;
   double **Rhihihi_h;
   double **Rlohihi_h;
   double **Rhilohi_h;
   double **Rlolohi_h;
   double **Rhihilo_h;
   double **Rlohilo_h;
   double **Rhilolo_h;
   double **Rlololo_h;
   double **Rhihihi_d;
   double **Rlohihi_d;
   double **Rhilohi_d;
   double **Rlolohi_d;
   double **Rhihilo_d;
   double **Rlohilo_d;
   double **Rhilolo_d;
   double **Rlololo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshihihi_h = new double*[degp1];
      rhslohihi_h = new double*[degp1];
      rhshilohi_h = new double*[degp1];
      rhslolohi_h = new double*[degp1];
      rhshihilo_h = new double*[degp1];
      rhslohilo_h = new double*[degp1];
      rhshilolo_h = new double*[degp1];
      rhslololo_h = new double*[degp1];
      urhshihihi_h = new double*[degp1];
      urhslohihi_h = new double*[degp1];
      urhshilohi_h = new double*[degp1];
      urhslolohi_h = new double*[degp1];
      urhshihilo_h = new double*[degp1];
      urhslohilo_h = new double*[degp1];
      urhshilolo_h = new double*[degp1];
      urhslololo_h = new double*[degp1];
      solhihihi_h = new double*[degp1];
      sollohihi_h = new double*[degp1];
      solhilohi_h = new double*[degp1];
      sollolohi_h = new double*[degp1];
      solhihilo_h = new double*[degp1];
      sollohilo_h = new double*[degp1];
      solhilolo_h = new double*[degp1];
      sollololo_h = new double*[degp1];
      Qhihihi_h = new double*[dim];
      Qlohihi_h = new double*[dim];
      Qhilohi_h = new double*[dim];
      Qlolohi_h = new double*[dim];
      Qhihilo_h = new double*[dim];
      Qlohilo_h = new double*[dim];
      Qhilolo_h = new double*[dim];
      Qlololo_h = new double*[dim];
      Rhihihi_h = new double*[dim];
      Rlohihi_h = new double*[dim];
      Rhilohi_h = new double*[dim];
      Rlolohi_h = new double*[dim];
      Rhihilo_h = new double*[dim];
      Rlohilo_h = new double*[dim];
      Rhilolo_h = new double*[dim];
      Rlololo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhshihihi_d = new double*[degp1];
      rhslohihi_d = new double*[degp1];
      rhshilohi_d = new double*[degp1];
      rhslolohi_d = new double*[degp1];
      rhshihilo_d = new double*[degp1];
      rhslohilo_d = new double*[degp1];
      rhshilolo_d = new double*[degp1];
      rhslololo_d = new double*[degp1];
      urhshihihi_d = new double*[degp1];
      urhslohihi_d = new double*[degp1];
      urhshilohi_d = new double*[degp1];
      urhslolohi_d = new double*[degp1];
      urhshihilo_d = new double*[degp1];
      urhslohilo_d = new double*[degp1];
      urhshilolo_d = new double*[degp1];
      urhslololo_d = new double*[degp1];
      solhihihi_d = new double*[degp1];
      sollohihi_d = new double*[degp1];
      solhilohi_d = new double*[degp1];
      sollolohi_d = new double*[degp1];
      solhihilo_d = new double*[degp1];
      sollohilo_d = new double*[degp1];
      solhilolo_d = new double*[degp1];
      sollololo_d = new double*[degp1];
      Qhihihi_d = new double*[dim];
      Qlohihi_d = new double*[dim];
      Qhilohi_d = new double*[dim];
      Qlolohi_d = new double*[dim];
      Qhihilo_d = new double*[dim];
      Qlohilo_d = new double*[dim];
      Qhilolo_d = new double*[dim];
      Qlololo_d = new double*[dim];
      Rhihihi_d = new double*[dim];
      Rlohihi_d = new double*[dim];
      Rhilohi_d = new double*[dim];
      Rlolohi_d = new double*[dim];
      Rhihilo_d = new double*[dim];
      Rlohilo_d = new double*[dim];
      Rhilolo_d = new double*[dim];
      Rlololo_d = new double*[dim];
   }
   dbl8_allocate_rhsqrsol(dim,deg,mode,
      rhshihihi_h,rhslohihi_h,rhshilohi_h,rhslolohi_h,
      rhshihilo_h,rhslohilo_h,rhshilolo_h,rhslololo_h,
      rhshihihi_d,rhslohihi_d,rhshilohi_d,rhslolohi_d,
      rhshihilo_d,rhslohilo_d,rhshilolo_d,rhslololo_d,
      urhshihihi_h,urhslohihi_h,urhshilohi_h,urhslolohi_h,
      urhshihilo_h,urhslohilo_h,urhshilolo_h,urhslololo_h,
      urhshihihi_d,urhslohihi_d,urhshilohi_d,urhslolohi_d,
      urhshihilo_d,urhslohilo_d,urhshilolo_d,urhslololo_d,
      Qhihihi_h,Qlohihi_h,Qhilohi_h,Qlolohi_h,
      Qhihilo_h,Qlohilo_h,Qhilolo_h,Qlololo_h,
      Qhihihi_d,Qlohihi_d,Qhilohi_d,Qlolohi_d,
      Qhihilo_d,Qlohilo_d,Qhilolo_d,Qlololo_d,
      Rhihihi_h,Rlohihi_h,Rhilohi_h,Rlolohi_h,
      Rhihilo_h,Rlohilo_h,Rhilolo_h,Rlololo_h,
      Rhihihi_d,Rlohihi_d,Rhilohi_d,Rlolohi_d,
      Rhihilo_d,Rlohilo_d,Rhilolo_d,Rlololo_d,
      solhihihi_h,sollohihi_h,solhilohi_h,sollolohi_h,
      solhihilo_h,sollohilo_h,solhilolo_h,sollololo_h,
      solhihihi_d,sollohihi_d,solhilohi_d,sollolohi_d,
      solhihilo_d,sollohilo_d,solhilolo_d,sollololo_d);
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolhihihi = new double*[dim];
   double **testsollohihi = new double*[dim];
   double **testsolhilohi = new double*[dim];
   double **testsollolohi = new double*[dim];
   double **testsolhihilo = new double*[dim];
   double **testsollohilo = new double*[dim];
   double **testsolhilolo = new double*[dim];
   double **testsollololo = new double*[dim];
   double **mbrhshihihi = new double*[dim];
   double **mbrhslohihi = new double*[dim];
   double **mbrhshilohi = new double*[dim];
   double **mbrhslolohi = new double*[dim];
   double **mbrhshihilo = new double*[dim];
   double **mbrhslohilo = new double*[dim];
   double **mbrhshilolo = new double*[dim];
   double **mbrhslololo = new double*[dim];

   dbl8_row_setup
      (dim,deg,nbr,nvr,idx,
       csthihihi,cstlohihi,csthilohi,cstlolohi,
       csthihilo,cstlohilo,csthilolo,cstlololo,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
       testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d,
       outputhihihi_h,outputlohihi_h,outputhilohi_h,outputlolohi_h,
       outputhihilo_h,outputlohilo_h,outputhilolo_h,outputlololo_h,
       outputhihihi_d,outputlohihi_d,outputhilohi_d,outputlolohi_d,
       outputhihilo_d,outputlohilo_d,outputhilolo_d,outputlololo_d,
       mode,vrblvl);

   // allocating extra work space

   double *workvechihihi = new double[dim];
   double *workveclohihi = new double[dim];
   double *workvechilohi = new double[dim];
   double *workveclolohi = new double[dim];
   double *workvechihilo = new double[dim];
   double *workveclohilo = new double[dim];
   double *workvechilolo = new double[dim];
   double *workveclololo = new double[dim];
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

      dbl8_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,nbr,nvr,idx,
          csthihihi,cstlohihi,csthilohi,cstlolohi,
          csthihilo,cstlohilo,csthilolo,cstlololo,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,dpr,
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
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
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

   dbl8_error_testsol
      (dim,deg,mode,
       testsolhihihi,testsollohihi,testsolhilohi,testsollolohi,
       testsolhihilo,testsollohilo,testsolhilolo,testsollololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       inputhihihi_d,inputlohihi_d,inputhilohi_d,inputlolohi_d,
       inputhihilo_d,inputlohilo_d,inputhilolo_d,inputlololo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
