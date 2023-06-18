// The file cmplx4_newton_method.cpp defines the functions with prototypes in
// the file cmplx4_newton_method.h.

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
#include "random_monomials.h"
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

int cmplx4_errors_funjacrhs
 ( int dim, int deg,
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
   double **rhsimhilo_d, double **rhsimlolo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-50;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU function values ... " << endl;
   errsum = cmplx4_error2sum(dim,degp1,
               funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
               funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
               funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
               funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
               "funval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU Jacobians ... " << endl;
   errsum = cmplx4_error3sum(degp1,dim,dim,
               jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
               jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
               jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
               jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
               "jacval",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU right hand sides ... " << endl;
   errsum = cmplx4_error2sum(degp1,dim,
               rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
               rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
               rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
               rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
               "rhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int cmplx4_errors_inurhsQRsol
 ( int dim, int deg,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d,
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
   double **solimhilo_d, double **solimlolo_d, int vrblvl )
{
   const int degp1 = deg+1;
   const double tol = 1.0e-50;
   int fail = 0;
   double errsum = 0.0;

   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices Q ... " << endl;
   errsum = cmplx4_error2sum(dim,dim,
               Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
               Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
               Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
               Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,"Q",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU matrices R ... " << endl;
   errsum = cmplx4_error2sum(dim,dim,
               Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
               Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
               Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
               Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,"R",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU updated rhs ... " << endl;
   errsum = cmplx4_error2sum(degp1,dim,
               urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
               urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
               urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
               urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
               "urhs",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU update to solutions ... " << endl;
   errsum = cmplx4_error2sum(degp1,dim,
               solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
               solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
               solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
               solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
               "sol",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);
   cout << scientific << setprecision(16);
   cout << "comparing CPU with GPU series ... " << endl;
   errsum = cmplx4_error2sum(dim,degp1,
               inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
               inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
               inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
               inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
               "input",vrblvl);
   cout << scientific << setprecision(2)
        << "sum of errors : " << errsum << endl;
   fail += (errsum > tol);

   return fail;
}

int cmplx4_update_newton_qrstep
 ( int szt, int nbt, int dim, int deg,
   int *tailidx_h, int *tailidx_d,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d,
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
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
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
            urhsrehihi_h[i][j] = rhsrehihi_h[i][j];
            urhsrehilo_h[i][j] = rhsrehilo_h[i][j];
            urhsimhihi_h[i][j] = rhsimhihi_h[i][j];
            urhsimhilo_h[i][j] = rhsimhilo_h[i][j];
            urhsrelohi_h[i][j] = rhsrelohi_h[i][j];
            urhsrelolo_h[i][j] = rhsrelolo_h[i][j];
            urhsimlohi_h[i][j] = rhsimlohi_h[i][j];
            urhsimlolo_h[i][j] = rhsimlolo_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solrehihi_h[i][j] = 0.0; solimhihi_h[i][j] = 0.0;
            solrelohi_h[i][j] = 0.0; solimlohi_h[i][j] = 0.0;
            solrehilo_h[i][j] = 0.0; solimhilo_h[i][j] = 0.0;
            solrelolo_h[i][j] = 0.0; solimlolo_h[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling CPU_cmplx4_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_cmplx4_qrbs_solve
         (dim,degp1,oldtail,
          jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
          jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
          urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
          urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
          solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
          solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
          Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
          Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
          Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
          Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
          workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
          workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,
          zeroQ_h,noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_cmplx4_linear_residue ..." << endl;

         CPU_cmplx4_linear_residue
            (dim,degp1,*tailidx_h-1,
             jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
             jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
             rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
             rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
             solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
             solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
             resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
             resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
             resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhihi << endl;
      }
      cmplx4_update_series
         (dim,degp1,*tailidx_h-1,
          inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
          inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
          solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
          solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
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
            solrehihi_d[i][j] = 0.0; solimhihi_d[i][j] = 0.0;
            solrelohi_d[i][j] = 0.0; solimlohi_d[i][j] = 0.0;
            solrehilo_d[i][j] = 0.0; solimhilo_d[i][j] = 0.0;
            solrelolo_d[i][j] = 0.0; solimlolo_d[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling GPU_cmplx4_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_cmplx4_bals_solve
         (dim,degp1,szt,nbt,oldtail,
          jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
          jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
          Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
          Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
          Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
          Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
          urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
          urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
          solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
          solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,zeroQ_d,noqr_d,
          upidx_d,bsidx_d,&newtail,totqrlapsedms,totqtblapsedms,
          totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_cmplx4_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_cmplx4_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,
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

         cout << scientific << setprecision(3) 
              << "maximum residual : " << *resmaxhihi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;

         *totreslapsedms += elapsedms;
      }
      cmplx4_update_series
         (dim,degp1,*tailidx_d-1,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
          solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
          solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,vrblvl);
   }
   return 0;
}

int cmplx4_column_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehihi, double **mbrelohi, double **mbrehilo, double **mbrelolo,
   double **mbimhihi, double **mbimlohi, double **mbimhilo, double **mbimlolo,
   double dpr,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **accrehihi, double **accrelohi,
   double **accrehilo, double **accrelolo,
   double **accimhihi, double **accimlohi,
   double **accimhilo, double **accimlolo,
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
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
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
         cout << "calling CPU_cmplx4_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx4_unit_series_vector
            (dim,deg,cffrehihi[0],cffrelohi[0],cffrehilo[0],cffrelolo[0],
                     cffimhihi[0],cffimlohi[0],cffimhilo[0],cffimlolo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_cmplx4_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehihi[0],cffrelohi[0],cffrehilo[0],cffrelolo[0],
             cffimhihi[0],cffimlohi[0],cffimhilo[0],cffimlolo[0],
             accrehihi[0],accrelohi[0],accrehilo[0],accrelolo[0],
             accimhihi[0],accimlohi[0],accimhilo[0],accimlolo[0],
             inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
             inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
             outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
             outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
             vrblvl);
      }
      else
         CPU_cmplx4_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo,
             accrehihi,accrelohi,accrehilo,accrelolo,
             accimhihi,accimlohi,accimhilo,accimlolo,
             inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
             inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
             funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
             funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
             jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
             jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
             vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_cmplx4_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx4_unit_series_vector
            (dim,deg,cffrehihi[0],cffrelohi[0],cffrehilo[0],cffrelolo[0],
                     cffimhihi[0],cffimlohi[0],cffimhilo[0],cffimlolo[0]);
         // reset the coefficients
/*
         GPU_cmplx4_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehihi[0],cffrelohi[0],cffrehilo[0],cffrelolo[0],
             cffimhihi[0],cffimlohi[0],cffimhilo[0],cffimlolo[0],
             accrehihi[0],accrelohi[0],accrehilo[0],accrelolo[0],
             accimhihi[0],accimlohi[0],accimhilo[0],accimlolo[0],
             inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
             inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
             outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
             outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
             totcnvlapsedms,vrblvl);
 */
         GPU_cmplx4vectorized_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehihi[0],cffrelohi[0],cffrehilo[0],cffrelolo[0],
             cffimhihi[0],cffimlohi[0],cffimhilo[0],cffimlolo[0],
             accrehihi[0],accrelohi[0],accrehilo[0],accrelolo[0],
             accimhihi[0],accimlohi[0],accimhilo[0],accimlolo[0],
             inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
             inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
             outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
             outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
             totcnvlapsedms,vrblvl);
      }
      else
         GPU_cmplx4_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,
             cffrehihi,cffrelohi,cffrehilo,cffrelolo,
             cffimhihi,cffimlohi,cffimhilo,cffimlolo,
             inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
             inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
             outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
             outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
             funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
             funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
             jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
             jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
             totcnvlapsedms,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
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

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum << endl;
   }
   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         cmplx4_define_rhs
            (dim,degp1,
             mbrehihi,mbrelohi,mbrehilo,mbrelolo,
             mbimhihi,mbimlohi,mbimhilo,mbimlolo,
             funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
             funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
             rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
             rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalrehihi_h[i][j][k] = 0.0;
                  jacvalrelohi_h[i][j][k] = 0.0;
                  jacvalrehilo_h[i][j][k] = 0.0;
                  jacvalrelolo_h[i][j][k] = 0.0;
                  jacvalimhihi_h[i][j][k] = 0.0;
                  jacvalimlohi_h[i][j][k] = 0.0;
                  jacvalimhilo_h[i][j][k] = 0.0;
                  jacvalimlolo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         cmplx4_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],
             mbrehihi,mbrelohi,mbrehilo,mbrelolo,
             mbimhihi,mbimlohi,mbimhilo,mbimlolo,dpr,
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
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         cmplx4_define_rhs
            (dim,degp1,
             mbrehihi,mbrelohi,mbrehilo,mbrelolo,
             mbimhihi,mbimlohi,mbimhilo,mbimlolo,
             funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
             funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
             rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
             rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalrehihi_d[i][j][k] = 0.0; jacvalimhihi_d[i][j][k] = 0.0;
                  jacvalrelohi_d[i][j][k] = 0.0; jacvalimlohi_d[i][j][k] = 0.0;
                  jacvalrehilo_d[i][j][k] = 0.0; jacvalimhilo_d[i][j][k] = 0.0;
                  jacvalrelolo_d[i][j][k] = 0.0; jacvalimlolo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         cmplx4_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],
             mbrehihi,mbrelohi,mbrehilo,mbrelolo,
             mbimhihi,mbimlohi,mbimhilo,mbimlolo,dpr,
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
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = cmplx4_errors_funjacrhs(dim,deg,
                    funvalrehihi_h,funvalrelohi_h,
                    funvalrelohi_h,funvalrelolo_h,
                    funvalimhihi_h,funvalimlohi_h,
                    funvalimlohi_h,funvalimlolo_h,
                    funvalrehihi_d,funvalrelohi_d,
                    funvalrelohi_d,funvalrelolo_d,
                    funvalimhihi_d,funvalimlohi_d,
                    funvalimlohi_d,funvalimlolo_d,
                    jacvalrehihi_h,jacvalrelohi_h,
                    jacvalrelohi_h,jacvalrelolo_h,
                    jacvalimhihi_h,jacvalimlohi_h,
                    jacvalimlohi_h,jacvalimlolo_h,
                    jacvalrehihi_d,jacvalrelohi_d,
                    jacvalrelohi_d,jacvalrelolo_d,
                    jacvalimhihi_d,jacvalimlohi_d,
                    jacvalimlohi_d,jacvalimlolo_d,
                    rhsrehihi_h,rhsrelohi_h,rhsrelohi_h,rhsrelolo_h,
                    rhsimhihi_h,rhsimlohi_h,rhsimlohi_h,rhsimlolo_h,
                    rhsrehihi_d,rhsrelohi_d,rhsrelohi_d,rhsrelolo_d,
                    rhsimhihi_d,rhsimlohi_d,rhsimlohi_d,rhsimlolo_d,vrblvl);
   }
   cmplx4_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
       inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
       inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
       inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
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
       workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
       workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,
       resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
       resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
       resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return cmplx4_errors_inurhsQRsol(dim,deg,
                inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
                inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
                inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
                inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
                Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
                Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
                Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
                Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
                Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
                Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
                Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
                Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
                urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
                urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
                urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
                urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
                solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
                solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
                solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
                solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
                vrblvl);
   }
   else 
      return 0;
}

int cmplx4_row_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int *tailidx_h, int *tailidx_d,
   int *nbr, int **nvr, int ***idx,
   double **cstrehihi, double **cstrelohi,
   double **cstrehilo, double **cstrelolo,
   double **cstimhihi, double **cstimlohi,
   double **cstimhilo, double **cstimlolo,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double dpr,
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
   double *workvecrehihi, double *workvecrelohi,
   double *workvecrehilo, double *workvecrelolo,
   double *workvecimhihi, double *workvecimlohi,
   double *workvecimhilo, double *workvecimlolo,
   double **resvecrehihi, double **resvecrelohi,
   double **resvecrehilo, double **resvecrelolo,
   double **resvecimhihi, double **resvecimlohi,
   double **resvecimhilo, double **resvecimlolo,
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
         cout << "calling CPU_cmplx4_poly_evaldiff ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_cmplx4_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],
            cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
            cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
            cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
            cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
            inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
            inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
            outputrehihi_h[i],outputrelohi_h[i],
            outputrehilo_h[i],outputrelolo_h[i],
            outputimhihi_h[i],outputimlohi_h[i],
            outputimhilo_h[i],outputimlolo_h[i],&lapsed,0); // no output

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
         cout << "calling GPU_cmplx4_poly_evaldiff ..." << endl;

      double timelapsed_d = 0.0;

      for(int i=0; i<dim; i++)
      {
         double cnvlapms,addlapms,timelapms_d,walltimes_d;

         ComplexConvolutionJobs cnvjobs(dim);
         ComplexIncrementJobs incjobs(cnvjobs,false);
         ComplexAdditionJobs addjobs(dim,nbr[i]);

         make_all_complex_jobs
            (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,false);

         GPU_cmplx4vectorized_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
             cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
             cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
             cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
             inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
             inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
             outputrehihi_d[i],outputrelohi_d[i],
             outputrehilo_d[i],outputrelolo_d[i],
             outputimhihi_d[i],outputimlohi_d[i],
             outputimhilo_d[i],outputimlolo_d[i],cnvjobs,incjobs,addjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0); // vrblvl);

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

      errsum = cmplx4_error3sum(dim,dim+1,degp1,
                  outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
                  outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
                  outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
                  outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
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
      cmplx4_map_evaldiff_output(dim,deg,
         outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
         outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
         funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
         funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
         jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
         jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsrehihi_h[i][j] = -funvalrehihi_h[j][i];
            rhsrelohi_h[i][j] = -funvalrelohi_h[j][i];
            rhsrehilo_h[i][j] = -funvalrehilo_h[j][i];
            rhsrelolo_h[i][j] = -funvalrelolo_h[j][i];
            rhsimhihi_h[i][j] = -funvalimhihi_h[j][i];
            rhsimlohi_h[i][j] = -funvalimlohi_h[j][i];
            rhsimhilo_h[i][j] = -funvalimhilo_h[j][i];
            rhsimlolo_h[i][j] = -funvalimlolo_h[j][i];
         }
   }
   if((mode == 0) || (mode == 2))
   {
      cmplx4_map_evaldiff_output(dim,deg,
         outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
         outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,
         funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
         funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
         jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
         jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,vrblvl);

      for(int i=0; i<degp1; i++)
         for(int j=0; j<dim; j++)
         {
            rhsrehihi_d[i][j] = -funvalrehihi_d[j][i];
            rhsrelohi_d[i][j] = -funvalrelohi_d[j][i];
            rhsrehilo_d[i][j] = -funvalrehilo_d[j][i];
            rhsrelolo_d[i][j] = -funvalrelolo_d[j][i];
            rhsimhihi_d[i][j] = -funvalimhihi_d[j][i];
            rhsimlohi_d[i][j] = -funvalimlohi_d[j][i];
            rhsimhilo_d[i][j] = -funvalimhilo_d[j][i];
            rhsimlolo_d[i][j] = -funvalimlolo_d[j][i];
         }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      int fail = cmplx4_errors_funjacrhs(dim,deg,
                    funvalrehihi_h,funvalrelohi_h,
                    funvalrelohi_h,funvalrelolo_h,
                    funvalimhihi_h,funvalimlohi_h,
                    funvalimlohi_h,funvalimlolo_h,
                    funvalrehihi_d,funvalrelohi_d,
                    funvalrelohi_d,funvalrelolo_d,
                    funvalimhihi_d,funvalimlohi_d,
                    funvalimlohi_d,funvalimlolo_d,
                    jacvalrehihi_h,jacvalrelohi_h,
                    jacvalrelohi_h,jacvalrelolo_h,
                    jacvalimhihi_h,jacvalimlohi_h,
                    jacvalimlohi_h,jacvalimlolo_h,
                    jacvalrehihi_d,jacvalrelohi_d,
                    jacvalrelohi_d,jacvalrelolo_d,
                    jacvalimhihi_d,jacvalimlohi_d,
                    jacvalimlohi_d,jacvalimlolo_d,
                    rhsrehihi_h,rhsrelohi_h,rhsrelohi_h,rhsrelolo_h,
                    rhsimhihi_h,rhsimlohi_h,rhsimlohi_h,rhsimlolo_h,
                    rhsrehihi_d,rhsrelohi_d,rhsrelohi_d,rhsrelolo_d,
                    rhsimhihi_d,rhsimlohi_d,rhsimlohi_d,rhsimlolo_d,vrblvl);
   }
   cmplx4_update_newton_qrstep
      (szt,nbt,dim,deg,tailidx_h,tailidx_d,
       inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
       inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
       inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
       inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
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
       workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
       workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,
       resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
       resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
       resmaxhihi,resmaxlohi,resmaxhilo,resmaxlolo,
       zeroQ_h,noqr_h,zeroQ_d,noqr_d,upidx_h,bsidx_h,upidx_d,bsidx_d,
       totqrlapsedms,totqtblapsedms,totbslapsedms,
       totupdlapsedms,totreslapsedms,vrblvl,mode);

   if((vrblvl > 0) && (mode == 2))
   {
      return cmplx4_errors_inurhsQRsol(dim,deg,
                inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
                inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
                inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
                inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
                Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
                Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
                Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
                Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
                Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
                Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
                Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
                Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
                urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
                urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
                urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
                urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
                solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
                solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
                solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
                solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d,
                vrblvl);
   }
   else 
      return 0;
}

int cmplx4_allocate_inoutfunjac
 ( int dim, int deg, int mode,
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
   double ***jacvalimhilo_d, double ***jacvalimlolo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
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

         outputrehihi_h[i] = new double*[dim+1];
         outputrelohi_h[i] = new double*[dim+1];
         outputrehilo_h[i] = new double*[dim+1];
         outputrelolo_h[i] = new double*[dim+1];
         outputimhihi_h[i] = new double*[dim+1];
         outputimlohi_h[i] = new double*[dim+1];
         outputimhilo_h[i] = new double*[dim+1];
         outputimlolo_h[i] = new double*[dim+1];

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
         }
         funvalrehihi_h[i] = new double[degp1];
         funvalrelohi_h[i] = new double[degp1];
         funvalrehilo_h[i] = new double[degp1];
         funvalrelolo_h[i] = new double[degp1];
         funvalimhihi_h[i] = new double[degp1];
         funvalimlohi_h[i] = new double[degp1];
         funvalimhilo_h[i] = new double[degp1];
         funvalimlolo_h[i] = new double[degp1];
      }
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
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<dim; i++)
      {
         inputrehihi_d[i] = new double[degp1];
         inputrelohi_d[i] = new double[degp1];
         inputrehilo_d[i] = new double[degp1];
         inputrelolo_d[i] = new double[degp1];
         inputimhihi_d[i] = new double[degp1];
         inputimlohi_d[i] = new double[degp1];
         inputimhilo_d[i] = new double[degp1];
         inputimlolo_d[i] = new double[degp1];

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
            outputrehihi_d[i][j] = new double[degp1];
            outputrelohi_d[i][j] = new double[degp1];
            outputrehilo_d[i][j] = new double[degp1];
            outputrelolo_d[i][j] = new double[degp1];
            outputimhihi_d[i][j] = new double[degp1];
            outputimlohi_d[i][j] = new double[degp1];
            outputimhilo_d[i][j] = new double[degp1];
            outputimlolo_d[i][j] = new double[degp1];
         }
         funvalrehihi_d[i] = new double[degp1];
         funvalrelohi_d[i] = new double[degp1];
         funvalrehilo_d[i] = new double[degp1];
         funvalrelolo_d[i] = new double[degp1];
         funvalimhihi_d[i] = new double[degp1];
         funvalimlohi_d[i] = new double[degp1];
         funvalimhilo_d[i] = new double[degp1];
         funvalimlolo_d[i] = new double[degp1];
      }
      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
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
   }
   return 0;
}

int cmplx4_allocate_rhsqrsol
 ( int dim, int deg, int mode,
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
   double **solrehihi_h, double **solrelohi_h,
   double **solrehilo_h, double **solrelolo_h,
   double **solimhihi_h, double **solimlohi_h,
   double **solimhilo_h, double **solimlolo_h,
   double **solrehihi_d, double **solrelohi_d,
   double **solrehilo_d, double **solrelolo_d,
   double **solimhihi_d, double **solimlohi_d,
   double **solimhilo_d, double **solimlolo_d )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
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
         urhsrehihi_h[i] = new double[dim];
         urhsrelohi_h[i] = new double[dim];
         urhsrehilo_h[i] = new double[dim];
         urhsrelolo_h[i] = new double[dim];
         urhsimhihi_h[i] = new double[dim];
         urhsimlohi_h[i] = new double[dim];
         urhsimhilo_h[i] = new double[dim];
         urhsimlolo_h[i] = new double[dim];
         solrehihi_h[i] = new double[dim];
         solrelohi_h[i] = new double[dim];
         solrehilo_h[i] = new double[dim];
         solrelolo_h[i] = new double[dim];
         solimhihi_h[i] = new double[dim];
         solimlohi_h[i] = new double[dim];
         solimhilo_h[i] = new double[dim];
         solimlolo_h[i] = new double[dim];
      }
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
         Rrehihi_h[i] = new double[dim];
         Rrelohi_h[i] = new double[dim];
         Rrehilo_h[i] = new double[dim];
         Rrelolo_h[i] = new double[dim];
         Rimhihi_h[i] = new double[dim];
         Rimlohi_h[i] = new double[dim];
         Rimhilo_h[i] = new double[dim];
         Rimlolo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) 
      {
         rhsrehihi_d[i] = new double[dim];
         rhsrelohi_d[i] = new double[dim];
         rhsrehilo_d[i] = new double[dim];
         rhsrelolo_d[i] = new double[dim];
         rhsimhihi_d[i] = new double[dim];
         rhsimlohi_d[i] = new double[dim];
         rhsimhilo_d[i] = new double[dim];
         rhsimlolo_d[i] = new double[dim];
         urhsrehihi_d[i] = new double[dim];
         urhsrelohi_d[i] = new double[dim];
         urhsrehilo_d[i] = new double[dim];
         urhsrelolo_d[i] = new double[dim];
         urhsimhihi_d[i] = new double[dim];
         urhsimlohi_d[i] = new double[dim];
         urhsimhilo_d[i] = new double[dim];
         urhsimlolo_d[i] = new double[dim];
         solrehihi_d[i] = new double[dim];
         solrelohi_d[i] = new double[dim];
         solrehilo_d[i] = new double[dim];
         solrelolo_d[i] = new double[dim];
         solimhihi_d[i] = new double[dim];
         solimlohi_d[i] = new double[dim];
         solimhilo_d[i] = new double[dim];
         solimlolo_d[i] = new double[dim];
      }
      for(int i=0; i<dim; i++)
      {
         Qrehihi_d[i] = new double[dim];
         Qrelohi_d[i] = new double[dim];
         Qrehilo_d[i] = new double[dim];
         Qrelolo_d[i] = new double[dim];
         Qimhihi_d[i] = new double[dim];
         Qimlohi_d[i] = new double[dim];
         Qimhilo_d[i] = new double[dim];
         Qimlolo_d[i] = new double[dim];
         Rrehihi_d[i] = new double[dim];
         Rrelohi_d[i] = new double[dim];
         Rrehilo_d[i] = new double[dim];
         Rrelolo_d[i] = new double[dim];
         Rimhihi_d[i] = new double[dim];
         Rimlohi_d[i] = new double[dim];
         Rimhilo_d[i] = new double[dim];
         Rimlolo_d[i] = new double[dim];
      }
   }
   return 0;
}

void cmplx4_start_setup
 ( int dim, int deg,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d, int mode, int vrblvl )
{
   double *start0rehihi = new double[dim];
   double *start0relohi = new double[dim];
   double *start0rehilo = new double[dim];
   double *start0relolo = new double[dim];
   double *start0imhihi = new double[dim];
   double *start0imlohi = new double[dim];
   double *start0imhilo = new double[dim];
   double *start0imlolo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0rehihi[i] = testsolrehihi[i][0];
      start0relohi[i] = testsolrelohi[i][0];
      start0rehilo[i] = testsolrehilo[i][0];
      start0relolo[i] = testsolrelolo[i][0];
      start0imhihi[i] = testsolimhihi[i][0]; 
      start0imlohi[i] = testsolimlohi[i][0]; 
      start0imhilo[i] = testsolimhilo[i][0]; 
      start0imlolo[i] = testsolimlolo[i][0]; 
   }
   if((mode == 1) || (mode == 2))
      cmplx4_start_series_vector
         (dim,deg,start0rehihi,start0relohi,start0rehilo,start0relolo,
                  start0imhihi,start0imlohi,start0imhilo,start0imlolo,
          inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
          inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h);
   else
      cmplx4_start_series_vector
         (dim,deg,start0rehihi,start0relohi,start0rehilo,start0relolo,
                  start0imhihi,start0imlohi,start0imhilo,start0imlolo,
          inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
          inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d);

   if(mode == 2)
   {
      for(int i=0; i<dim; i++)
         for(int j=0; j<=deg; j++)
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
   }
   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;

      if((mode == 1) || (mode == 2))
      {
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
      else
      {
         for(int i=0; i<dim; i++)
            cout << i << " : "
                 << inputrehihi_d[i][0] << "  "
                 << inputrelohi_d[i][0] << endl << "  "
                 << inputrehilo_d[i][0] << "  "
                 << inputrelolo_d[i][0] << endl << "  "
                 << inputimhihi_d[i][0] << "  "
                 << inputimlohi_d[i][0] << endl << "  "
                 << inputimhilo_d[i][0] << "  "
                 << inputimlolo_d[i][0] << endl;
      }
   }
}

void cmplx4_column_setup
 ( int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **rowsA,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
   double **mbrhsrehihi, double **mbrhsrelohi,
   double **mbrhsrehilo, double **mbrhsrelolo,
   double **mbrhsimhihi, double **mbrhsimlohi,
   double **mbrhsimhilo, double **mbrhsimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolrehihi[i] = new double[degp1];
      testsolrelohi[i] = new double[degp1];
      testsolrehilo[i] = new double[degp1];
      testsolrelolo[i] = new double[degp1];
      testsolimhihi[i] = new double[degp1];
      testsolimlohi[i] = new double[degp1];
      testsolimhilo[i] = new double[degp1];
      testsolimlolo[i] = new double[degp1];
   }
   make_complex4_exponentials
      (dim,deg,testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
               testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo);

   // compute the right hand sides via evaluation

   for(int i=0; i<dim; i++)
   {
      mbrhsrehihi[i] = new double[degp1];
      mbrhsrelohi[i] = new double[degp1];
      mbrhsrehilo[i] = new double[degp1];
      mbrhsrelolo[i] = new double[degp1];
      mbrhsimhihi[i] = new double[degp1];
      mbrhsimlohi[i] = new double[degp1];
      mbrhsimhilo[i] = new double[degp1];
      mbrhsimlolo[i] = new double[degp1];

      mbrhsrehihi[i][0] = 1.0;     // initialize product to one
      mbrhsrelohi[i][0] = 0.0;
      mbrhsrehilo[i][0] = 0.0; 
      mbrhsrelolo[i][0] = 0.0;
      mbrhsimhihi[i][0] = 0.0;
      mbrhsimlohi[i][0] = 0.0;
      mbrhsimhilo[i][0] = 0.0;
      mbrhsimlolo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhsrehihi[i][k] = 0.0; mbrhsrelohi[i][k] = 0.0;
         mbrhsrehilo[i][k] = 0.0; mbrhsrelolo[i][k] = 0.0;
         mbrhsimhihi[i][k] = 0.0; mbrhsimlohi[i][k] = 0.0;
         mbrhsimhilo[i][k] = 0.0; mbrhsimlolo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_complex4_monomials
         (dim,deg,rowsA,
          testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
          testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
          mbrhsrehihi,mbrhsrelohi,mbrhsrehilo,mbrhsrelolo,
          mbrhsimhihi,mbrhsimlohi,mbrhsimhilo,mbrhsimlolo);
   else
      evaluate_complex4_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,
          testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
          testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
          mbrhsrehihi,mbrhsrelohi,mbrhsrehilo,mbrhsrelolo,
          mbrhsimhihi,mbrhsimlohi,mbrhsimhilo,mbrhsimlolo,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsrehihi[i][j] << "  " << mbrhsrelohi[i][j] << endl
                 << "  "
                 << mbrhsrehilo[i][j] << "  " << mbrhsrelolo[i][j] << endl
                 << "  "
                 << mbrhsimhihi[i][j] << "  " << mbrhsimlohi[i][j] << endl
                 << "  "
                 << mbrhsimhilo[i][j] << "  " << mbrhsimlolo[i][j] << endl;
   }
   cmplx4_start_setup(dim,deg,
      testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
      testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
      inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
      inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
      inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
      inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,mode,vrblvl);
}

void cmplx4_row_setup
 ( int dim, int deg, int *nbr, int **nvr, int ***idx,
   double **cstrehihi, double **cstrelohi, 
   double **cstrehilo, double **cstrelolo,
   double **cstimhihi, double **cstimlohi,
   double **cstimhilo, double **cstimlolo,
   double ***cffrehihi, double ***cffrelohi,
   double ***cffrehilo, double ***cffrelolo,
   double ***cffimhihi, double ***cffimlohi,
   double ***cffimhilo, double ***cffimlolo,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
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
   double ***outputimhilo_d, double ***outputimlolo_d, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   for(int i=0; i<dim; i++)
   {
      testsolrehihi[i] = new double[degp1];
      testsolrelohi[i] = new double[degp1];
      testsolrehilo[i] = new double[degp1];
      testsolrelolo[i] = new double[degp1];
      testsolimhihi[i] = new double[degp1];
      testsolimlohi[i] = new double[degp1];
      testsolimhilo[i] = new double[degp1];
      testsolimlolo[i] = new double[degp1];
   }
   make_complex4_exponentials(dim,deg,
      testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
      testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo);

   if(mode == 1)
   {
      if(vrblvl > 0)
         cout << "Evaluating test solution on the host ..." << endl;

      double timelapsed_h = 0.0;

      for(int i=0; i<dim; i++)
      {
         double lapsed;

         CPU_cmplx4_poly_evaldiff
           (dim,nbr[i],deg,nvr[i],idx[i],
            cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
            cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
            cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
            cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
            testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
            testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
            outputrehihi_h[i],outputrelohi_h[i],
            outputrehilo_h[i],outputrelolo_h[i],
            outputimhihi_h[i],outputimlohi_h[i],
            outputimhilo_h[i],outputimlolo_h[i],&lapsed,0); // no output

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
         for(int j=0; j<=deg; j++)
         {
            // cstre[i][j] = cstre[i][j] - outputre_h[i][dim][j];
            qdf_dec(&cstrehihi[i][j],&cstrelohi[i][j],
                    &cstrehilo[i][j],&cstrelolo[i][j],
                    outputrehihi_h[i][dim][j],outputrelohi_h[i][dim][j],
                    outputrehilo_h[i][dim][j],outputrelolo_h[i][dim][j]);
            // cstim[i][j] = cstim[i][j] - outputim_h[i][dim][j];
            qdf_dec(&cstimhihi[i][j],&cstimlohi[i][j],
                    &cstimhilo[i][j],&cstimlolo[i][j],
                    outputimhihi_h[i][dim][j],outputimlohi_h[i][dim][j],
                    outputimhilo_h[i][dim][j],outputimlolo_h[i][dim][j]);
         }
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double lapsed;

            CPU_cmplx4_poly_evaldiff
              (dim,nbr[i],deg,nvr[i],idx[i],
               cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
               cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
               cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
               cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
               testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
               testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
               outputrehihi_h[i],outputrelohi_h[i],
               outputrehilo_h[i],outputrelolo_h[i],
               outputimhihi_h[i],outputimlohi_h[i],
               outputimhilo_h[i],outputimlolo_h[i],&lapsed,0); // no output

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
               errsum = errsum + outputrehihi_h[i][dim][j]
                               + outputrelohi_h[i][dim][j]
                               + outputrehilo_h[i][dim][j]
                               + outputrelolo_h[i][dim][j]
                               + outputimhihi_h[i][dim][j]
                               + outputimlohi_h[i][dim][j]
                               + outputimhilo_h[i][dim][j]
                               + outputimlolo_h[i][dim][j]; 

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

         ComplexConvolutionJobs cnvjobs(dim);
         ComplexIncrementJobs incjobs(cnvjobs,false);
         ComplexAdditionJobs addjobs(dim,nbr[i]);

         make_all_complex_jobs
            (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrb);

         GPU_cmplx4vectorized_poly_evaldiff
            (degp1,dim,nbr[i],deg,nvr[i],idx[i],
             cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
             cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
             cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
             cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
             testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
             testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
             outputrehihi_d[i],outputrelohi_d[i],
             outputrehilo_d[i],outputrelolo_d[i],
             outputimhihi_d[i],outputimlohi_d[i],
             outputimhilo_d[i],outputimlolo_d[i],cnvjobs,incjobs,addjobs,
             &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0); // vrblvl);

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
         for(int j=0; j<=deg; j++)
         {
            // cstre[i][j] = cstre[i][j] - outputre_d[i][dim][j];
            qdf_dec(&cstrehihi[i][j],&cstrelohi[i][j],
                    &cstrehilo[i][j],&cstrelolo[i][j],
                    outputrehihi_d[i][dim][j],outputrelohi_d[i][dim][j],
                    outputrehilo_d[i][dim][j],outputrelolo_d[i][dim][j]);
            // cstim[i][j] = cstim[i][j] - outputim_d[i][dim][j];
            qdf_dec(&cstimhihi[i][j],&cstimlohi[i][j],
                    &cstimhilo[i][j],&cstimlolo[i][j],
                    outputimhihi_d[i][dim][j],outputimlohi_d[i][dim][j],
                    outputimhilo_d[i][dim][j],outputimlolo_d[i][dim][j]);
         }
      }
      if(vrblvl > 1) // evaluate again to test
      {
         cout << "Evaluating again to compute the residual ..." << endl;

         for(int i=0; i<dim; i++)
         {
            double cnvlapms,addlapms,timelapms_d,walltimes_d;
   
            ComplexConvolutionJobs cnvjobs(dim);
            ComplexIncrementJobs incjobs(cnvjobs,false);
            ComplexAdditionJobs addjobs(dim,nbr[i]);

            make_all_complex_jobs
               (dim,nbr[i],nvr[i],idx[i],&cnvjobs,&incjobs,&addjobs,vrb);

            GPU_cmplx4vectorized_poly_evaldiff
               (degp1,dim,nbr[i],deg,nvr[i],idx[i],
                cstrehihi[i],cstrelohi[i],cstrehilo[i],cstrelolo[i],
                cstimhihi[i],cstimlohi[i],cstimhilo[i],cstimlolo[i],
                cffrehihi[i],cffrelohi[i],cffrehilo[i],cffrelolo[i],
                cffimhihi[i],cffimlohi[i],cffimhilo[i],cffimlolo[i],
                testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
                testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
                outputrehihi_d[i],outputrelohi_d[i],
                outputrehilo_d[i],outputrelolo_d[i],
                outputimhihi_d[i],outputimlohi_d[i],
                outputimhilo_d[i],outputimlolo_d[i],cnvjobs,incjobs,addjobs,
                &cnvlapms,&addlapms,&timelapms_d,&walltimes_d,0); // vrblvl);

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
               errsum = errsum + outputrehihi_d[i][dim][j]
                               + outputrelohi_d[i][dim][j]
                               + outputrehilo_d[i][dim][j]
                               + outputrelolo_d[i][dim][j]
                               + outputimhihi_d[i][dim][j]
                               + outputimlohi_d[i][dim][j]
                               + outputimhilo_d[i][dim][j]
                               + outputimlolo_d[i][dim][j];

         cout << scientific << setprecision(2)
              << "Residual of test solution : " << errsum << endl;
      }
   }
   cmplx4_start_setup(dim,deg,
      testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
      testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
      inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
      inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
      inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
      inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,mode,vrblvl);
}

int cmplx4_error_testsol
 ( int dim, int deg, int mode,
   double **testsolrehihi, double **testsolrelohi,
   double **testsolrehilo, double **testsolrelolo,
   double **testsolimhihi, double **testsolimlohi,
   double **testsolimhilo, double **testsolimlolo,
   double **inputrehihi_h, double **inputrelohi_h,
   double **inputrehilo_h, double **inputrelolo_h,
   double **inputimhihi_h, double **inputimlohi_h,
   double **inputimhilo_h, double **inputimlolo_h,
   double **inputrehihi_d, double **inputrelohi_d,
   double **inputrehilo_d, double **inputrelolo_d,
   double **inputimhihi_d, double **inputimlohi_d,
   double **inputimhilo_d, double **inputimlolo_d )
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
                        << testsolrehihi[i][j] << "  "
                        << testsolrelohi[i][j] << endl << "  "
                        << testsolrehilo[i][j] << "  "
                        << testsolrelolo[i][j] << endl << "  "
                        << testsolimhihi[i][j] << "  "
                        << testsolimlohi[i][j] << endl << "  "
                        << testsolimhilo[i][j] << "  "
                        << testsolimlolo[i][j] << endl;

         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehihi_d[i][j] << "  "
                           << inputrelohi_d[i][j] << endl << "  "
                           << inputrehilo_d[i][j] << "  "
                           << inputrelolo_d[i][j] << endl << "  "
                           << inputimhihi_d[i][j] << "  "
                           << inputimlohi_d[i][j] << endl << "  "
                           << inputimhilo_d[i][j] << "  "
                           << inputimlolo_d[i][j] << endl;

            errsum += abs(testsolrehihi[i][j] - inputrehihi_d[i][j])
                    + abs(testsolrelohi[i][j] - inputrelohi_d[i][j])
                    + abs(testsolrehilo[i][j] - inputrehilo_d[i][j])
                    + abs(testsolrelolo[i][j] - inputrelolo_d[i][j])
                    + abs(testsolimhihi[i][j] - inputimhihi_d[i][j])
                    + abs(testsolimlohi[i][j] - inputimlohi_d[i][j])
                    + abs(testsolimhilo[i][j] - inputimhilo_d[i][j])
                    + abs(testsolimlolo[i][j] - inputimlolo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehihi_h[i][j] << "  "
                           << inputrelohi_h[i][j] << endl << "  "
                           << inputrehilo_h[i][j] << "  "
                           << inputrelolo_h[i][j] << endl << "  "
                           << inputimhihi_h[i][j] << "  "
                           << inputimlohi_h[i][j] << endl << "  "
                           << inputimhilo_h[i][j] << "  "
                           << inputimlolo_h[i][j] << endl;

            errsum += abs(testsolrehihi[i][j] - inputrehihi_h[i][j])
                    + abs(testsolrelohi[i][j] - inputrelohi_h[i][j])
                    + abs(testsolrehilo[i][j] - inputrehilo_h[i][j])
                    + abs(testsolrelolo[i][j] - inputrelolo_h[i][j])
                    + abs(testsolimhihi[i][j] - inputimhihi_h[i][j])
                    + abs(testsolimlohi[i][j] - inputimlohi_h[i][j])
                    + abs(testsolimhilo[i][j] - inputimhilo_h[i][j])
                    + abs(testsolimlolo[i][j] - inputimlolo_h[i][j]);
         }
      }
   }
   cout << "error : " << errsum << endl;

   return (errsum > 1.0e-50);
}

int test_cmplx4_column_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;

   double **accrehihi = new double*[dim+1]; // accumulated series
   double **accrelohi = new double*[dim+1]; // in one column
   double **accrehilo = new double*[dim+1];
   double **accrelolo = new double*[dim+1];
   double **accimhihi = new double*[dim+1];
   double **accimlohi = new double*[dim+1];
   double **accimhilo = new double*[dim+1];
   double **accimlolo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehihi[i] = new double[degp1];
      accrelohi[i] = new double[degp1];
      accrehilo[i] = new double[degp1];
      accrelolo[i] = new double[degp1];
      accimhihi[i] = new double[degp1];
      accimlohi[i] = new double[degp1];
      accimhilo[i] = new double[degp1];
      accimlolo[i] = new double[degp1];
   }
   double ***cffrehihi = new double**[nbrcol]; // coefficients of monomials
   double ***cffrelohi = new double**[nbrcol];
   double ***cffrehilo = new double**[nbrcol]; 
   double ***cffrelolo = new double**[nbrcol];
   double ***cffimhihi = new double**[nbrcol]; 
   double ***cffimlohi = new double**[nbrcol]; 
   double ***cffimhilo = new double**[nbrcol]; 
   double ***cffimlolo = new double**[nbrcol]; 

   for(int i=0; i<nbrcol; i++)
   {
      cffrehihi[i] = new double*[dim];
      cffrelohi[i] = new double*[dim];
      cffrehilo[i] = new double*[dim];
      cffrelolo[i] = new double*[dim];
      cffimhihi[i] = new double*[dim];
      cffimlohi[i] = new double*[dim];
      cffimhilo[i] = new double*[dim];
      cffimlolo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffrehihi[i][j] = new double[degp1];
         cffrelohi[i][j] = new double[degp1];
         cffrehilo[i][j] = new double[degp1];
         cffrelolo[i][j] = new double[degp1];
         cffimhihi[i][j] = new double[degp1];
         cffimlohi[i][j] = new double[degp1];
         cffimhilo[i][j] = new double[degp1];
         cffimlolo[i][j] = new double[degp1];
      }
   }
   if(nbrcol != 1) // generate coefficients for the columns
      make_complex4_coefficients
         (nbrcol,dim,cffrehihi,cffrelohi,cffrehilo,cffrelolo,
                     cffimhihi,cffimlohi,cffimhilo,cffimlolo);

   double **inputrehihi_h;
   double **inputrelohi_h;
   double **inputrehilo_h;
   double **inputrelolo_h;
   double **inputimhihi_h;
   double **inputimlohi_h;
   double **inputimhilo_h;
   double **inputimlolo_h;
   double **inputrehihi_d;
   double **inputrelohi_d;
   double **inputrehilo_d;
   double **inputrelolo_d;
   double **inputimhihi_d;
   double **inputimlohi_d;
   double **inputimhilo_d;
   double **inputimlolo_d;
   double ***outputrehihi_h;
   double ***outputrelohi_h;
   double ***outputrehilo_h;
   double ***outputrelolo_h;
   double ***outputimhihi_h;
   double ***outputimlohi_h;
   double ***outputimhilo_h;
   double ***outputimlolo_h;
   double ***outputrehihi_d;
   double ***outputrelohi_d;
   double ***outputrehilo_d;
   double ***outputrelolo_d;
   double ***outputimhihi_d;
   double ***outputimlohi_d;
   double ***outputimhilo_d;
   double ***outputimlolo_d;
   double **funvalrehihi_h;
   double **funvalrelohi_h;
   double **funvalrehilo_h;
   double **funvalrelolo_h;
   double **funvalimhihi_h;
   double **funvalimlohi_h;
   double **funvalimhilo_h;
   double **funvalimlolo_h;
   double **funvalrehihi_d;
   double **funvalrelohi_d;
   double **funvalrehilo_d;
   double **funvalrelolo_d;
   double **funvalimhihi_d;
   double **funvalimlohi_d;
   double **funvalimhilo_d;
   double **funvalimlolo_d;
   double ***jacvalrehihi_h;
   double ***jacvalrelohi_h;
   double ***jacvalrehilo_h;
   double ***jacvalrelolo_h;
   double ***jacvalimhihi_h;
   double ***jacvalimlohi_h;
   double ***jacvalimhilo_h;
   double ***jacvalimlolo_h;
   double ***jacvalrehihi_d;
   double ***jacvalrelohi_d;
   double ***jacvalrehilo_d;
   double ***jacvalrelolo_d;
   double ***jacvalimhihi_d;
   double ***jacvalimlohi_d;
   double ***jacvalimhilo_d;
   double ***jacvalimlolo_d;

   if((mode == 1) || (mode == 2))
   {
      inputrehihi_h = new double*[dim];
      inputrelohi_h = new double*[dim];
      inputrehilo_h = new double*[dim];
      inputrelolo_h = new double*[dim];
      inputimhihi_h = new double*[dim];
      inputimlohi_h = new double*[dim];
      inputimhilo_h = new double*[dim];
      inputimlolo_h = new double*[dim];
      outputrehihi_h = new double**[dim];
      outputrelohi_h = new double**[dim];
      outputrehilo_h = new double**[dim];
      outputrelolo_h = new double**[dim];
      outputimhihi_h = new double**[dim];
      outputimlohi_h = new double**[dim];
      outputimhilo_h = new double**[dim];
      outputimlolo_h = new double**[dim];
      funvalrehihi_h = new double*[dim];
      funvalrelohi_h = new double*[dim];
      funvalrehilo_h = new double*[dim];
      funvalrelolo_h = new double*[dim];
      funvalimhihi_h = new double*[dim];
      funvalimlohi_h = new double*[dim];
      funvalimhilo_h = new double*[dim];
      funvalimlolo_h = new double*[dim];
      jacvalrehihi_h = new double**[degp1];
      jacvalrelohi_h = new double**[degp1];
      jacvalrehilo_h = new double**[degp1];
      jacvalrelolo_h = new double**[degp1];
      jacvalimhihi_h = new double**[degp1];
      jacvalimlohi_h = new double**[degp1];
      jacvalimhilo_h = new double**[degp1];
      jacvalimlolo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputrehihi_d = new double*[dim];
      inputrelohi_d = new double*[dim];
      inputrehilo_d = new double*[dim];
      inputrelolo_d = new double*[dim];
      inputimhihi_d = new double*[dim];
      inputimlohi_d = new double*[dim];
      inputimhilo_d = new double*[dim];
      inputimlolo_d = new double*[dim];
      outputrehihi_d = new double**[dim];
      outputrelohi_d = new double**[dim];
      outputrehilo_d = new double**[dim];
      outputrelolo_d = new double**[dim];
      outputimhihi_d = new double**[dim];
      outputimlohi_d = new double**[dim];
      outputimhilo_d = new double**[dim];
      outputimlolo_d = new double**[dim];
      funvalrehihi_d = new double*[dim];
      funvalrelohi_d = new double*[dim];
      funvalrehilo_d = new double*[dim];
      funvalrelolo_d = new double*[dim];
      funvalimhihi_d = new double*[dim];
      funvalimlohi_d = new double*[dim];
      funvalimhilo_d = new double*[dim];
      funvalimlolo_d = new double*[dim];
      jacvalrehihi_d = new double**[degp1];
      jacvalrelohi_d = new double**[degp1];
      jacvalrehilo_d = new double**[degp1];
      jacvalrelolo_d = new double**[degp1];
      jacvalimhihi_d = new double**[degp1];
      jacvalimlohi_d = new double**[degp1];
      jacvalimhilo_d = new double**[degp1];
      jacvalimlolo_d = new double**[degp1];
   }
   cmplx4_allocate_inoutfunjac(dim,deg,mode,
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
       jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhsrehihi_h;
   double **rhsrelohi_h;
   double **rhsrehilo_h;
   double **rhsrelolo_h;
   double **rhsimhihi_h;
   double **rhsimlohi_h;
   double **rhsimhilo_h;
   double **rhsimlolo_h;
   double **rhsrehihi_d;
   double **rhsrelohi_d;
   double **rhsrehilo_d;
   double **rhsrelolo_d;
   double **rhsimhihi_d;
   double **rhsimlohi_d;
   double **rhsimhilo_d;
   double **rhsimlolo_d;
   double **urhsrehihi_h;
   double **urhsrelohi_h;
   double **urhsrehilo_h;
   double **urhsrelolo_h;
   double **urhsimhihi_h;
   double **urhsimlohi_h;
   double **urhsimhilo_h;
   double **urhsimlolo_h;
   double **urhsrehihi_d;
   double **urhsrelohi_d;
   double **urhsrehilo_d;
   double **urhsrelolo_d;
   double **urhsimhihi_d;
   double **urhsimlohi_d;
   double **urhsimhilo_d;
   double **urhsimlolo_d;
   double **solrehihi_h;
   double **solrelohi_h;
   double **solrehilo_h;
   double **solrelolo_h;
   double **solimhihi_h;
   double **solimlohi_h;
   double **solimhilo_h;
   double **solimlolo_h;
   double **solrehihi_d;
   double **solrelohi_d;
   double **solrehilo_d;
   double **solrelolo_d;
   double **solimhihi_d;
   double **solimlohi_d;
   double **solimhilo_d;
   double **solimlolo_d;
   double **Qrehihi_h;
   double **Qrelohi_h;
   double **Qrehilo_h;
   double **Qrelolo_h;
   double **Qimhihi_h;
   double **Qimlohi_h;
   double **Qimhilo_h;
   double **Qimlolo_h;
   double **Qrehihi_d;
   double **Qrelohi_d;
   double **Qrehilo_d;
   double **Qrelolo_d;
   double **Qimhihi_d;
   double **Qimlohi_d;
   double **Qimhilo_d;
   double **Qimlolo_d;
   double **Rrehihi_h;
   double **Rrelohi_h;
   double **Rrehilo_h;
   double **Rrelolo_h;
   double **Rimhihi_h;
   double **Rimlohi_h;
   double **Rimhilo_h;
   double **Rimlolo_h;
   double **Rrehihi_d;
   double **Rrelohi_d;
   double **Rrehilo_d;
   double **Rrelolo_d;
   double **Rimhihi_d;
   double **Rimlohi_d;
   double **Rimhilo_d;
   double **Rimlolo_d;

   if((mode == 1) || (mode == 2))
   {
      rhsrehihi_h = new double*[degp1];
      rhsrelohi_h = new double*[degp1];
      rhsrehilo_h = new double*[degp1];
      rhsrelolo_h = new double*[degp1];
      rhsimhihi_h = new double*[degp1];
      rhsimlohi_h = new double*[degp1];
      rhsimhilo_h = new double*[degp1];
      rhsimlolo_h = new double*[degp1];
      urhsrehihi_h = new double*[degp1];
      urhsrelohi_h = new double*[degp1];
      urhsrehilo_h = new double*[degp1];
      urhsrelolo_h = new double*[degp1];
      urhsimhihi_h = new double*[degp1];
      urhsimlohi_h = new double*[degp1];
      urhsimhilo_h = new double*[degp1];
      urhsimlolo_h = new double*[degp1];
      solrehihi_h = new double*[degp1];
      solrelohi_h = new double*[degp1];
      solrehilo_h = new double*[degp1];
      solrelolo_h = new double*[degp1];
      solimhihi_h = new double*[degp1];
      solimlohi_h = new double*[degp1];
      solimhilo_h = new double*[degp1];
      solimlolo_h = new double*[degp1];
      Qrehihi_h = new double*[dim];
      Qrelohi_h = new double*[dim];
      Qrehilo_h = new double*[dim];
      Qrelolo_h = new double*[dim];
      Qimhihi_h = new double*[dim];
      Qimlohi_h = new double*[dim];
      Qimhilo_h = new double*[dim];
      Qimlolo_h = new double*[dim];
      Rrehihi_h = new double*[dim];
      Rrelohi_h = new double*[dim];
      Rrehilo_h = new double*[dim];
      Rrelolo_h = new double*[dim];
      Rimhihi_h = new double*[dim];
      Rimlohi_h = new double*[dim];
      Rimhilo_h = new double*[dim];
      Rimlolo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhsrehihi_d = new double*[degp1];
      rhsrelohi_d = new double*[degp1];
      rhsrehilo_d = new double*[degp1];
      rhsrelolo_d = new double*[degp1];
      rhsimhihi_d = new double*[degp1];
      rhsimlohi_d = new double*[degp1];
      rhsimhilo_d = new double*[degp1];
      rhsimlolo_d = new double*[degp1];
      urhsrehihi_d = new double*[degp1];
      urhsrelohi_d = new double*[degp1];
      urhsrehilo_d = new double*[degp1];
      urhsrelolo_d = new double*[degp1];
      urhsimhihi_d = new double*[degp1];
      urhsimlohi_d = new double*[degp1];
      urhsimhilo_d = new double*[degp1];
      urhsimlolo_d = new double*[degp1];
      solrehihi_d = new double*[degp1];
      solrelohi_d = new double*[degp1];
      solrehilo_d = new double*[degp1];
      solrelolo_d = new double*[degp1];
      solimhihi_d = new double*[degp1];
      solimlohi_d = new double*[degp1];
      solimhilo_d = new double*[degp1];
      solimlolo_d = new double*[degp1];
      Qrehihi_d = new double*[dim];
      Qrelohi_d = new double*[dim];
      Qrehilo_d = new double*[dim];
      Qrelolo_d = new double*[dim];
      Qimhihi_d = new double*[dim];
      Qimlohi_d = new double*[dim];
      Qimhilo_d = new double*[dim];
      Qimlolo_d = new double*[dim];
      Rrehihi_d = new double*[dim];
      Rrelohi_d = new double*[dim];
      Rrehilo_d = new double*[dim];
      Rrelolo_d = new double*[dim];
      Rimhihi_d = new double*[dim];
      Rimlohi_d = new double*[dim];
      Rimhilo_d = new double*[dim];
      Rimlolo_d = new double*[dim];
   }
   cmplx4_allocate_rhsqrsol(dim,deg,mode,
       rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
       rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
       rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
       rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
       urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
       urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
       urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
       urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
       Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
       Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
       Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
       Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
       Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
       Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
       Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
       Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
       solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
       solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
       solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
       solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d);

   double *workvecrehihi = new double[dim];
   double *workvecrelohi = new double[dim];
   double *workvecrehilo = new double[dim];
   double *workvecrelolo = new double[dim];
   double *workvecimhihi = new double[dim];
   double *workvecimlohi = new double[dim];
   double *workvecimhilo = new double[dim];
   double *workvecimlolo = new double[dim];

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
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolrehihi = new double*[dim];
   double **testsolrelohi = new double*[dim];
   double **testsolrehilo = new double*[dim];
   double **testsolrelolo = new double*[dim];
   double **testsolimhihi = new double*[dim];
   double **testsolimlohi = new double*[dim];
   double **testsolimhilo = new double*[dim];
   double **testsolimlolo = new double*[dim];
   double **mbrhsrehihi = new double*[dim];
   double **mbrhsrelohi = new double*[dim];
   double **mbrhsrehilo = new double*[dim];
   double **mbrhsrelolo = new double*[dim];
   double **mbrhsimhihi = new double*[dim];
   double **mbrhsimlohi = new double*[dim];
   double **mbrhsimhilo = new double*[dim];
   double **mbrhsimlolo = new double*[dim];

   cmplx4_column_setup
      (dim,deg,nbrcol,nvr,idx,rowsA,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
       testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
       mbrhsrehihi,mbrhsrelohi,mbrhsrehilo,mbrhsrelolo,
       mbrhsimhihi,mbrhsimlohi,mbrhsimhilo,mbrhsimlolo,
       inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
       inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
       inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
       inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,mode,vrblvl);

   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = false;
   bool zeroQ_d = false;
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

      cmplx4_column_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,
          &tailidx_h,&tailidx_d,nvr,idx,exp,nbrfac,expfac,
          mbrhsrehihi,mbrhsrelohi,mbrhsrehilo,mbrhsrelolo,
          mbrhsimhihi,mbrhsimlohi,mbrhsimhilo,mbrhsimlolo,dpr,
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
          workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
          workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,
          resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
          resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
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

   cmplx4_error_testsol(dim,deg,mode,
      testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
      testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
      inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
      inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
      inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
      inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}

int test_cmplx4_row_newton
 ( int szt, int nbt, int dim, int deg, int *nbr, int **nvr, int ***idx,
   double dpr, int nbsteps, int mode, int vrblvl )
{
   const int degp1 = deg+1;

   double **cstrehihi = new double*[dim];
   double **cstrelohi = new double*[dim];
   double **cstrehilo = new double*[dim]; 
   double **cstrelolo = new double*[dim];
   double **cstimhihi = new double*[dim]; 
   double **cstimlohi = new double*[dim]; 
   double **cstimhilo = new double*[dim]; 
   double **cstimlolo = new double*[dim]; 
   double ***cffrehihi = new double**[dim];
   double ***cffrelohi = new double**[dim];
   double ***cffrehilo = new double**[dim]; 
   double ***cffrelolo = new double**[dim];
   double ***cffimhihi = new double**[dim]; 
   double ***cffimlohi = new double**[dim]; 
   double ***cffimhilo = new double**[dim]; 
   double ***cffimlolo = new double**[dim]; 

   cmplx4_make_coefficients(dim,deg,nbr,nvr,idx,
      cstrehihi,cstrelohi,cstrehilo,cstrelolo,
      cstimhihi,cstimlohi,cstimhilo,cstimlolo,
      cffrehihi,cffrelohi,cffrehilo,cffrelolo,
      cffimhihi,cffimlohi,cffimhilo,cffimlolo,vrblvl);

   double **inputrehihi_h;
   double **inputrelohi_h;
   double **inputrehilo_h;
   double **inputrelolo_h;
   double **inputimhihi_h;
   double **inputimlohi_h;
   double **inputimhilo_h;
   double **inputimlolo_h;
   double **inputrehihi_d;
   double **inputrelohi_d;
   double **inputrehilo_d;
   double **inputrelolo_d;
   double **inputimhihi_d;
   double **inputimlohi_d;
   double **inputimhilo_d;
   double **inputimlolo_d;
   double ***outputrehihi_h;
   double ***outputrelohi_h;
   double ***outputrehilo_h;
   double ***outputrelolo_h;
   double ***outputimhihi_h;
   double ***outputimlohi_h;
   double ***outputimhilo_h;
   double ***outputimlolo_h;
   double ***outputrehihi_d;
   double ***outputrelohi_d;
   double ***outputrehilo_d;
   double ***outputrelolo_d;
   double ***outputimhihi_d;
   double ***outputimlohi_d;
   double ***outputimhilo_d;
   double ***outputimlolo_d;
   double **funvalrehihi_h;
   double **funvalrelohi_h;
   double **funvalrehilo_h;
   double **funvalrelolo_h;
   double **funvalimhihi_h;
   double **funvalimlohi_h;
   double **funvalimhilo_h;
   double **funvalimlolo_h;
   double **funvalrehihi_d;
   double **funvalrelohi_d;
   double **funvalrehilo_d;
   double **funvalrelolo_d;
   double **funvalimhihi_d;
   double **funvalimlohi_d;
   double **funvalimhilo_d;
   double **funvalimlolo_d;
   double ***jacvalrehihi_h;
   double ***jacvalrelohi_h;
   double ***jacvalrehilo_h;
   double ***jacvalrelolo_h;
   double ***jacvalimhihi_h;
   double ***jacvalimlohi_h;
   double ***jacvalimhilo_h;
   double ***jacvalimlolo_h;
   double ***jacvalrehihi_d;
   double ***jacvalrelohi_d;
   double ***jacvalrehilo_d;
   double ***jacvalrelolo_d;
   double ***jacvalimhihi_d;
   double ***jacvalimlohi_d;
   double ***jacvalimhilo_d;
   double ***jacvalimlolo_d;

   if((mode == 1) || (mode == 2))
   {
      inputrehihi_h = new double*[dim];
      inputrelohi_h = new double*[dim];
      inputrehilo_h = new double*[dim];
      inputrelolo_h = new double*[dim];
      inputimhihi_h = new double*[dim];
      inputimlohi_h = new double*[dim];
      inputimhilo_h = new double*[dim];
      inputimlolo_h = new double*[dim];
      outputrehihi_h = new double**[dim];
      outputrelohi_h = new double**[dim];
      outputrehilo_h = new double**[dim];
      outputrelolo_h = new double**[dim];
      outputimhihi_h = new double**[dim];
      outputimlohi_h = new double**[dim];
      outputimhilo_h = new double**[dim];
      outputimlolo_h = new double**[dim];
      funvalrehihi_h = new double*[dim];
      funvalrelohi_h = new double*[dim];
      funvalrehilo_h = new double*[dim];
      funvalrelolo_h = new double*[dim];
      funvalimhihi_h = new double*[dim];
      funvalimlohi_h = new double*[dim];
      funvalimhilo_h = new double*[dim];
      funvalimlolo_h = new double*[dim];
      jacvalrehihi_h = new double**[degp1];
      jacvalrelohi_h = new double**[degp1];
      jacvalrehilo_h = new double**[degp1];
      jacvalrelolo_h = new double**[degp1];
      jacvalimhihi_h = new double**[degp1];
      jacvalimlohi_h = new double**[degp1];
      jacvalimhilo_h = new double**[degp1];
      jacvalimlolo_h = new double**[degp1];
   }
   if((mode == 0) || (mode == 2))
   {
      inputrehihi_d = new double*[dim];
      inputrelohi_d = new double*[dim];
      inputrehilo_d = new double*[dim];
      inputrelolo_d = new double*[dim];
      inputimhihi_d = new double*[dim];
      inputimlohi_d = new double*[dim];
      inputimhilo_d = new double*[dim];
      inputimlolo_d = new double*[dim];
      outputrehihi_d = new double**[dim];
      outputrelohi_d = new double**[dim];
      outputrehilo_d = new double**[dim];
      outputrelolo_d = new double**[dim];
      outputimhihi_d = new double**[dim];
      outputimlohi_d = new double**[dim];
      outputimhilo_d = new double**[dim];
      outputimlolo_d = new double**[dim];
      funvalrehihi_d = new double*[dim];
      funvalrelohi_d = new double*[dim];
      funvalrehilo_d = new double*[dim];
      funvalrelolo_d = new double*[dim];
      funvalimhihi_d = new double*[dim];
      funvalimlohi_d = new double*[dim];
      funvalimhilo_d = new double*[dim];
      funvalimlolo_d = new double*[dim];
      jacvalrehihi_d = new double**[degp1];
      jacvalrelohi_d = new double**[degp1];
      jacvalrehilo_d = new double**[degp1];
      jacvalrelolo_d = new double**[degp1];
      jacvalimhihi_d = new double**[degp1];
      jacvalimlohi_d = new double**[degp1];
      jacvalimhilo_d = new double**[degp1];
      jacvalimlolo_d = new double**[degp1];
   }
   cmplx4_allocate_inoutfunjac(dim,deg,mode,
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
       jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d);
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **rhsrehihi_h;
   double **rhsrelohi_h;
   double **rhsrehilo_h;
   double **rhsrelolo_h;
   double **rhsimhihi_h;
   double **rhsimlohi_h;
   double **rhsimhilo_h;
   double **rhsimlolo_h;
   double **rhsrehihi_d;
   double **rhsrelohi_d;
   double **rhsrehilo_d;
   double **rhsrelolo_d;
   double **rhsimhihi_d;
   double **rhsimlohi_d;
   double **rhsimhilo_d;
   double **rhsimlolo_d;
   double **urhsrehihi_h;
   double **urhsrelohi_h;
   double **urhsrehilo_h;
   double **urhsrelolo_h;
   double **urhsimhihi_h;
   double **urhsimlohi_h;
   double **urhsimhilo_h;
   double **urhsimlolo_h;
   double **urhsrehihi_d;
   double **urhsrelohi_d;
   double **urhsrehilo_d;
   double **urhsrelolo_d;
   double **urhsimhihi_d;
   double **urhsimlohi_d;
   double **urhsimhilo_d;
   double **urhsimlolo_d;
   double **solrehihi_h;
   double **solrelohi_h;
   double **solrehilo_h;
   double **solrelolo_h;
   double **solimhihi_h;
   double **solimlohi_h;
   double **solimhilo_h;
   double **solimlolo_h;
   double **solrehihi_d;
   double **solrelohi_d;
   double **solrehilo_d;
   double **solrelolo_d;
   double **solimhihi_d;
   double **solimlohi_d;
   double **solimhilo_d;
   double **solimlolo_d;
   double **Qrehihi_h;
   double **Qrelohi_h;
   double **Qrehilo_h;
   double **Qrelolo_h;
   double **Qimhihi_h;
   double **Qimlohi_h;
   double **Qimhilo_h;
   double **Qimlolo_h;
   double **Qrehihi_d;
   double **Qrelohi_d;
   double **Qrehilo_d;
   double **Qrelolo_d;
   double **Qimhihi_d;
   double **Qimlohi_d;
   double **Qimhilo_d;
   double **Qimlolo_d;
   double **Rrehihi_h;
   double **Rrelohi_h;
   double **Rrehilo_h;
   double **Rrelolo_h;
   double **Rimhihi_h;
   double **Rimlohi_h;
   double **Rimhilo_h;
   double **Rimlolo_h;
   double **Rrehihi_d;
   double **Rrelohi_d;
   double **Rrehilo_d;
   double **Rrelolo_d;
   double **Rimhihi_d;
   double **Rimlohi_d;
   double **Rimhilo_d;
   double **Rimlolo_d;

   if((mode == 1) || (mode == 2))
   {
      rhsrehihi_h = new double*[degp1];
      rhsrelohi_h = new double*[degp1];
      rhsrehilo_h = new double*[degp1];
      rhsrelolo_h = new double*[degp1];
      rhsimhihi_h = new double*[degp1];
      rhsimlohi_h = new double*[degp1];
      rhsimhilo_h = new double*[degp1];
      rhsimlolo_h = new double*[degp1];
      urhsrehihi_h = new double*[degp1];
      urhsrelohi_h = new double*[degp1];
      urhsrehilo_h = new double*[degp1];
      urhsrelolo_h = new double*[degp1];
      urhsimhihi_h = new double*[degp1];
      urhsimlohi_h = new double*[degp1];
      urhsimhilo_h = new double*[degp1];
      urhsimlolo_h = new double*[degp1];
      solrehihi_h = new double*[degp1];
      solrelohi_h = new double*[degp1];
      solrehilo_h = new double*[degp1];
      solrelolo_h = new double*[degp1];
      solimhihi_h = new double*[degp1];
      solimlohi_h = new double*[degp1];
      solimhilo_h = new double*[degp1];
      solimlolo_h = new double*[degp1];
      Qrehihi_h = new double*[dim];
      Qrelohi_h = new double*[dim];
      Qrehilo_h = new double*[dim];
      Qrelolo_h = new double*[dim];
      Qimhihi_h = new double*[dim];
      Qimlohi_h = new double*[dim];
      Qimhilo_h = new double*[dim];
      Qimlolo_h = new double*[dim];
      Rrehihi_h = new double*[dim];
      Rrelohi_h = new double*[dim];
      Rrehilo_h = new double*[dim];
      Rrelolo_h = new double*[dim];
      Rimhihi_h = new double*[dim];
      Rimlohi_h = new double*[dim];
      Rimhilo_h = new double*[dim];
      Rimlolo_h = new double*[dim];
   }
   if((mode == 0) || (mode == 2))
   {
      rhsrehihi_d = new double*[degp1];
      rhsrelohi_d = new double*[degp1];
      rhsrehilo_d = new double*[degp1];
      rhsrelolo_d = new double*[degp1];
      rhsimhihi_d = new double*[degp1];
      rhsimlohi_d = new double*[degp1];
      rhsimhilo_d = new double*[degp1];
      rhsimlolo_d = new double*[degp1];
      urhsrehihi_d = new double*[degp1];
      urhsrelohi_d = new double*[degp1];
      urhsrehilo_d = new double*[degp1];
      urhsrelolo_d = new double*[degp1];
      urhsimhihi_d = new double*[degp1];
      urhsimlohi_d = new double*[degp1];
      urhsimhilo_d = new double*[degp1];
      urhsimlolo_d = new double*[degp1];
      solrehihi_d = new double*[degp1];
      solrelohi_d = new double*[degp1];
      solrehilo_d = new double*[degp1];
      solrelolo_d = new double*[degp1];
      solimhihi_d = new double*[degp1];
      solimlohi_d = new double*[degp1];
      solimhilo_d = new double*[degp1];
      solimlolo_d = new double*[degp1];
      Qrehihi_d = new double*[dim];
      Qrelohi_d = new double*[dim];
      Qrehilo_d = new double*[dim];
      Qrelolo_d = new double*[dim];
      Qimhihi_d = new double*[dim];
      Qimlohi_d = new double*[dim];
      Qimhilo_d = new double*[dim];
      Qimlolo_d = new double*[dim];
      Rrehihi_d = new double*[dim];
      Rrelohi_d = new double*[dim];
      Rrehilo_d = new double*[dim];
      Rrelolo_d = new double*[dim];
      Rimhihi_d = new double*[dim];
      Rimlohi_d = new double*[dim];
      Rimhilo_d = new double*[dim];
      Rimlolo_d = new double*[dim];
   }
   cmplx4_allocate_rhsqrsol(dim,deg,mode,
       rhsrehihi_h,rhsrelohi_h,rhsrehilo_h,rhsrelolo_h,
       rhsimhihi_h,rhsimlohi_h,rhsimhilo_h,rhsimlolo_h,
       rhsrehihi_d,rhsrelohi_d,rhsrehilo_d,rhsrelolo_d,
       rhsimhihi_d,rhsimlohi_d,rhsimhilo_d,rhsimlolo_d,
       urhsrehihi_h,urhsrelohi_h,urhsrehilo_h,urhsrelolo_h,
       urhsimhihi_h,urhsimlohi_h,urhsimhilo_h,urhsimlolo_h,
       urhsrehihi_d,urhsrelohi_d,urhsrehilo_d,urhsrelolo_d,
       urhsimhihi_d,urhsimlohi_d,urhsimhilo_d,urhsimlolo_d,
       Qrehihi_h,Qrelohi_h,Qrehilo_h,Qrelolo_h,
       Qimhihi_h,Qimlohi_h,Qimhilo_h,Qimlolo_h,
       Qrehihi_d,Qrelohi_d,Qrehilo_d,Qrelolo_d,
       Qimhihi_d,Qimlohi_d,Qimhilo_d,Qimlolo_d,
       Rrehihi_h,Rrelohi_h,Rrehilo_h,Rrelolo_h,
       Rimhihi_h,Rimlohi_h,Rimhilo_h,Rimlolo_h,
       Rrehihi_d,Rrelohi_d,Rrehilo_d,Rrelolo_d,
       Rimhihi_d,Rimlohi_d,Rimhilo_d,Rimlolo_d,
       solrehihi_h,solrelohi_h,solrehilo_h,solrelolo_h,
       solimhihi_h,solimlohi_h,solimhilo_h,solimlolo_h,
       solrehihi_d,solrelohi_d,solrehilo_d,solrelolo_d,
       solimhihi_d,solimlohi_d,solimhilo_d,solimlolo_d);

   if(vrblvl > 0) cout << "Setting up the test solution ..." << endl;

   double **testsolrehihi = new double*[dim];
   double **testsolrelohi = new double*[dim];
   double **testsolrehilo = new double*[dim];
   double **testsolrelolo = new double*[dim];
   double **testsolimhihi = new double*[dim];
   double **testsolimlohi = new double*[dim];
   double **testsolimhilo = new double*[dim];
   double **testsolimlolo = new double*[dim];
   double **mbrhsrehihi = new double*[dim];
   double **mbrhsrelohi = new double*[dim];
   double **mbrhsrehilo = new double*[dim];
   double **mbrhsrelolo = new double*[dim];
   double **mbrhsimhihi = new double*[dim];
   double **mbrhsimlohi = new double*[dim];
   double **mbrhsimhilo = new double*[dim];
   double **mbrhsimlolo = new double*[dim];

   cmplx4_row_setup(dim,deg,nbr,nvr,idx,
       cstrehihi,cstrelohi,cstrehilo,cstrelolo,
       cstimhihi,cstimlohi,cstimhilo,cstimlolo,
       cffrehihi,cffrelohi,cffrehilo,cffrelolo,
       cffimhihi,cffimlohi,cffimhilo,cffimlolo,
       testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
       testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
       inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
       inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
       inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
       inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d,
       outputrehihi_h,outputrelohi_h,outputrehilo_h,outputrelolo_h,
       outputimhihi_h,outputimlohi_h,outputimhilo_h,outputimlolo_h,
       outputrehihi_d,outputrelohi_d,outputrehilo_d,outputrelolo_d,
       outputimhihi_d,outputimlohi_d,outputimhilo_d,outputimlolo_d,mode,vrblvl);

   double *workvecrehihi = new double[dim];
   double *workvecrelohi = new double[dim];
   double *workvecrehilo = new double[dim];
   double *workvecrelolo = new double[dim];
   double *workvecimhihi = new double[dim];
   double *workvecimlohi = new double[dim];
   double *workvecimhilo = new double[dim];
   double *workvecimlolo = new double[dim];

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

   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
   bool zeroQ_h = false;
   bool zeroQ_d = false;
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

      cmplx4_row_newton_qrstep
         (szt,nbt,dim,wrkdeg,&tailidx_h,&tailidx_d,nbr,nvr,idx,
          cstrehihi,cstrelohi,cstrehilo,cstrelolo,
          cstimhihi,cstimlohi,cstimhilo,cstimlolo,
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,dpr,
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
          workvecrehihi,workvecrelohi,workvecrehilo,workvecrelolo,
          workvecimhihi,workvecimlohi,workvecimhilo,workvecimlolo,
          resvecrehihi,resvecrelohi,resvecrehilo,resvecrelolo,
          resvecimhihi,resvecimlohi,resvecimhilo,resvecimlolo,
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

   cmplx4_error_testsol(dim,deg,mode,
      testsolrehihi,testsolrelohi,testsolrehilo,testsolrelolo,
      testsolimhihi,testsolimlohi,testsolimhilo,testsolimlolo,
      inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
      inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
      inputrehihi_d,inputrelohi_d,inputrehilo_d,inputrelolo_d,
      inputimhihi_d,inputimlohi_d,inputimhilo_d,inputimlolo_d);

   write_newton_times
      (stepcnt,walltimesec,totcnvlapsedms,totqrlapsedms,totqtblapsedms,
       totbslapsedms,totupdlapsedms,totreslapsedms);

   return 0;
}
