// The file dbl2_newton_method.cpp defines the functions with prototypes in
// the file dbl2_newton_method.h.

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
#include "double_double_functions.h"
#include "dbl2_convolutions_host.h"
#include "dbl2_monomials_host.h"
#include "dbl2_factorizations.h"
#include "dbl2_monomial_systems.h"
#include "dbl2_bals_host.h"
#include "dbl2_bals_kernels.h"
#include "dbl2_tail_kernels.h"
#include "dbl2_systems_host.h"
#include "dbl2_systems_kernels.h"
#include "dbl_bals_flopcounts.h"
#include "dbl2_newton_testers.h"

using namespace std;

void dbl2_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbhi, double **mblo, double dpr,
   double ***cffhi, double ***cfflo, double **acchi, double **acclo,
   double **inputhi_h, double **inputlo_h,
   double **inputhi_d, double **inputlo_d,
   double ***outputhi_h, double ***outputlo_h,
   double ***outputhi_d, double ***outputlo_d,
   double **funvalhi_h, double **funvallo_h,
   double **funvalhi_d, double **funvallo_d,
   double ***jacvalhi_h, double ***jacvallo_h,
   double ***jacvalhi_d, double ***jacvallo_d,
   double **rhshi_h, double **rhslo_h, double **rhshi_d, double **rhslo_d,
   double **urhshi_h, double **urhslo_h, double **urhshi_d, double **urhslo_d,
   double **solhi_h, double **sollo_h, double **solhi_d, double **sollo_d,
   double **Qhi_h, double **Qlo_h, double **Qhi_d, double **Qlo_d,
   double **Rhi_h, double **Rlo_h, double **Rhi_d, double **Rlo_d,
   double *workvechi, double *workveclo,
   double **resvechi, double **resveclo, double *resmaxhi, double *resmaxlo,
   bool *noqr_h, bool *noqr_d,
   int *upidx_h, int *bsidx_h, int *upidx_d, int *bsidx_d,
   double *totcnvlapsedms, double *totqrlapsedms, double *totqtblapsedms,
   double *totbslapsedms, double *totupdlapsedms, double *totreslapsedms,
   int vrblvl, int mode )
{
   const int degp1 = deg+1;

   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling CPU_dbl2_evaluate_monomials ..." << endl;
 
      if(nbrcol == 1)
      {
         dbl2_unit_series_vector(dim,deg,cffhi[0],cfflo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_dbl2_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhi[0],cfflo[0],acchi[0],acclo[0],
             inputhi_h,inputlo_h,outputhi_h,outputlo_h,vrblvl);
      }
      else
         CPU_dbl2_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,cffhi,cfflo,acchi,acclo,
             inputhi_h,inputlo_h,funvalhi_h,funvallo_h,
             jacvalhi_h,jacvallo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_dbl2_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         dbl2_unit_series_vector(dim,deg,cffhi[0],cfflo[0]);
         // reset coefficients

         GPU_dbl2_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffhi[0],cfflo[0],acchi[0],acclo[0],
             inputhi_d,inputlo_d,outputhi_d,outputlo_d,totcnvlapsedms,vrblvl);
      }
      else
         GPU_dbl2_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cffhi,cfflo,
             inputhi_d,inputlo_d,outputhi_d,outputlo_d,
             funvalhi_d,funvallo_d,jacvalhi_d,jacvallo_d,
             totcnvlapsedms,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;
      double errsum = 0.0;

      errsum = dbl2_error3sum
                  (dim,dim+1,degp1,outputhi_h,outputlo_h,
                                   outputhi_d,outputlo_d,"output",vrblvl);
      // first dim is the number of monomials,
      // dim+1 is number of variables for each derivative,
      // plus the last component with the function value

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum << endl;
   }
   if(vrblvl > 0) cout << "initializing the Jacobian ..." << endl;
   cout << flush;

   if((mode == 1) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl2_define_rhs
            (dim,degp1,mbhi,mblo,funvalhi_h,funvallo_h,
             rhshi_h,rhslo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  cout << "  i = " << i
                       << "  j = " << j
                       << "  k = " << k << endl << flush;

                  jacvalhi_h[i][j][k] = 0.0; jacvallo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbhi,mblo,dpr,outputhi_h,outputlo_h,
             funvalhi_h,funvallo_h,rhshi_h,rhslo_h,
             jacvalhi_h,jacvallo_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         dbl2_define_rhs
            (dim,degp1,mbhi,mblo,funvalhi_d,funvallo_d,
             rhshi_d,rhslo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalhi_d[i][j][k] = 0.0; jacvallo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         dbl2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbhi,mblo,dpr,outputhi_d,outputlo_d,
             funvalhi_d,funvallo_d,rhshi_d,rhslo_d,jacvalhi_d,jacvallo_d,
             vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      cout << scientific << setprecision(3);
      cout << "comparing CPU with GPU function values ... " << endl;
      double errsum = 0.0;
      errsum = dbl2_error2sum
                  (dim,degp1,funvalhi_h,funvallo_h,
                             funvalhi_d,funvallo_d,"funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = dbl2_error3sum
                  (degp1,dim,dim,jacvalhi_h,jacvallo_h,
                                 jacvalhi_d,jacvallo_d,"jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = dbl2_error2sum
                  (degp1,dim,rhshi_h,rhslo_h,rhshi_d,rhslo_d,"rhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
   if((mode == 1) || (mode == 2))
   {
      if(vrblvl > 0) cout << "saving the original rhs ..." << endl;

      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhshi_h[i][j] = rhshi_h[i][j]; urhslo_h[i][j] = rhslo_h[i][j];
         }

      if(vrblvl > 0) cout << "initializing the solution ..." << endl;

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhi_h[i][j] = 0.0; sollo_h[i][j] = 0.0;
         }
      if(vrblvl > 0) cout << "calling CPU_dbl2_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_dbl2_qrbs_solve
         (dim,degp1,oldtail,jacvalhi_h,jacvallo_h,
          urhshi_h,urhslo_h,solhi_h,sollo_h,
          Qhi_h,Qlo_h,Rhi_h,Rlo_h,workvechi,workveclo,
          noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;

      if(vrblvl > 0)
      {
         cout << "calling CPU_dbl2_linear_residue ..." << endl;

         CPU_dbl2_linear_residue
            (dim,degp1,*tailidx_h-1,jacvalhi_h,jacvallo_h,rhshi_h,rhslo_h,
             solhi_h,sollo_h,resvechi,resveclo,resmaxhi,resmaxlo,vrblvl);
   
         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi << endl;
      }
      dbl2_update_series
         (dim,degp1,*tailidx_h-1,inputhi_h,inputlo_h,solhi_h,sollo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0) cout << "saving the original rhs ..." << endl;

      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhshi_d[i][j] = rhshi_d[i][j]; urhslo_d[i][j] = rhslo_d[i][j];
         }

      if(vrblvl > 0) cout << "initializing the solution ..." << endl;

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solhi_d[i][j] = 0.0; sollo_d[i][j] = 0.0;
         }

      if(vrblvl > 0) cout << "calling GPU_dbl2_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_dbl2_bals_solve
         (dim,degp1,szt,nbt,oldtail,jacvalhi_d,jacvallo_d,
          Qhi_d,Qlo_d,Rhi_d,Rlo_d,urhshi_d,urhslo_d,solhi_d,sollo_d,noqr_d,
          upidx_d,bsidx_d,&newtail,totqrlapsedms,totqtblapsedms,
          totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_dbl2_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_dbl2_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,jacvalhi_d,jacvallo_d,
             rhshi_d,rhslo_d,solhi_d,sollo_d,resvechi,resveclo,
             resmaxhi,resmaxlo,&elapsedms,&addcnt,&mulcnt,vrblvl);
   
         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;

         *totreslapsedms += elapsedms;
      }
      dbl2_update_series
         (dim,degp1,*tailidx_d-1,inputhi_d,inputlo_d,solhi_d,sollo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << scientific << setprecision(3);
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      errsum = dbl2_error2sum(dim,dim,Qhi_h,Qlo_h,Qhi_d,Qlo_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = dbl2_error2sum(dim,dim,Rhi_h,Rlo_h,Rhi_d,Rlo_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = dbl2_error2sum(degp1,dim,urhshi_h,urhslo_h,
                                        urhshi_d,urhslo_d,"urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = dbl2_error2sum(degp1,dim,solhi_h,sollo_h,
                                        solhi_d,sollo_d,"sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = dbl2_error2sum(dim,degp1,inputhi_h,inputlo_h,
                                        inputhi_d,inputlo_d,"input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl2_real_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputhi_h = new double*[dim];
   double **inputlo_h = new double*[dim];
   double **inputhi_d = new double*[dim];
   double **inputlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      inputhi_h[i] = new double[degp1];
      inputlo_h[i] = new double[degp1];
      inputhi_d[i] = new double[degp1];
      inputlo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **acchi = new double*[dim+1]; // accumulated power series
   double **acclo = new double*[dim+1]; // in one column
   for(int i=0; i<=dim; i++)
   {
      acchi[i] = new double[degp1];
      acclo[i] = new double[degp1];
   }
   double ***cffhi = new double**[nbrcol]; // coefficients of monomials
   double ***cfflo = new double**[nbrcol];
   for(int i=0; i<nbrcol; i++)
   {
      cffhi[i] = new double*[dim];
      cfflo[i] = new double*[dim];
      for(int j=0; j<dim; j++)
      {
         cffhi[i][j] = new double[degp1];
         cfflo[i][j] = new double[degp1];
      }
   }
   double ***outputhi_h;
   double ***outputlo_h;
   double ***outputhi_d;
   double ***outputlo_d;

   if((mode == 1) || (mode == 2))
   {
      outputhi_h = new double**[dim];
      outputlo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhi_h[i] = new double*[dim+1];
         outputlo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhi_h[i][j] = new double[degp1];
            outputlo_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputhi_d = new double**[dim];
      outputlo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputhi_d[i] = new double*[dim+1];
         outputlo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputhi_d[i][j] = new double[degp1];
            outputlo_d[i][j] = new double[degp1];
         }
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalhi_h;
   double **funvallo_h;
   double **funvalhi_d;
   double **funvallo_d;

   if((mode == 1) || (mode == 2))
   {
      funvalhi_h = new double*[dim];
      funvallo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhi_h[i] = new double[degp1];
         funvallo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalhi_d = new double*[dim];
      funvallo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalhi_d[i] = new double[degp1];
         funvallo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h;
   double ***jacvallo_h;
   double ***jacvalhi_d;
   double ***jacvallo_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalhi_h = new double**[degp1];
      jacvallo_h = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhi_h[i] = new double*[dim];
         jacvallo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhi_h[i][j] = new double[dim];
            jacvallo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      jacvalhi_d = new double**[degp1];
      jacvallo_d = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalhi_d[i] = new double*[dim];
         jacvallo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalhi_d[i][j] = new double[dim];
            jacvallo_d[i][j] = new double[dim];
         }
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.

   double **solhi_h;
   double **sollo_h;
   double **solhi_d;
   double **sollo_d;

   if((mode == 1) || (mode == 2))
   {
      solhi_h = new double*[degp1];
      sollo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhi_h[i] = new double[dim];
         sollo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solhi_d = new double*[degp1];
      sollo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         solhi_d[i] = new double[dim];
         sollo_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshi_h;
   double **rhslo_h;
   double **rhshi_d;
   double **rhslo_d;

   if((mode == 1) || (mode == 2))
   {
      rhshi_h = new double*[degp1];
      rhslo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshi_h[i] = new double[dim];
         rhslo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhshi_d = new double*[degp1];
      rhslo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhshi_d[i] = new double[dim];
         rhslo_d[i] = new double[dim];
      }
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshi_h;
   double **urhslo_h;
   double **urhshi_d;
   double **urhslo_d;

   if((mode == 1) || (mode == 2))
   {
      urhshi_h = new double*[degp1];
      urhslo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshi_h[i] = new double[dim];
         urhslo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhshi_d = new double*[degp1];
      urhslo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhshi_d[i] = new double[dim];
         urhslo_d[i] = new double[dim];
      }
   }
   double *workvechi = new double[dim];
   double *workveclo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **workrhshi = new double*[degp1];
   double **workrhslo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      workrhshi[i] = new double[dim];
      workrhslo[i] = new double[dim];
   }
   double **resvechi = new double*[degp1];
   double **resveclo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvechi[i] = new double[dim];
      resveclo[i] = new double[dim];
   }
   double resmaxhi,resmaxlo;

   double **Qhi_h;
   double **Qlo_h;
   double **Qhi_d;
   double **Qlo_d;

   if((mode == 1) || (mode == 2))
   {
      Qhi_h = new double*[dim];
      Qlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhi_h[i] = new double[dim];
         Qlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qhi_d = new double*[dim];
      Qlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qhi_d[i] = new double[dim];
         Qlo_d[i] = new double[dim];
      }
   }
   double **Rhi_h;
   double **Rlo_h;
   double **Rhi_d;
   double **Rlo_d;

   if((mode == 1) || (mode == 2))
   {
      Rhi_h = new double*[dim];
      Rlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhi_h[i] = new double[dim];
         Rlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rhi_d = new double*[dim];
      Rlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rhi_d[i] = new double[dim];
         Rlo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **solhi = new double*[dim];
   double **sollo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solhi[i] = new double[degp1];
      sollo[i] = new double[degp1];
   }
   make_real2_exponentials(dim,deg,solhi,sollo);
   if(nbrcol != 1) // generate coefficients for the columns
      make_real2_coefficients(nbrcol,dim,cffhi,cfflo);

   // compute the right hand sides via evaluation

   double **mbrhshi = new double*[dim];
   double **mbrhslo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhshi[i] = new double[degp1];
      mbrhslo[i] = new double[degp1];

      mbrhshi[i][0] = 1.0;     // initialize product to one
      mbrhslo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhshi[i][k] = 0.0;
         mbrhslo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_real2_monomials(dim,deg,rowsA,solhi,sollo,mbrhshi,mbrhslo);
   else
      evaluate_real2_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffhi,cfflo,solhi,sollo,
          mbrhshi,mbrhslo,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhshi[i][j] << "  " << mbrhslo[i][j] << endl;
   }
   double *start0hi = new double[dim];
   double *start0lo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0hi[i] = solhi[i][0];
      start0lo[i] = sollo[i][0];
   }
   real2_start_series_vector(dim,deg,start0hi,start0lo,inputhi_h,inputlo_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhi_d[i][j] = inputhi_h[i][j];
         inputlo_d[i][j] = inputlo_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : " << inputhi_h[i][0] << "  "
                            << inputlo_h[i][0] << endl;
   }
   if(vrblvl > 0) cout << scientific << setprecision(16);

   int upidx_h = 0;
   int bsidx_h = 0;
   int upidx_d = 0;
   int bsidx_d = 0;
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

      dbl2_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,
          mbrhshi,mbrhslo,dpr,cffhi,cfflo,acchi,acclo,
          inputhi_h,inputlo_h,inputhi_d,inputlo_d,
          outputhi_h,outputlo_h,outputhi_d,outputlo_d,
          funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,
          jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,
          rhshi_h,rhslo_h,rhshi_d,rhslo_d,urhshi_h,urhslo_h,urhshi_d,urhslo_d,
          solhi_h,sollo_h,solhi_d,sollo_d,
          Qhi_h,Qlo_h,Qhi_d,Qlo_d,Rhi_h,Rlo_h,Rhi_d,Rlo_d,
          workvechi,workveclo, resvechi,resveclo,&resmaxhi,&resmaxlo,
          &noqr_h,&noqr_d,&upidx_h,&bsidx_h,&upidx_d,&bsidx_d,
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

   double errsum = 0.0;

   cout << scientific << setprecision(16); // just in case vrblvl == 0
   cout << "The solution series : " << endl;
   for(int j=0; j<degp1; j++)
   {
      cout << "coefficient of degree " << j << " :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "sol[" << i << "][" << j << "] : "
                        << solhi[i][j] << "  "
                        << sollo[i][j] << endl;
         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputhi_d[i][j] << "  "
                           << inputlo_d[i][j] << endl;
            errsum += abs(solhi[i][j] - inputhi_d[i][j])
                    + abs(sollo[i][j] - inputlo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputhi_h[i][j] << "  "
                           << inputlo_h[i][j] << endl;
            errsum += abs(solhi[i][j] - inputhi_h[i][j])
                    + abs(sollo[i][j] - inputlo_h[i][j]);
         }
      }
   }
   cout << "error : " << errsum << endl;

   cout << "Wall clock time on all " << stepcnt << " Newton steps : ";
   cout << fixed << setprecision(3) 
        << walltimesec << " seconds." << endl;
   cout << "     Time spent by all convolution kernels : "
        << totcnvlapsedms << " milliseconds." << endl;
   cout << "  Time spent by all Householder QR kernels : "
        << totqrlapsedms << " milliseconds." << endl;
   cout << "     Time spent by all Q times rhs kernels : "
        << totqtblapsedms << " milliseconds." << endl;
   cout << "Time spent by all backsubstitution kernels : "
        << totbslapsedms << " milliseconds." << endl;
   cout << "          Time spent by all update kernels : "
        << totupdlapsedms << " milliseconds." << endl;
   cout << "        Time spent by all residual kernels : "
        << totreslapsedms << " milliseconds." << endl;

   double totkerneltime = totcnvlapsedms + totqrlapsedms + totqtblapsedms
                        + totbslapsedms + totupdlapsedms + totreslapsedms;

   cout << "           Total time spent by all kernels : "
        << totkerneltime << " milliseconds." << endl;

   return 0;
}
