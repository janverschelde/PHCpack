// The file cmplx2_newton_method.cpp defines the functions with prototypes in
// the file cmplx2_newton_method.h.

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

void cmplx2_newton_qrstep
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int *tailidx_h, int *tailidx_d,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac,
   double **mbrehi, double **mbrelo, double **mbimhi, double **mbimlo,
   double dpr,
   double ***cffrehi, double ***cffrelo, double ***cffimhi, double ***cffimlo,
   double **accrehi, double **accrelo, double **accimhi, double **accimlo,
   double **inputrehi_h, double **inputrelo_h,
   double **inputimhi_h, double **inputimlo_h,
   double **inputrehi_d, double **inputrelo_d,
   double **inputimhi_d, double **inputimlo_d,
   double ***outputrehi_h, double ***outputrelo_h,
   double ***outputimhi_h, double ***outputimlo_h,
   double ***outputrehi_d, double ***outputrelo_d,
   double ***outputimhi_d, double ***outputimlo_d,
   double **funvalrehi_h, double **funvalrelo_h,
   double **funvalimhi_h, double **funvalimlo_h,
   double **funvalrehi_d, double **funvalrelo_d,
   double **funvalimhi_d, double **funvalimlo_d,
   double ***jacvalrehi_h, double ***jacvalrelo_h,
   double ***jacvalimhi_h, double ***jacvalimlo_h,
   double ***jacvalrehi_d, double ***jacvalrelo_d,
   double ***jacvalimhi_d, double ***jacvalimlo_d,
   double **rhsrehi_h, double **rhsrelo_h,
   double **rhsimhi_h, double **rhsimlo_h,
   double **rhsrehi_d, double **rhsrelo_d, 
   double **rhsimhi_d, double **rhsimlo_d,
   double **urhsrehi_h, double **urhsrelo_h,
   double **urhsimhi_h, double **urhsimlo_h,
   double **urhsrehi_d, double **urhsrelo_d,
   double **urhsimhi_d, double **urhsimlo_d,
   double **solrehi_h, double **solrelo_h,
   double **solimhi_h, double **solimlo_h, 
   double **solrehi_d, double **solrelo_d, 
   double **solimhi_d, double **solimlo_d,
   double **Qrehi_h, double **Qrelo_h, double **Qimhi_h, double **Qimlo_h,
   double **Qrehi_d, double **Qrelo_d, double **Qimhi_d, double **Qimlo_d, 
   double **Rrehi_h, double **Rrelo_h, double **Rimhi_h, double **Rimlo_h, 
   double **Rrehi_d, double **Rrelo_d, double **Rimhi_d, double **Rimlo_d,
   double *workvecrehi, double *workvecrelo,
   double *workvecimhi, double *workvecimlo,
   double **resvecrehi, double **resvecrelo,
   double **resvecimhi, double **resvecimlo,
   double *resmaxhi, double *resmaxlo,
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
         cout << "calling CPU_cmplx2_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx2_unit_series_vector
            (dim,deg,cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0]);
         // The series coefficients accumulate common factors,
         // initially the coefficients are set to one.

         CPU_cmplx2_evaluate_monomials
            (dim,deg,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0],
             accrehi[0],accrelo[0],accimhi[0],accimlo[0],
             inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
             outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,vrblvl);
      }
      else
         CPU_cmplx2_evaluate_columns
            (dim,deg,nbrcol,nvr,idx,
             cffrehi,cffrelo,cffimhi,cffimlo,
             accrehi,accrelo,accimhi,accimlo,
             inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
             funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      if(vrblvl > 0)
         cout << "calling GPU_cmplx2_evaluate_monomials ..." << endl;

      if(nbrcol == 1)
      {
         cmplx2_unit_series_vector
            (dim,deg,cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0]);
         // reset the coefficients

         GPU_cmplx2_evaluate_monomials
            (dim,deg,szt,nbt,nvr[0],idx[0],exp,nbrfac,expfac,
             cffrehi[0],cffrelo[0],cffimhi[0],cffimlo[0],
             accrehi[0],accrelo[0],accimhi[0],accimlo[0],
             inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             totcnvlapsedms,vrblvl);
      }
      else
         GPU_cmplx2_evaluate_columns
            (dim,deg,nbrcol,szt,nbt,nvr,idx,cffrehi,cffrelo,cffimhi,cffimlo,
             inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
             totcnvlapsedms,vrblvl);

   }
   if((vrblvl > 0) && (mode == 2) && (nbrcol == 1))
   {
      cout << "comparing CPU with GPU evaluations ... " << endl;

      double errsum = cmplx2_error3sum
                         (dim,dim+1,degp1,
                          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
                          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
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
         cmplx2_define_rhs
            (dim,degp1,mbrehi,mbrelo,mbimhi,mbimlo,
             funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalrehi_h[i][j][k] = 0.0; jacvalimhi_h[i][j][k] = 0.0;
                  jacvalrelo_h[i][j][k] = 0.0; jacvalimlo_h[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         cmplx2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbrehi,mbrelo,mbimhi,mbimlo,dpr,
             outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
             funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,vrblvl);
      }
   }
   if((mode == 0) || (mode == 2))
   {
      if(nbrcol != 1)
         cmplx2_define_rhs
            (dim,degp1,mbrehi,mbrelo,mbimhi,mbimlo,
             funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,vrblvl);
      else
      {
         for(int i=0; i<degp1; i++) // initialize the Jacobian to zero
            for(int j=0; j<dim; j++) 
               for(int k=0; k<dim; k++)
               {
                  jacvalrehi_d[i][j][k] = 0.0; jacvalimhi_d[i][j][k] = 0.0;
                  jacvalrelo_d[i][j][k] = 0.0; jacvalimlo_d[i][j][k] = 0.0;
               }

         if(vrblvl > 0) cout << "linearizing the output ..." << endl;
         cmplx2_linearize_evaldiff_output
            (dim,degp1,nvr[0],idx[0],mbrehi,mbrelo,mbimhi,mbimlo,dpr,
             outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
             funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,vrblvl);
      }
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << scientific << setprecision(3);
      cout << "comparing CPU with GPU function values ... " << endl;
      errsum = cmplx2_error2sum(dim,degp1,
                  funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
                  funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
                  "funval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU Jacobians ... " << endl;
      errsum = cmplx2_error3sum(degp1,dim,dim,
                  jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
                  jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
                  "jacval",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU right hand sides ... " << endl;
      errsum = cmplx2_error2sum(degp1,dim,
                  rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
                  rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,"rhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }

   if((mode == 1) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhsrehi_h[i][j] = rhsrehi_h[i][j];
            urhsimhi_h[i][j] = rhsimhi_h[i][j];
            urhsrelo_h[i][j] = rhsrelo_h[i][j];
            urhsimlo_h[i][j] = rhsimlo_h[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solrehi_h[i][j] = 0.0; solimhi_h[i][j] = 0.0;
            solrelo_h[i][j] = 0.0; solimlo_h[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling CPU_cmplx2_qrbs_solve ..." << endl;

      int oldtail = *tailidx_h;
      int newtail = oldtail;

      CPU_cmplx2_qrbs_solve
         (dim,degp1,oldtail,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
          urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,
          noqr_h,upidx_h,bsidx_h,&newtail,vrblvl);

      *tailidx_h = newtail;
 
      if(vrblvl > 0)
      {
         cout << "calling CPU_cmplx2_linear_residue ..." << endl;

         CPU_cmplx2_linear_residue
            (dim,degp1,*tailidx_h-1,
             jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
             rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,solrehi_h,
             solrelo_h,solimhi_h,solimlo_h,
             resvecrehi,resvecrelo,resvecimhi,resvecimlo,
             resmaxhi,resmaxlo,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi << endl;
      }
      cmplx2_update_series
         (dim,degp1,*tailidx_h-1,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,vrblvl);
   }
   if((mode == 0) || (mode == 2))
   {
      for(int i=0; i<degp1; i++) // save original rhs for residual
         for(int j=0; j<dim; j++)
         {
            urhsrehi_d[i][j] = rhsrehi_d[i][j];
            urhsimhi_d[i][j] = rhsimhi_d[i][j];
            urhsrelo_d[i][j] = rhsrelo_d[i][j];
            urhsimlo_d[i][j] = rhsimlo_d[i][j];
         }

      for(int i=0; i<degp1; i++) // initialize the solution to zero
         for(int j=0; j<dim; j++)
         {
            solrehi_d[i][j] = 0.0; solimhi_d[i][j] = 0.0;
            solrelo_d[i][j] = 0.0; solimlo_d[i][j] = 0.0;
         }

      if(vrblvl > 0)
         cout << "calling GPU_cmplx2_bals_solve ..." << endl;

      int oldtail = *tailidx_d;
      int newtail = oldtail;

      GPU_cmplx2_bals_solve
         (dim,degp1,szt,nbt,oldtail,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
          Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,noqr_d,
          upidx_d,bsidx_d,&newtail,totqrlapsedms,totqtblapsedms,
          totbslapsedms,totupdlapsedms,vrblvl);

      *tailidx_d = newtail;

      if(vrblvl > 0)
      {
         cout << "calling GPU_cmplx2_linear_residue ..." << endl;

         double elapsedms = 0.0;
         long long int addcnt = 0;
         long long int mulcnt = 0;

         GPU_cmplx2_linear_residue
            (dim,degp1,szt,nbt,*tailidx_d-1,
             jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
             rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
             solrehi_d,solrelo_d,solimhi_d,solimlo_d,
             resvecrehi,resvecrelo,resvecimhi,resvecimlo,
             resmaxhi,resmaxlo,&elapsedms,&addcnt,&mulcnt,vrblvl);

         cout << scientific << setprecision(3)
              << "maximum residual : " << *resmaxhi;
         cout << fixed << setprecision(3)
              << "  total kernel time : " << elapsedms
              << " milliseconds" << endl;

         *totreslapsedms += elapsedms;
      }
      cmplx2_update_series
         (dim,degp1,*tailidx_d-1,
          inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,vrblvl);
   }
   if((vrblvl > 0) && (mode == 2))
   {
      double errsum = 0.0;
      cout << scientific << setprecision(3);
      cout << "comparing CPU with GPU matrices Q ... " << endl;
      errsum = cmplx2_error2sum(dim,dim,
                  Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,
                  Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,"Q",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU matrices R ... " << endl;
      errsum = cmplx2_error2sum(dim,dim,
                  Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,
                  Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,"R",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU updated rhs ... " << endl;
      errsum = cmplx2_error2sum(degp1,dim,
                  urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
                  urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,"urhs",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU update to solutions ... " << endl;
      errsum = cmplx2_error2sum(degp1,dim,
                  solrehi_h,solrelo_h,solimhi_h,solimlo_h,
                  solrehi_d,solrelo_d,solimhi_d,solimlo_d,"sol",vrblvl);
      cout << "sum of errors : " << errsum << endl;
      cout << "comparing CPU with GPU series ... " << endl;
      errsum = cmplx2_error2sum(dim,degp1,
                  inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
                  inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
                  "input",vrblvl);
      cout << "sum of errors : " << errsum << endl;
   }
}

int test_dbl2_complex_newton
 ( int szt, int nbt, int dim, int deg, int nbrcol,
   int **nvr, int ***idx, int **exp, int *nbrfac, int **expfac, int **rowsA,
   double dpr, int nbsteps, int mode, int vrblvl )
{
/*
 * 1. allocating input and output space for evaluation and differentiation
 */
   const int degp1 = deg+1;
   double **inputrehi_h = new double*[dim];
   double **inputrelo_h = new double*[dim];
   double **inputimhi_h = new double*[dim];
   double **inputimlo_h = new double*[dim];
   double **inputrehi_d = new double*[dim];
   double **inputrelo_d = new double*[dim];
   double **inputimhi_d = new double*[dim];
   double **inputimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputrehi_h[i] = new double[degp1];
       inputrelo_h[i] = new double[degp1];
       inputimhi_h[i] = new double[degp1];
       inputimlo_h[i] = new double[degp1];
       inputrehi_d[i] = new double[degp1];
       inputrelo_d[i] = new double[degp1];
       inputimhi_d[i] = new double[degp1];
       inputimlo_d[i] = new double[degp1];
   }
   // allocate memory for coefficients and the output
   double **accrehi = new double*[dim+1]; // accumulated power series
   double **accrelo = new double*[dim+1]; // in one column
   double **accimhi = new double*[dim+1];
   double **accimlo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehi[i] = new double[degp1];
      accrelo[i] = new double[degp1];
      accimhi[i] = new double[degp1];
      accimlo[i] = new double[degp1];
   }
   double ***cffrehi = new double**[nbrcol]; // coefficients of monomials
   double ***cffrelo = new double**[nbrcol];
   double ***cffimhi = new double**[nbrcol]; 
   double ***cffimlo = new double**[nbrcol]; 

   for(int i=0; i<nbrcol; i++)
   {
      cffrehi[i] = new double*[dim];
      cffrelo[i] = new double*[dim];
      cffimhi[i] = new double*[dim];
      cffimlo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffrehi[i][j] = new double[degp1];
         cffrelo[i][j] = new double[degp1];
         cffimhi[i][j] = new double[degp1];
         cffimlo[i][j] = new double[degp1];
      }
   }
   double ***outputrehi_h;
   double ***outputrelo_h;
   double ***outputimhi_h;
   double ***outputimlo_h;
   double ***outputrehi_d;
   double ***outputrelo_d;
   double ***outputimhi_d;
   double ***outputimlo_d;

   if((mode == 1) || (mode == 2))
   {
      outputrehi_h = new double**[dim];
      outputrelo_h = new double**[dim];
      outputimhi_h = new double**[dim];
      outputimlo_h = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputrehi_h[i] = new double*[dim+1];
         outputrelo_h[i] = new double*[dim+1];
         outputimhi_h[i] = new double*[dim+1];
         outputimlo_h[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputrehi_h[i][j] = new double[degp1];
            outputrelo_h[i][j] = new double[degp1];
            outputimhi_h[i][j] = new double[degp1];
            outputimlo_h[i][j] = new double[degp1];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      outputrehi_d = new double**[dim];
      outputrelo_d = new double**[dim];
      outputimhi_d = new double**[dim];
      outputimlo_d = new double**[dim];

      for(int i=0; i<dim; i++)
      {
         outputrehi_d[i] = new double*[dim+1];
         outputrelo_d[i] = new double*[dim+1];
         outputimhi_d[i] = new double*[dim+1];
         outputimlo_d[i] = new double*[dim+1];

         for(int j=0; j<=dim; j++)
         {
            outputrehi_d[i][j] = new double[degp1];
            outputrelo_d[i][j] = new double[degp1];
            outputimhi_d[i][j] = new double[degp1];
            outputimlo_d[i][j] = new double[degp1];
         }
      }
   }
   // The function values are power series truncated at degree deg.
   double **funvalrehi_h;
   double **funvalrelo_h;
   double **funvalimhi_h;
   double **funvalimlo_h;
   double **funvalrehi_d;
   double **funvalrelo_d;
   double **funvalimhi_d;
   double **funvalimlo_d;

   if((mode == 1) || (mode == 2))
   {
      funvalrehi_h = new double*[dim];
      funvalrelo_h = new double*[dim];
      funvalimhi_h = new double*[dim];
      funvalimlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalrehi_h[i] = new double[degp1];
         funvalrelo_h[i] = new double[degp1];
         funvalimhi_h[i] = new double[degp1];
         funvalimlo_h[i] = new double[degp1];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      funvalrehi_d = new double*[dim];
      funvalrelo_d = new double*[dim];
      funvalimhi_d = new double*[dim];
      funvalimlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         funvalrehi_d[i] = new double[degp1];
         funvalrelo_d[i] = new double[degp1];
         funvalimhi_d[i] = new double[degp1];
         funvalimlo_d[i] = new double[degp1];
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalrehi_h;
   double ***jacvalrelo_h;
   double ***jacvalimhi_h;
   double ***jacvalimlo_h;
   double ***jacvalrehi_d;
   double ***jacvalrelo_d;
   double ***jacvalimhi_d;
   double ***jacvalimlo_d;

   if((mode == 1) || (mode == 2))
   {
      jacvalrehi_h = new double**[degp1];
      jacvalrelo_h = new double**[degp1];
      jacvalimhi_h = new double**[degp1];
      jacvalimlo_h = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalrehi_h[i] = new double*[dim];
         jacvalrelo_h[i] = new double*[dim];
         jacvalimhi_h[i] = new double*[dim];
         jacvalimlo_h[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalrehi_h[i][j] = new double[dim];
            jacvalrelo_h[i][j] = new double[dim];
            jacvalimhi_h[i][j] = new double[dim];
            jacvalimlo_h[i][j] = new double[dim];
         }
      }
   }
   if((mode == 0) || (mode == 2))
   {
      jacvalrehi_d = new double**[degp1];
      jacvalrelo_d = new double**[degp1];
      jacvalimhi_d = new double**[degp1];
      jacvalimlo_d = new double**[degp1];

      for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
      {
         jacvalrehi_d[i] = new double*[dim];
         jacvalrelo_d[i] = new double*[dim];
         jacvalimhi_d[i] = new double*[dim];
         jacvalimlo_d[i] = new double*[dim];

         for(int j=0; j<dim; j++)
         {
            jacvalrehi_d[i][j] = new double[dim];
            jacvalrelo_d[i][j] = new double[dim];
            jacvalimhi_d[i][j] = new double[dim];
            jacvalimlo_d[i][j] = new double[dim];
         }
      }
   }
/*
 * 2. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solrehi_h;
   double **solrelo_h;
   double **solimhi_h;
   double **solimlo_h;
   double **solrehi_d;
   double **solrelo_d;
   double **solimhi_d;
   double **solimlo_d;

   if((mode == 1) || (mode == 2))
   {
      solrehi_h = new double*[degp1];
      solrelo_h = new double*[degp1];
      solimhi_h = new double*[degp1];
      solimlo_h = new double*[degp1];

      for(int i=0; i<degp1; i++) 
      {
         solrehi_h[i] = new double[dim];
         solrelo_h[i] = new double[dim];
         solimhi_h[i] = new double[dim];
         solimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      solrehi_d = new double*[degp1];
      solrelo_d = new double*[degp1];
      solimhi_d = new double*[degp1];
      solimlo_d = new double*[degp1];

      for(int i=0; i<degp1; i++) 
      {
         solrehi_d[i] = new double[dim];
         solrelo_d[i] = new double[dim];
         solimhi_d[i] = new double[dim];
         solimlo_d[i] = new double[dim];
      }
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhsrehi_h;
   double **rhsrelo_h;
   double **rhsimhi_h;
   double **rhsimlo_h;
   double **rhsrehi_d;
   double **rhsrelo_d;
   double **rhsimhi_d;
   double **rhsimlo_d;

   if((mode == 1) || (mode == 2))
   {
      rhsrehi_h = new double*[degp1];
      rhsrelo_h = new double*[degp1];
      rhsimhi_h = new double*[degp1];
      rhsimlo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsrehi_h[i] = new double[dim];
         rhsrelo_h[i] = new double[dim];
         rhsimhi_h[i] = new double[dim];
         rhsimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      rhsrehi_d = new double*[degp1];
      rhsrelo_d = new double*[degp1];
      rhsimhi_d = new double*[degp1];
      rhsimlo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         rhsrehi_d[i] = new double[dim];
         rhsrelo_d[i] = new double[dim];
         rhsimhi_d[i] = new double[dim];
         rhsimlo_d[i] = new double[dim];
      }
   }
   double *workvecrehi = new double[dim];
   double *workvecrelo = new double[dim];
   double *workvecimhi = new double[dim];
   double *workvecimlo = new double[dim];
   // Copy the rhs vector into work space for inplace solver.
   double **urhsrehi_h;
   double **urhsrelo_h;
   double **urhsimhi_h;
   double **urhsimlo_h;
   double **urhsrehi_d;
   double **urhsrelo_d;
   double **urhsimhi_d;
   double **urhsimlo_d;

   if((mode == 1) || (mode == 2))
   {
      urhsrehi_h = new double*[degp1];
      urhsrelo_h = new double*[degp1];
      urhsimhi_h = new double*[degp1];
      urhsimlo_h = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsrehi_h[i] = new double[dim];
         urhsrelo_h[i] = new double[dim];
         urhsimhi_h[i] = new double[dim];
         urhsimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      urhsrehi_d = new double*[degp1];
      urhsrelo_d = new double*[degp1];
      urhsimhi_d = new double*[degp1];
      urhsimlo_d = new double*[degp1];

      for(int i=0; i<degp1; i++)
      {
         urhsrehi_d[i] = new double[dim];
         urhsrelo_d[i] = new double[dim];
         urhsimhi_d[i] = new double[dim];
         urhsimlo_d[i] = new double[dim];
      }
   }
   double **resvecrehi = new double*[degp1];
   double **resvecrelo = new double*[degp1];
   double **resvecimhi = new double*[degp1];
   double **resvecimlo = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      resvecrehi[i] = new double[dim];
      resvecrelo[i] = new double[dim];
      resvecimhi[i] = new double[dim];
      resvecimlo[i] = new double[dim];
   }
   double resmaxhi,resmaxlo;

   double **Qrehi_h;
   double **Qrelo_h;
   double **Qimhi_h;
   double **Qimlo_h;
   double **Qrehi_d;
   double **Qrelo_d;
   double **Qimhi_d;
   double **Qimlo_d;

   if((mode == 1) || (mode == 2))
   {
      Qrehi_h = new double*[dim];
      Qrelo_h = new double*[dim];
      Qimhi_h = new double*[dim];
      Qimlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qrehi_h[i] = new double[dim];
         Qrelo_h[i] = new double[dim];
         Qimhi_h[i] = new double[dim];
         Qimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Qrehi_d = new double*[dim];
      Qrelo_d = new double*[dim];
      Qimhi_d = new double*[dim];
      Qimlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Qrehi_d[i] = new double[dim];
         Qrelo_d[i] = new double[dim];
         Qimhi_d[i] = new double[dim];
         Qimlo_d[i] = new double[dim];
      }
   }
   double **Rrehi_h;
   double **Rrelo_h;
   double **Rimhi_h;
   double **Rimlo_h;
   double **Rrehi_d;
   double **Rrelo_d;
   double **Rimhi_d;
   double **Rimlo_d;

   if((mode == 1) || (mode == 2))
   {
      Rrehi_h = new double*[dim];
      Rrelo_h = new double*[dim];
      Rimhi_h = new double*[dim];
      Rimlo_h = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rrehi_h[i] = new double[dim];
         Rrelo_h[i] = new double[dim];
         Rimhi_h[i] = new double[dim];
         Rimlo_h[i] = new double[dim];
      }
   }
   if((mode == 0) || (mode == 2))
   {
      Rrehi_d = new double*[dim];
      Rrelo_d = new double*[dim];
      Rimhi_d = new double*[dim];
      Rimlo_d = new double*[dim];

      for(int i=0; i<dim; i++)
      {
         Rrehi_d[i] = new double[dim];
         Rrelo_d[i] = new double[dim];
         Rimhi_d[i] = new double[dim];
         Rimlo_d[i] = new double[dim];
      }
   }
/*
 * 3. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the test solution and the start series.

   if(vrblvl > 0) cout << "setting up the test system ..." << endl;

   double **solrehi = new double*[dim];
   double **solrelo = new double*[dim];
   double **solimhi = new double*[dim];
   double **solimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solrehi[i] = new double[degp1];
      solrelo[i] = new double[degp1];
      solimhi[i] = new double[degp1];
      solimlo[i] = new double[degp1];
   }
   make_complex2_exponentials(dim,deg,solrehi,solrelo,solimhi,solimlo);
   if(nbrcol != 1) // generate coefficients for the columns
      make_complex2_coefficients(nbrcol,dim,cffrehi,cffrelo,cffimhi,cffimlo);

   // compute the right hand sides via evaluation

   double **mbrhsrehi = new double*[dim];
   double **mbrhsrelo = new double*[dim];
   double **mbrhsimhi = new double*[dim];
   double **mbrhsimlo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      mbrhsrehi[i] = new double[degp1];
      mbrhsrelo[i] = new double[degp1];
      mbrhsimhi[i] = new double[degp1];
      mbrhsimlo[i] = new double[degp1];

      mbrhsrehi[i][0] = 1.0;     // initialize product to one
      mbrhsrelo[i][0] = 0.0;
      mbrhsimhi[i][0] = 0.0;
      mbrhsimlo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         mbrhsrehi[i][k] = 0.0;
         mbrhsrelo[i][k] = 0.0;
         mbrhsimhi[i][k] = 0.0;
         mbrhsimlo[i][k] = 0.0;
      }
   }
   if(nbrcol == 1)
      evaluate_complex2_monomials
         (dim,deg,rowsA,solrehi,solrelo,solimhi,solimlo,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo);
   else
      evaluate_complex2_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,cffrehi,cffrelo,cffimhi,cffimlo,
          solrehi,solrelo,solimhi,solimlo,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,vrblvl);

   if(vrblvl > 1)
   {
      cout << "the right hand side series :" << endl;
      cout << scientific << setprecision(16);
      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
            cout << "rhs[" << i << "][" << j << "] : "
                 << mbrhsrehi[i][j] << "  " << mbrhsrelo[i][j] << endl
                 << "  "
                 << mbrhsimhi[i][j] << "  " << mbrhsimlo[i][j] << endl;
   }
   double *start0rehi = new double[dim];
   double *start0relo = new double[dim];
   double *start0imhi = new double[dim];
   double *start0imlo = new double[dim];

   for(int i=0; i<dim; i++)  // compute start vector
   {
      start0rehi[i] = solrehi[i][0];
      start0relo[i] = solrelo[i][0];
      start0imhi[i] = solimhi[i][0]; 
      start0imlo[i] = solimlo[i][0]; 
   }
   cmplx2_start_series_vector
      (dim,deg,start0rehi,start0relo,start0imhi,start0imlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h);

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehi_d[i][j] = inputrehi_h[i][j];
         inputrelo_d[i][j] = inputrelo_h[i][j];
         inputimhi_d[i][j] = inputimhi_h[i][j];
         inputimlo_d[i][j] = inputimlo_h[i][j];
      }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
         cout << i << " : "
              << inputrehi_h[i][0] << "  "
              << inputrelo_h[i][0] << endl << "  "
              << inputimhi_h[i][0] << "  "
              << inputimlo_h[i][0] << endl;
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

      cmplx2_newton_qrstep
         (szt,nbt,dim,wrkdeg,nbrcol,&tailidx_h,&tailidx_d,
          nvr,idx,exp,nbrfac,expfac,
          mbrhsrehi,mbrhsrelo,mbrhsimhi,mbrhsimlo,dpr,
          cffrehi,cffrelo,cffimhi,cffimlo,accrehi,accrelo,accimhi,accimlo,
          inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
          inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
          outputrehi_h,outputrelo_h,outputimhi_h,outputimlo_h,
          outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
          funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
          funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
          jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
          jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
          rhsrehi_h,rhsrelo_h,rhsimhi_h,rhsimlo_h,
          rhsrehi_d,rhsrelo_d,rhsimhi_d,rhsimlo_d,
          urhsrehi_h,urhsrelo_h,urhsimhi_h,urhsimlo_h,
          urhsrehi_d,urhsrelo_d,urhsimhi_d,urhsimlo_d,
          solrehi_h,solrelo_h,solimhi_h,solimlo_h,
          solrehi_d,solrelo_d,solimhi_d,solimlo_d,
          Qrehi_h,Qrelo_h,Qimhi_h,Qimlo_h,Qrehi_d,Qrelo_d,Qimhi_d,Qimlo_d,
          Rrehi_h,Rrelo_h,Rimhi_h,Rimlo_h,Rrehi_d,Rrelo_d,Rimhi_d,Rimlo_d,
          workvecrehi,workvecrelo,workvecimhi,workvecimlo,
          resvecrehi,resvecrelo,resvecimhi,resvecimlo,&resmaxhi,&resmaxlo,
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
                        << solrehi[i][j] << "  "
                        << solrelo[i][j] << endl << "  "
                        << solimhi[i][j] << "  "
                        << solimlo[i][j] << endl;
         if((mode == 0) || (mode == 2))
         {
            cout << "x_d[" << i << "][" << j << "] : "
                           << inputrehi_d[i][j] << "  "
                           << inputrelo_d[i][j] << endl << "  "
                           << inputimhi_d[i][j] << "  "
                           << inputimlo_d[i][j] << endl;
            errsum += abs(solrehi[i][j] - inputrehi_d[i][j])
                    + abs(solrelo[i][j] - inputrelo_d[i][j])
                    + abs(solimhi[i][j] - inputimhi_d[i][j])
                    + abs(solimlo[i][j] - inputimlo_d[i][j]);
         }
         if((mode == 1) || (mode == 2))
         {
            cout << "x_h[" << i << "][" << j << "] : "
                           << inputrehi_h[i][j] << "  "
                           << inputrelo_h[i][j] << endl << "  "
                           << inputimhi_h[i][j] << "  "
                           << inputimlo_h[i][j] << endl;
            errsum += abs(solrehi[i][j] - inputrehi_h[i][j])
                    + abs(solrelo[i][j] - inputrelo_h[i][j])
                    + abs(solimhi[i][j] - inputimhi_h[i][j])
                    + abs(solimlo[i][j] - inputimlo_h[i][j]);
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
