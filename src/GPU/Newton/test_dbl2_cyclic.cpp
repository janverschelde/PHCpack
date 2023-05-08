/* Tests evaluation and differentiation of the column representation of 
 * the cyclic n-roots system in double double precision.
 * In this column representation, the polynomial system is a sequence
 * of monomial systems, each column defines one monomial system.
 * Each column can be evaluated and differentiated independently,
 * just as each monomial in every column. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "dbl2_systems_host.h"
#include "dbl2_monomial_systems.h"
#include "dbl2_systems_kernels.h"
#include "dbl2_newton_testers.h"
#include "cyclic_columns.h"

using namespace std;

int dbl2_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns
 *   on real data in double double precision, both on CPU and GPU.
 *
 * ON ENTRY :
 *   nbt      number of tiles;
 *   szt      size of each tile;
 *   deg      degree of the series;
 *   dim      dimension of the problem;
 *   nvr      number of variables for each monomial;
 *   idx      index representation of the monomials;
 *   vrblvl   is the verbose level. */

int dbl2_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns,
 *   on complex data in double double precision, both on CPU and GPU.
 *
 * ON ENTRY :
 *   nbt      number of tiles;
 *   szt      size of each tile;
 *   deg      degree of the series;
 *   dim      dimension of the problem;
 *   nvr      number of variables for each monomial;
 *   idx      index representation of the monomials;
 *   vrblvl   is the verbose level. */

int main ( void )
{
   cout << "testing evaluation of cyclic n-roots columns ..." << endl;

   cout << "-> give the dimension : ";
   int dim; cin >> dim;
   cout << "-> give the degree : ";
   int deg; cin >> deg;
   // cout << "-> give the number of tiles : ";
   // int nbt; cin >> nbt;
   const int nbt = 1;
   const int szt = deg+1;
   // cout << "-> give the size of each tile : "; 
   // int szt; cin >> szt;
   cout << "-> give the verbose level (0 for silent) : "; 
   int vrblvl; cin >> vrblvl;

   int **nvr = new int*[dim];
   for(int i=0; i<dim; i++) nvr[i] = new int[dim];

   make_cyclic_variables(dim,dim,nvr);

   int ***idx = new int**[dim]; // we have dim columns and
   for(int i=0; i<dim; i++)     // dim monomials in each column
   {
      idx[i] = new int*[dim];
      for(int j=0; j<dim; j++) idx[i][j] = new int[nvr[i][j]];
   }
   make_cyclic_columns(dim,dim,nvr,idx);

   write_cyclic_columns(dim,dim,nvr,idx);

   int fail = 0;

   fail += dbl2_real_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);
   fail += dbl2_complex_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);

   if(fail == 0) cout << "-> Test passed." << endl;

   return 0;
}

int dbl2_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cffhi = new double**[dim];
   double ***cfflo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffhi[i] = new double*[dim]; // the coefficients of monomials
      cfflo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffhi[i][j] = new double[degp1];
         cfflo[i][j] = new double[degp1];
      }
   }
   double **solhi = new double*[dim];
   double **sollo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solhi[i] = new double[degp1];
      sollo[i] = new double[degp1];
   }
   make_real2_exponentials(dim,deg,solhi,sollo);
   make_real2_coefficients(dim,dim,cffhi,cfflo);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++)
         {
            cffhi[i][j][k] = 0.0;
            cfflo[i][j][k] = 0.0;
         }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << solhi[i][j] << "  "
                 << sollo[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cffhi[i][j][k] << "  "
                    << cfflo[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

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
   double ***outputhi_h;
   double ***outputlo_h;
   double ***outputhi_d;
   double ***outputlo_d;

   outputhi_h = new double**[dim];
   outputlo_h = new double**[dim];
   outputhi_d = new double**[dim];
   outputlo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhi_h[i] = new double*[dim+1];
      outputlo_h[i] = new double*[dim+1];
      outputhi_d[i] = new double*[dim+1];
      outputlo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhi_h[i][j] = new double[degp1];
         outputlo_h[i][j] = new double[degp1];
         outputhi_d[i][j] = new double[degp1];
         outputlo_d[i][j] = new double[degp1];
      }
   }
   double **funvalhi_h;  // function values on host
   double **funvallo_h;
   double **funvalhi_d;  // function values on device
   double **funvallo_d;

   funvalhi_h = new double*[dim];
   funvallo_h = new double*[dim];
   funvalhi_d = new double*[dim];
   funvallo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhi_h[i] = new double[degp1];
      funvallo_h[i] = new double[degp1];
      funvalhi_d[i] = new double[degp1];
      funvallo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhi_h;
   double ***jacvallo_h;
   double ***jacvalhi_d;
   double ***jacvallo_d;

   jacvalhi_h = new double**[degp1];
   jacvallo_h = new double**[degp1];
   jacvalhi_d = new double**[degp1];
   jacvallo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhi_h[i] = new double*[dim];
      jacvallo_h[i] = new double*[dim];
      jacvalhi_d[i] = new double*[dim];
      jacvallo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhi_h[i][j] = new double[dim];
         jacvallo_h[i][j] = new double[dim];
         jacvalhi_d[i][j] = new double[dim];
         jacvallo_d[i][j] = new double[dim];
      }
   }
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **acchi = new double*[dim+1]; // accumulate series in one column
   double **acclo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchi[i] = new double[degp1];
      acclo[i] = new double[degp1];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhi_h[i][j] = solhi[i][j];
         inputlo_h[i][j] = sollo[i][j];
         inputhi_d[i][j] = inputhi_h[i][j];
         inputlo_d[i][j] = inputlo_h[i][j];
      }

   CPU_dbl2_evaluate_columns
      (dim,deg,dim,nvr,idx,cffhi,cfflo,acchi,acclo,inputhi_h,inputlo_h,
       funvalhi_h,funvallo_h,jacvalhi_h,jacvallo_h,0);

   double totcnvlapsedms;

   GPU_dbl2_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,cffhi,cfflo,inputhi_d,inputlo_d,
       outputhi_d,outputlo_d,funvalhi_d,funvallo_d,jacvalhi_d,jacvallo_d,
       &totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = dbl2_error2sum(dim,degp1,
                funvalhi_h,funvallo_h,funvalhi_d,funvallo_d,"funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = dbl2_error3sum(degp1,dim,dim,
                jacvalhi_h,jacvallo_h,jacvalhi_d,jacvallo_d,"jacval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-12);
}

int dbl2_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cffrehi = new double**[dim];
   double ***cffrelo = new double**[dim];
   double ***cffimhi = new double**[dim];
   double ***cffimlo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffrehi[i] = new double*[dim]; // the coefficients of monomials
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
   make_complex2_coefficients(dim,dim,cffrehi,cffrelo,cffimhi,cffimlo);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++)
         {
            cffrehi[i][j][k] = 0.0; cffrelo[i][j][k] = 0.0;
            cffimhi[i][j][k] = 0.0; cffimlo[i][j][k] = 0.0;
         }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << solrehi[i][j] << "  "
                 << solrelo[i][j] << endl
                 << solimhi[i][j] << "  "
                 << solimlo[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cffrehi[i][j][k] << "  "
                    << cffrelo[i][j][k] << endl
                    << cffimhi[i][j][k] << "  "
                    << cffimlo[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

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
   double ***outputrehi_h;
   double ***outputrelo_h;
   double ***outputimhi_h;
   double ***outputimlo_h;
   double ***outputrehi_d;
   double ***outputrelo_d;
   double ***outputimhi_d;
   double ***outputimlo_d;

   outputrehi_h = new double**[dim];
   outputrelo_h = new double**[dim];
   outputimhi_h = new double**[dim];
   outputimlo_h = new double**[dim];
   outputrehi_d = new double**[dim];
   outputrelo_d = new double**[dim];
   outputimhi_d = new double**[dim];
   outputimlo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehi_h[i] = new double*[dim+1];
      outputrelo_h[i] = new double*[dim+1];
      outputimhi_h[i] = new double*[dim+1];
      outputimlo_h[i] = new double*[dim+1];
      outputrehi_d[i] = new double*[dim+1];
      outputrelo_d[i] = new double*[dim+1];
      outputimhi_d[i] = new double*[dim+1];
      outputimlo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehi_h[i][j] = new double[degp1];
         outputrelo_h[i][j] = new double[degp1];
         outputimhi_h[i][j] = new double[degp1];
         outputimlo_h[i][j] = new double[degp1];
         outputrehi_d[i][j] = new double[degp1];
         outputrelo_d[i][j] = new double[degp1];
         outputimhi_d[i][j] = new double[degp1];
         outputimlo_d[i][j] = new double[degp1];
      }
   }
   double **funvalrehi_h;  // function values on host
   double **funvalrelo_h;
   double **funvalimhi_h;
   double **funvalimlo_h;
   double **funvalrehi_d;  // function values on device
   double **funvalrelo_d;
   double **funvalimhi_d;
   double **funvalimlo_d;

   funvalrehi_h = new double*[dim];
   funvalrelo_h = new double*[dim];
   funvalimhi_h = new double*[dim];
   funvalimlo_h = new double*[dim];
   funvalrehi_d = new double*[dim];
   funvalrelo_d = new double*[dim];
   funvalimhi_d = new double*[dim];
   funvalimlo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalrehi_h[i] = new double[degp1];
      funvalrelo_h[i] = new double[degp1];
      funvalimhi_h[i] = new double[degp1];
      funvalimlo_h[i] = new double[degp1];
      funvalrehi_d[i] = new double[degp1];
      funvalrelo_d[i] = new double[degp1];
      funvalimhi_d[i] = new double[degp1];
      funvalimlo_d[i] = new double[degp1];

      for(int j=0; j<degp1; j++)
      {
         funvalrehi_h[i][j] = 0.0; funvalrelo_h[i][j] = 0.0;
         funvalimhi_h[i][j] = 0.0; funvalimlo_h[i][j] = 0.0;
         funvalrehi_d[i][j] = 0.0; funvalrelo_d[i][j] = 0.0;
         funvalimhi_d[i][j] = 0.0; funvalimlo_d[i][j] = 0.0;
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

   jacvalrehi_h = new double**[degp1];
   jacvalrelo_h = new double**[degp1];
   jacvalimhi_h = new double**[degp1];
   jacvalimlo_h = new double**[degp1];
   jacvalrehi_d = new double**[degp1];
   jacvalrelo_d = new double**[degp1];
   jacvalimhi_d = new double**[degp1];
   jacvalimlo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalrehi_h[i] = new double*[dim];
      jacvalrelo_h[i] = new double*[dim];
      jacvalimhi_h[i] = new double*[dim];
      jacvalimlo_h[i] = new double*[dim];
      jacvalrehi_d[i] = new double*[dim];
      jacvalrelo_d[i] = new double*[dim];
      jacvalimhi_d[i] = new double*[dim];
      jacvalimlo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalrehi_h[i][j] = new double[dim];
         jacvalrelo_h[i][j] = new double[dim];
         jacvalimhi_h[i][j] = new double[dim];
         jacvalimlo_h[i][j] = new double[dim];
         jacvalrehi_d[i][j] = new double[dim];
         jacvalrelo_d[i][j] = new double[dim];
         jacvalimhi_d[i][j] = new double[dim];
         jacvalimlo_d[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            jacvalrehi_h[i][j][k] = 0.0; jacvalrelo_h[i][j][k] = 0.0;
            jacvalimhi_h[i][j][k] = 0.0; jacvalimlo_h[i][j][k] = 0.0;
            jacvalrehi_d[i][j][k] = 0.0; jacvalrelo_d[i][j][k] = 0.0;
            jacvalimhi_d[i][j][k] = 0.0; jacvalimlo_d[i][j][k] = 0.0;
         }
      }
   }
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **accrehi = new double*[dim+1]; // accumulate series in one column
   double **accrelo = new double*[dim+1];
   double **accimhi = new double*[dim+1];
   double **accimlo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehi[i] = new double[degp1];
      accrelo[i] = new double[degp1];
      accimhi[i] = new double[degp1];
      accimlo[i] = new double[degp1];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehi_h[i][j] = solrehi[i][j];
         inputrelo_h[i][j] = solrelo[i][j];
         inputimhi_h[i][j] = solimhi[i][j];
         inputimlo_h[i][j] = solimlo[i][j];
         inputrehi_d[i][j] = inputrehi_h[i][j];
         inputrelo_d[i][j] = inputrelo_h[i][j];
         inputimhi_d[i][j] = inputimhi_h[i][j];
         inputimlo_d[i][j] = inputimlo_h[i][j];
      }

   CPU_cmplx2_evaluate_columns
      (dim,deg,dim,nvr,idx,cffrehi,cffrelo,cffimhi,cffimlo,
       accrehi,accrelo,accimhi,accimlo,
       inputrehi_h,inputrelo_h,inputimhi_h,inputimlo_h,
       funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
       jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,0);

   double totcnvlapsedms;

   GPU_cmplx2_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,cffrehi,cffrelo,cffimhi,cffimlo,
       inputrehi_d,inputrelo_d,inputimhi_d,inputimlo_d,
       outputrehi_d,outputrelo_d,outputimhi_d,outputimlo_d,
       funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
       jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
       &totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = cmplx2_error2sum
                (dim,degp1,
                 funvalrehi_h,funvalrelo_h,funvalimhi_h,funvalimlo_h,
                 funvalrehi_d,funvalrelo_d,funvalimhi_d,funvalimlo_d,
                 "funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = cmplx2_error3sum
                (degp1,dim,dim,
                 jacvalrehi_h,jacvalrelo_h,jacvalimhi_h,jacvalimlo_h,
                 jacvalrehi_d,jacvalrelo_d,jacvalimhi_d,jacvalimlo_d,
                 "jacval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-12);
}
