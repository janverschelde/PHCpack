/* Tests evaluation and differentiation of the column representation of 
 * the cyclic n-roots system in double precision.
 * In this column representation, the polynomial system is a sequence
 * of monomial systems, each column defines one monomial system.
 * Each column can be evaluated and differentiated independently,
 * just as each monomial in every column. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "dbl_systems_host.h"
#include "dbl_monomial_systems.h"
#include "dbl_systems_kernels.h"
#include "dbl_newton_testers.h"
#include "cyclic_columns.h"

using namespace std;

int dbl_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns
 *   on real data, both on CPU and GPU.
 *
 * ON ENTRY :
 *   nbt      number of tiles;
 *   szt      size of each tile;
 *   deg      degree of the series;
 *   dim      dimension of the problem;
 *   nvr      number of variables for each monomial;
 *   idx      index representation of the monomials;
 *   vrblvl   is the verbose level. */

int dbl_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns,
 *   on complex data, both on CPU and GPU.
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

   fail += dbl_real_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);
   fail += dbl_complex_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);

   if(fail == 0) cout << "-> Test passed." << endl;

   return 0;
}

int dbl_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cff = new double**[dim];
   for(int i=0; i<dim; i++)
   {
      cff[i] = new double*[dim]; // the coefficients of monomials
      for(int j=0; j<dim; j++) cff[i][j] = new double[degp1];
   }
   double **sol = new double*[dim];

   for(int i=0; i<dim; i++) sol[i] = new double[degp1];

   make_real_exponentials(dim,deg,sol);
   make_real_coefficients(dim,dim,cff);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++) cff[i][j][k] = 0.0;

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << sol[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cff[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

   double **input_h = new double*[dim];
   double **input_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
       input_h[i] = new double[degp1];
       input_d[i] = new double[degp1];
   }
   double ***output_h;
   double ***output_d;

   output_h = new double**[dim];
   output_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      output_h[i] = new double*[dim+1];
      output_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         output_h[i][j] = new double[degp1];
         output_d[i][j] = new double[degp1];
      }
   }
   double **funval_h;  // function values on host
   double **funval_d;  // function values on device

   funval_h = new double*[dim];
   funval_d = new double*[dim];
   for(int i=0; i<dim; i++)
   {
      funval_h[i] = new double[degp1];
      funval_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacval_h;
   double ***jacval_d;

   jacval_h = new double**[degp1];
   jacval_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacval_h[i] = new double*[dim];
      jacval_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacval_h[i][j] = new double[dim];
         jacval_d[i][j] = new double[dim];
      }
   }
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **acc = new double*[dim+1]; // accumulate series in one column
   for(int i=0; i<=dim; i++) acc[i] = new double[degp1];

   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         input_h[i][j] = sol[i][j];
         input_d[i][j] = input_h[i][j];
      }

   CPU_dbl_evaluate_columns
      (dim,deg,dim,nvr,idx,cff,acc,input_h,funval_h,jacval_h,0);

   double totcnvlapsedms;

   GPU_dbl_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,cff,input_d,output_d,
       funval_d,jacval_d,&totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = dbl_error2sum(dim,degp1,funval_h,funval_d,"funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = dbl_error3sum(degp1,dim,dim,jacval_h,jacval_d,"jacval",vrblvl);
   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-12);
}

int dbl_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cffre = new double**[dim];
   double ***cffim = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffre[i] = new double*[dim]; // the coefficients of monomials
      cffim[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffre[i][j] = new double[degp1];
         cffim[i][j] = new double[degp1];
      }
   }
   double *angles = new double[dim];
   double **solre = new double*[dim];
   double **solim = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solre[i] = new double[degp1];
      solim[i] = new double[degp1];
   }
   make_complex_exponentials(dim,deg,angles,solre,solim);
   make_complex_coefficients(dim,dim,cffre,cffim);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++)
         {
            cffre[i][j][k] = 0.0;
            cffim[i][j][k] = 0.0;
         }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << solre[i][j] << "  "
                 << solim[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cffre[i][j][k] << "  "
                    << cffim[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

   double **inputre_h = new double*[dim];
   double **inputim_h = new double*[dim];
   double **inputre_d = new double*[dim];
   double **inputim_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputre_h[i] = new double[degp1];
       inputim_h[i] = new double[degp1];
       inputre_d[i] = new double[degp1];
       inputim_d[i] = new double[degp1];
   }
   double ***outputre_h;
   double ***outputim_h;
   double ***outputre_d;
   double ***outputim_d;

   outputre_h = new double**[dim];
   outputim_h = new double**[dim];
   outputre_d = new double**[dim];
   outputim_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputre_h[i] = new double*[dim+1];
      outputim_h[i] = new double*[dim+1];
      outputre_d[i] = new double*[dim+1];
      outputim_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputre_h[i][j] = new double[degp1];
         outputim_h[i][j] = new double[degp1];
         outputre_d[i][j] = new double[degp1];
         outputim_d[i][j] = new double[degp1];
      }
   }
   double **funvalre_h;  // function values on host
   double **funvalim_h;
   double **funvalre_d;  // function values on device
   double **funvalim_d;

   funvalre_h = new double*[dim];
   funvalim_h = new double*[dim];
   funvalre_d = new double*[dim];
   funvalim_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalre_h[i] = new double[degp1];
      funvalim_h[i] = new double[degp1];
      funvalre_d[i] = new double[degp1];
      funvalim_d[i] = new double[degp1];

      for(int j=0; j<degp1; j++)
      {
         funvalre_h[i][j] = 0.0; funvalim_h[i][j] = 0.0;
         funvalre_d[i][j] = 0.0; funvalim_d[i][j] = 0.0;
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalre_h;
   double ***jacvalim_h;
   double ***jacvalre_d;
   double ***jacvalim_d;

   jacvalre_h = new double**[degp1];
   jacvalim_h = new double**[degp1];
   jacvalre_d = new double**[degp1];
   jacvalim_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalre_h[i] = new double*[dim];
      jacvalim_h[i] = new double*[dim];
      jacvalre_d[i] = new double*[dim];
      jacvalim_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalre_h[i][j] = new double[dim];
         jacvalim_h[i][j] = new double[dim];
         jacvalre_d[i][j] = new double[dim];
         jacvalim_d[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            jacvalre_h[i][j][k] = 0.0; jacvalim_h[i][j][k] = 0.0;
            jacvalre_d[i][j][k] = 0.0; jacvalim_d[i][j][k] = 0.0;
         }
      }
   }
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **accre = new double*[dim+1]; // accumulate series in one column
   double **accim = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accre[i] = new double[degp1];
      accim[i] = new double[degp1];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputre_h[i][j] = solre[i][j];
         inputim_h[i][j] = solim[i][j];
         inputre_d[i][j] = inputre_h[i][j];
         inputim_d[i][j] = inputim_h[i][j];
      }

   CPU_cmplx_evaluate_columns
      (dim,deg,dim,nvr,idx,cffre,cffim,accre,accim,inputre_h,inputim_h,
       funvalre_h,funvalim_h,jacvalre_h,jacvalim_h,0);

   double totcnvlapsedms;

   GPU_cmplx_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,cffre,cffim,inputre_d,inputim_d,
       outputre_d,outputim_d,funvalre_d,funvalim_d,
       jacvalre_d,jacvalim_d,&totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = cmplx_error2sum
                (dim,degp1,funvalre_h,funvalim_h,
                           funvalre_d,funvalim_d,"funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = cmplx_error3sum
                (degp1,dim,dim,jacvalre_h,jacvalim_h,
                               jacvalre_d,jacvalim_d,"jacval",vrblvl);
   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-12);
}
