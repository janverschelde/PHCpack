/* Tests evaluation and differentiation of the column representation of 
 * the cyclic n-roots system in octo double precision.
 * In this column representation, the polynomial system is a sequence
 * of monomial systems, each column defines one monomial system.
 * Each column can be evaluated and differentiated independently,
 * just as each monomial in every column. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "dbl8_systems_host.h"
#include "dbl8_monomial_systems.h"
#include "dbl8_systems_kernels.h"
#include "dbl8_newton_testers.h"
#include "cyclic_columns.h"

using namespace std;

int dbl8_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns
 *   on real data in octo double precision, both on CPU and GPU.
 *
 * ON ENTRY :
 *   nbt      number of tiles;
 *   szt      size of each tile;
 *   deg      degree of the series;
 *   dim      dimension of the problem;
 *   nvr      number of variables for each monomial;
 *   idx      index representation of the monomials;
 *   vrblvl   is the verbose level. */

int dbl8_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns,
 *   on complex data in octo double precision, both on CPU and GPU.
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

   fail += dbl8_real_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);
   fail += dbl8_complex_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);

   if(fail == 0) cout << "-> Test passed." << endl;

   return 0;
}

int dbl8_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cffhihihi = new double**[dim];
   double ***cfflohihi = new double**[dim];
   double ***cffhilohi = new double**[dim];
   double ***cfflolohi = new double**[dim];
   double ***cffhihilo = new double**[dim];
   double ***cfflohilo = new double**[dim];
   double ***cffhilolo = new double**[dim];
   double ***cfflololo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double*[dim]; // the coefficients of monomials
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
   double **solhihihi = new double*[dim];
   double **sollohihi = new double*[dim];
   double **solhilohi = new double*[dim];
   double **sollolohi = new double*[dim];
   double **solhihilo = new double*[dim];
   double **sollohilo = new double*[dim];
   double **solhilolo = new double*[dim];
   double **sollololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solhihihi[i] = new double[degp1];
      sollohihi[i] = new double[degp1];
      solhilohi[i] = new double[degp1];
      sollolohi[i] = new double[degp1];
      solhihilo[i] = new double[degp1];
      sollohilo[i] = new double[degp1];
      solhilolo[i] = new double[degp1];
      sollololo[i] = new double[degp1];
   }
   make_real8_exponentials
      (dim,deg,solhihihi,sollohihi,solhilohi,sollolohi,
               solhihilo,sollohilo,solhilolo,sollololo);

   make_real8_coefficients
      (dim,dim,cffhihihi,cfflohihi,cffhilohi,cfflolohi,
               cffhihilo,cfflohilo,cffhilolo,cfflololo);

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++)
         {
            cffhihihi[i][j][k] = 0.0; cfflohihi[i][j][k] = 0.0;
            cffhilohi[i][j][k] = 0.0; cfflolohi[i][j][k] = 0.0;
            cffhihilo[i][j][k] = 0.0; cfflohilo[i][j][k] = 0.0;
            cffhilolo[i][j][k] = 0.0; cfflololo[i][j][k] = 0.0;
         }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << solhihihi[i][j] << "  "
                 << sollohihi[i][j] << endl
                 << solhilohi[i][j] << "  "
                 << sollolohi[i][j] << endl
                 << solhihilo[i][j] << "  "
                 << sollohilo[i][j] << endl
                 << solhilolo[i][j] << "  "
                 << sollololo[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cffhihihi[i][j][k] << "  "
                    << cfflohihi[i][j][k] << endl
                    << cffhilohi[i][j][k] << "  "
                    << cfflolohi[i][j][k] << endl
                    << cffhihilo[i][j][k] << "  "
                    << cfflohilo[i][j][k] << endl
                    << cffhilolo[i][j][k] << "  "
                    << cfflololo[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

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

   outputhihihi_h = new double**[dim];
   outputlohihi_h = new double**[dim];
   outputhilohi_h = new double**[dim];
   outputlolohi_h = new double**[dim];
   outputhihilo_h = new double**[dim];
   outputlohilo_h = new double**[dim];
   outputhilolo_h = new double**[dim];
   outputlololo_h = new double**[dim];
   outputhihihi_d = new double**[dim];
   outputlohihi_d = new double**[dim];
   outputhilohi_d = new double**[dim];
   outputlolohi_d = new double**[dim];
   outputhihilo_d = new double**[dim];
   outputlohilo_d = new double**[dim];
   outputhilolo_d = new double**[dim];
   outputlololo_d = new double**[dim];

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
   double **funvalhihihi_h;  // function values on host
   double **funvallohihi_h;
   double **funvalhilohi_h;
   double **funvallolohi_h;
   double **funvalhihilo_h;
   double **funvallohilo_h;
   double **funvalhilolo_h;
   double **funvallololo_h;
   double **funvalhihihi_d;  // function values on device
   double **funvallohihi_d;
   double **funvalhilohi_d;
   double **funvallolohi_d;
   double **funvalhihilo_d;
   double **funvallohilo_d;
   double **funvalhilolo_d;
   double **funvallololo_d;

   funvalhihihi_h = new double*[dim];
   funvallohihi_h = new double*[dim];
   funvalhilohi_h = new double*[dim];
   funvallolohi_h = new double*[dim];
   funvalhihilo_h = new double*[dim];
   funvallohilo_h = new double*[dim];
   funvalhilolo_h = new double*[dim];
   funvallololo_h = new double*[dim];
   funvalhihihi_d = new double*[dim];
   funvallohihi_d = new double*[dim];
   funvalhilohi_d = new double*[dim];
   funvallolohi_d = new double*[dim];
   funvalhihilo_d = new double*[dim];
   funvallohilo_d = new double*[dim];
   funvalhilolo_d = new double*[dim];
   funvallololo_d = new double*[dim];

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

   jacvalhihihi_h = new double**[degp1];
   jacvallohihi_h = new double**[degp1];
   jacvalhilohi_h = new double**[degp1];
   jacvallolohi_h = new double**[degp1];
   jacvalhihilo_h = new double**[degp1];
   jacvallohilo_h = new double**[degp1];
   jacvalhilolo_h = new double**[degp1];
   jacvallololo_h = new double**[degp1];
   jacvalhihihi_d = new double**[degp1];
   jacvallohihi_d = new double**[degp1];
   jacvalhilohi_d = new double**[degp1];
   jacvallolohi_d = new double**[degp1];
   jacvalhihilo_d = new double**[degp1];
   jacvallohilo_d = new double**[degp1];
   jacvalhilolo_d = new double**[degp1];
   jacvallololo_d = new double**[degp1];

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
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **acchihihi = new double*[dim+1]; // to accumulate series
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
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihihi_h[i][j] = solhihihi[i][j];
         inputlohihi_h[i][j] = sollohihi[i][j];
         inputhilohi_h[i][j] = solhilohi[i][j];
         inputlolohi_h[i][j] = sollolohi[i][j];
         inputhihilo_h[i][j] = solhihilo[i][j];
         inputlohilo_h[i][j] = sollohilo[i][j];
         inputhilolo_h[i][j] = solhilolo[i][j];
         inputlololo_h[i][j] = sollololo[i][j];
         inputhihihi_d[i][j] = inputhihihi_h[i][j];
         inputlohihi_d[i][j] = inputlohihi_h[i][j];
         inputhilohi_d[i][j] = inputhilohi_h[i][j];
         inputlolohi_d[i][j] = inputlolohi_h[i][j];
         inputhihilo_d[i][j] = inputhihilo_h[i][j];
         inputlohilo_d[i][j] = inputlohilo_h[i][j];
         inputhilolo_d[i][j] = inputhilolo_h[i][j];
         inputlololo_d[i][j] = inputlololo_h[i][j];
      }

   CPU_dbl8_evaluate_columns
      (dim,deg,dim,nvr,idx,
       cffhihihi,cfflohihi,cffhilohi,cfflolohi,
       cffhihilo,cfflohilo,cffhilolo,cfflololo,
       acchihihi,acclohihi,acchilohi,acclolohi,
       acchihilo,acclohilo,acchilolo,acclololo,
       inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
       inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h,
       funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
       funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
       jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
       jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,0);

   double totcnvlapsedms;

   GPU_dbl8_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,
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
       &totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = dbl8_error2sum(dim,degp1,
                funvalhihihi_h,funvallohihi_h,funvalhilohi_h,funvallolohi_h,
                funvalhihilo_h,funvallohilo_h,funvalhilolo_h,funvallololo_h,
                funvalhihihi_d,funvallohihi_d,funvalhilohi_d,funvallolohi_d,
                funvalhihilo_d,funvallohilo_d,funvalhilolo_d,funvallololo_d,
                "funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = dbl8_error3sum(degp1,dim,dim,
                jacvalhihihi_h,jacvallohihi_h,jacvalhilohi_h,jacvallolohi_h,
                jacvalhihilo_h,jacvallohilo_h,jacvalhilolo_h,jacvallololo_h,
                jacvalhihihi_d,jacvallohihi_d,jacvalhilohi_d,jacvallolohi_d,
                jacvalhihilo_d,jacvallohilo_d,jacvalhilolo_d,jacvallololo_d,
                "jacval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-12);
}

int dbl8_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cffrehihihi = new double**[dim];
   double ***cffrelohihi = new double**[dim];
   double ***cffrehilohi = new double**[dim];
   double ***cffrelolohi = new double**[dim];
   double ***cffrehihilo = new double**[dim];
   double ***cffrelohilo = new double**[dim];
   double ***cffrehilolo = new double**[dim];
   double ***cffrelololo = new double**[dim];
   double ***cffimhihihi = new double**[dim];
   double ***cffimlohihi = new double**[dim];
   double ***cffimhilohi = new double**[dim];
   double ***cffimlolohi = new double**[dim];
   double ***cffimhihilo = new double**[dim];
   double ***cffimlohilo = new double**[dim];
   double ***cffimhilolo = new double**[dim];
   double ***cffimlololo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffrehihihi[i] = new double*[dim]; // the coefficients of monomials
      cffrelohihi[i] = new double*[dim];
      cffrehilohi[i] = new double*[dim];
      cffrelolohi[i] = new double*[dim];
      cffrehihilo[i] = new double*[dim];
      cffrelohilo[i] = new double*[dim];
      cffrehilolo[i] = new double*[dim];
      cffrelololo[i] = new double*[dim];
      cffimhihihi[i] = new double*[dim];
      cffimlohihi[i] = new double*[dim];
      cffimhilohi[i] = new double*[dim];
      cffimlolohi[i] = new double*[dim];
      cffimhihilo[i] = new double*[dim];
      cffimlohilo[i] = new double*[dim];
      cffimhilolo[i] = new double*[dim];
      cffimlololo[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         cffrehihihi[i][j] = new double[degp1];
         cffrelohihi[i][j] = new double[degp1];
         cffrehilohi[i][j] = new double[degp1];
         cffrelolohi[i][j] = new double[degp1];
         cffrehihilo[i][j] = new double[degp1];
         cffrelohilo[i][j] = new double[degp1];
         cffrehilolo[i][j] = new double[degp1];
         cffrelololo[i][j] = new double[degp1];
         cffimhihihi[i][j] = new double[degp1];
         cffimlohihi[i][j] = new double[degp1];
         cffimhilohi[i][j] = new double[degp1];
         cffimlolohi[i][j] = new double[degp1];
         cffimhihilo[i][j] = new double[degp1];
         cffimlohilo[i][j] = new double[degp1];
         cffimhilolo[i][j] = new double[degp1];
         cffimlololo[i][j] = new double[degp1];
      }
   }
   double **solrehihihi = new double*[dim];
   double **solrelohihi = new double*[dim];
   double **solrehilohi = new double*[dim];
   double **solrelolohi = new double*[dim];
   double **solrehihilo = new double*[dim];
   double **solrelohilo = new double*[dim];
   double **solrehilolo = new double*[dim];
   double **solrelololo = new double*[dim];
   double **solimhihihi = new double*[dim];
   double **solimlohihi = new double*[dim];
   double **solimhilohi = new double*[dim];
   double **solimlolohi = new double*[dim];
   double **solimhihilo = new double*[dim];
   double **solimlohilo = new double*[dim];
   double **solimhilolo = new double*[dim];
   double **solimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solrehihihi[i] = new double[degp1];
      solrelohihi[i] = new double[degp1];
      solrehilohi[i] = new double[degp1];
      solrelolohi[i] = new double[degp1];
      solrehihilo[i] = new double[degp1];
      solrelohilo[i] = new double[degp1];
      solrehilolo[i] = new double[degp1];
      solrelololo[i] = new double[degp1];
      solimhihihi[i] = new double[degp1];
      solimlohihi[i] = new double[degp1];
      solimhilohi[i] = new double[degp1];
      solimlolohi[i] = new double[degp1];
      solimhihilo[i] = new double[degp1];
      solimlohilo[i] = new double[degp1];
      solimhilolo[i] = new double[degp1];
      solimlololo[i] = new double[degp1];
   }
   make_complex8_exponentials
      (dim,deg,solrehihihi,solrelohihi,solrehilohi,solrelolohi,
               solrehihilo,solrelohilo,solrehilolo,solrelololo,
               solimhihihi,solimlohihi,solimhilohi,solimlolohi,
               solimhihilo,solimlohilo,solimhilolo,solimlololo);

   make_complex8_coefficients
      (dim,dim,cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
               cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
               cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
               cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo);

   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++)
         {
            cffrehihihi[i][j][k] = 0.0; cffrelohihi[i][j][k] = 0.0;
            cffrehilohi[i][j][k] = 0.0; cffrelolohi[i][j][k] = 0.0;
            cffrehihilo[i][j][k] = 0.0; cffrelohilo[i][j][k] = 0.0;
            cffrehilolo[i][j][k] = 0.0; cffrelololo[i][j][k] = 0.0;
            cffimhihihi[i][j][k] = 0.0; cffimlohihi[i][j][k] = 0.0;
            cffimhilohi[i][j][k] = 0.0; cffimlolohi[i][j][k] = 0.0;
            cffimhihilo[i][j][k] = 0.0; cffimlohilo[i][j][k] = 0.0;
            cffimhilolo[i][j][k] = 0.0; cffimlololo[i][j][k] = 0.0;
         }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << solrehihihi[i][j] << "  "
                 << solrelohihi[i][j] << endl
                 << solrehilohi[i][j] << "  "
                 << solrelolohi[i][j] << endl
                 << solrehihilo[i][j] << "  "
                 << solrelohilo[i][j] << endl
                 << solrehilolo[i][j] << "  "
                 << solrelololo[i][j] << endl
                 << solimhihihi[i][j] << "  "
                 << solimlohihi[i][j] << endl
                 << solimhilohi[i][j] << "  "
                 << solimlolohi[i][j] << endl
                 << solimhihilo[i][j] << "  "
                 << solimlohilo[i][j] << endl
                 << solimhilolo[i][j] << "  "
                 << solimlololo[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cffrehihihi[i][j][k] << "  "
                    << cffrelohihi[i][j][k] << endl
                    << cffrehilohi[i][j][k] << "  "
                    << cffrelolohi[i][j][k] << endl
                    << cffrehihilo[i][j][k] << "  "
                    << cffrelohilo[i][j][k] << endl
                    << cffrehilolo[i][j][k] << "  "
                    << cffrelololo[i][j][k] << endl
                    << cffimhihihi[i][j][k] << "  "
                    << cffimlohihi[i][j][k] << endl
                    << cffimhilohi[i][j][k] << "  "
                    << cffimlolohi[i][j][k] << endl
                    << cffimhihilo[i][j][k] << "  "
                    << cffimlohilo[i][j][k] << endl
                    << cffimhilolo[i][j][k] << "  "
                    << cffimlololo[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

   double **inputrehihihi_h = new double*[dim];
   double **inputrelohihi_h = new double*[dim];
   double **inputrehilohi_h = new double*[dim];
   double **inputrelolohi_h = new double*[dim];
   double **inputrehihilo_h = new double*[dim];
   double **inputrelohilo_h = new double*[dim];
   double **inputrehilolo_h = new double*[dim];
   double **inputrelololo_h = new double*[dim];
   double **inputimhihihi_h = new double*[dim];
   double **inputimlohihi_h = new double*[dim];
   double **inputimhilohi_h = new double*[dim];
   double **inputimlolohi_h = new double*[dim];
   double **inputimhihilo_h = new double*[dim];
   double **inputimlohilo_h = new double*[dim];
   double **inputimhilolo_h = new double*[dim];
   double **inputimlololo_h = new double*[dim];
   double **inputrehihihi_d = new double*[dim];
   double **inputrelohihi_d = new double*[dim];
   double **inputrehilohi_d = new double*[dim];
   double **inputrelolohi_d = new double*[dim];
   double **inputrehihilo_d = new double*[dim];
   double **inputrelohilo_d = new double*[dim];
   double **inputrehilolo_d = new double*[dim];
   double **inputrelololo_d = new double*[dim];
   double **inputimhihihi_d = new double*[dim];
   double **inputimlohihi_d = new double*[dim];
   double **inputimhilohi_d = new double*[dim];
   double **inputimlolohi_d = new double*[dim];
   double **inputimhihilo_d = new double*[dim];
   double **inputimlohilo_d = new double*[dim];
   double **inputimhilolo_d = new double*[dim];
   double **inputimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputrehihihi_h[i] = new double[degp1];
       inputrelohihi_h[i] = new double[degp1];
       inputrehilohi_h[i] = new double[degp1];
       inputrelolohi_h[i] = new double[degp1];
       inputrehihilo_h[i] = new double[degp1];
       inputrelohilo_h[i] = new double[degp1];
       inputrehilolo_h[i] = new double[degp1];
       inputrelololo_h[i] = new double[degp1];
       inputimhihihi_h[i] = new double[degp1];
       inputimlohihi_h[i] = new double[degp1];
       inputimhilohi_h[i] = new double[degp1];
       inputimlolohi_h[i] = new double[degp1];
       inputimhihilo_h[i] = new double[degp1];
       inputimlohilo_h[i] = new double[degp1];
       inputimhilolo_h[i] = new double[degp1];
       inputimlololo_h[i] = new double[degp1];
       inputrehihihi_d[i] = new double[degp1];
       inputrelohihi_d[i] = new double[degp1];
       inputrehilohi_d[i] = new double[degp1];
       inputrelolohi_d[i] = new double[degp1];
       inputrehihilo_d[i] = new double[degp1];
       inputrelohilo_d[i] = new double[degp1];
       inputrehilolo_d[i] = new double[degp1];
       inputrelololo_d[i] = new double[degp1];
       inputimhihihi_d[i] = new double[degp1];
       inputimlohihi_d[i] = new double[degp1];
       inputimhilohi_d[i] = new double[degp1];
       inputimlolohi_d[i] = new double[degp1];
       inputimhihilo_d[i] = new double[degp1];
       inputimlohilo_d[i] = new double[degp1];
       inputimhilolo_d[i] = new double[degp1];
       inputimlololo_d[i] = new double[degp1];
   }
   double ***outputrehihihi_h;
   double ***outputrelohihi_h;
   double ***outputrehilohi_h;
   double ***outputrelolohi_h;
   double ***outputrehihilo_h;
   double ***outputrelohilo_h;
   double ***outputrehilolo_h;
   double ***outputrelololo_h;
   double ***outputimhihihi_h;
   double ***outputimlohihi_h;
   double ***outputimhilohi_h;
   double ***outputimlolohi_h;
   double ***outputimhihilo_h;
   double ***outputimlohilo_h;
   double ***outputimhilolo_h;
   double ***outputimlololo_h;
   double ***outputrehihihi_d;
   double ***outputrelohihi_d;
   double ***outputrehilohi_d;
   double ***outputrelolohi_d;
   double ***outputrehihilo_d;
   double ***outputrelohilo_d;
   double ***outputrehilolo_d;
   double ***outputrelololo_d;
   double ***outputimhihihi_d;
   double ***outputimlohihi_d;
   double ***outputimhilohi_d;
   double ***outputimlolohi_d;
   double ***outputimhihilo_d;
   double ***outputimlohilo_d;
   double ***outputimhilolo_d;
   double ***outputimlololo_d;

   outputrehihihi_h = new double**[dim];
   outputrelohihi_h = new double**[dim];
   outputrehilohi_h = new double**[dim];
   outputrelolohi_h = new double**[dim];
   outputrehihilo_h = new double**[dim];
   outputrelohilo_h = new double**[dim];
   outputrehilolo_h = new double**[dim];
   outputrelololo_h = new double**[dim];
   outputimhihihi_h = new double**[dim];
   outputimlohihi_h = new double**[dim];
   outputimhilohi_h = new double**[dim];
   outputimlolohi_h = new double**[dim];
   outputimhihilo_h = new double**[dim];
   outputimlohilo_h = new double**[dim];
   outputimhilolo_h = new double**[dim];
   outputimlololo_h = new double**[dim];
   outputrehihihi_d = new double**[dim];
   outputrelohihi_d = new double**[dim];
   outputrehilohi_d = new double**[dim];
   outputrelolohi_d = new double**[dim];
   outputrehihilo_d = new double**[dim];
   outputrelohilo_d = new double**[dim];
   outputrehilolo_d = new double**[dim];
   outputrelololo_d = new double**[dim];
   outputimhihihi_d = new double**[dim];
   outputimlohihi_d = new double**[dim];
   outputimhilohi_d = new double**[dim];
   outputimlolohi_d = new double**[dim];
   outputimhihilo_d = new double**[dim];
   outputimlohilo_d = new double**[dim];
   outputimhilolo_d = new double**[dim];
   outputimlololo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputrehihihi_h[i] = new double*[dim+1];
      outputrelohihi_h[i] = new double*[dim+1];
      outputrehilohi_h[i] = new double*[dim+1];
      outputrelolohi_h[i] = new double*[dim+1];
      outputrehihilo_h[i] = new double*[dim+1];
      outputrelohilo_h[i] = new double*[dim+1];
      outputrehilolo_h[i] = new double*[dim+1];
      outputrelololo_h[i] = new double*[dim+1];
      outputimhihihi_h[i] = new double*[dim+1];
      outputimlohihi_h[i] = new double*[dim+1];
      outputimhilohi_h[i] = new double*[dim+1];
      outputimlolohi_h[i] = new double*[dim+1];
      outputimhihilo_h[i] = new double*[dim+1];
      outputimlohilo_h[i] = new double*[dim+1];
      outputimhilolo_h[i] = new double*[dim+1];
      outputimlololo_h[i] = new double*[dim+1];
      outputrehihihi_d[i] = new double*[dim+1];
      outputrelohihi_d[i] = new double*[dim+1];
      outputrehilohi_d[i] = new double*[dim+1];
      outputrelolohi_d[i] = new double*[dim+1];
      outputrehihilo_d[i] = new double*[dim+1];
      outputrelohilo_d[i] = new double*[dim+1];
      outputrehilolo_d[i] = new double*[dim+1];
      outputrelololo_d[i] = new double*[dim+1];
      outputimhihihi_d[i] = new double*[dim+1];
      outputimlohihi_d[i] = new double*[dim+1];
      outputimhilohi_d[i] = new double*[dim+1];
      outputimlolohi_d[i] = new double*[dim+1];
      outputimhihilo_d[i] = new double*[dim+1];
      outputimlohilo_d[i] = new double*[dim+1];
      outputimhilolo_d[i] = new double*[dim+1];
      outputimlololo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputrehihihi_h[i][j] = new double[degp1];
         outputrelohihi_h[i][j] = new double[degp1];
         outputrehilohi_h[i][j] = new double[degp1];
         outputrelolohi_h[i][j] = new double[degp1];
         outputrehihilo_h[i][j] = new double[degp1];
         outputrelohilo_h[i][j] = new double[degp1];
         outputrehilolo_h[i][j] = new double[degp1];
         outputrelololo_h[i][j] = new double[degp1];
         outputimhihihi_h[i][j] = new double[degp1];
         outputimlohihi_h[i][j] = new double[degp1];
         outputimhilohi_h[i][j] = new double[degp1];
         outputimlolohi_h[i][j] = new double[degp1];
         outputimhihilo_h[i][j] = new double[degp1];
         outputimlohilo_h[i][j] = new double[degp1];
         outputimhilolo_h[i][j] = new double[degp1];
         outputimlololo_h[i][j] = new double[degp1];
         outputrehihihi_d[i][j] = new double[degp1];
         outputrelohihi_d[i][j] = new double[degp1];
         outputrehilohi_d[i][j] = new double[degp1];
         outputrelolohi_d[i][j] = new double[degp1];
         outputrehihilo_d[i][j] = new double[degp1];
         outputrelohilo_d[i][j] = new double[degp1];
         outputrehilolo_d[i][j] = new double[degp1];
         outputrelololo_d[i][j] = new double[degp1];
         outputimhihihi_d[i][j] = new double[degp1];
         outputimlohihi_d[i][j] = new double[degp1];
         outputimhilohi_d[i][j] = new double[degp1];
         outputimlolohi_d[i][j] = new double[degp1];
         outputimhihilo_d[i][j] = new double[degp1];
         outputimlohilo_d[i][j] = new double[degp1];
         outputimhilolo_d[i][j] = new double[degp1];
         outputimlololo_d[i][j] = new double[degp1];
      }
   }
   double **funvalrehihihi_h;  // function values on host
   double **funvalrelohihi_h;
   double **funvalrehilohi_h;
   double **funvalrelolohi_h;
   double **funvalrehihilo_h;
   double **funvalrelohilo_h;
   double **funvalrehilolo_h;
   double **funvalrelololo_h;
   double **funvalimhihihi_h;
   double **funvalimlohihi_h;
   double **funvalimhilohi_h;
   double **funvalimlolohi_h;
   double **funvalimhihilo_h;
   double **funvalimlohilo_h;
   double **funvalimhilolo_h;
   double **funvalimlololo_h;
   double **funvalrehihihi_d;  // function values on device
   double **funvalrelohihi_d;
   double **funvalrehilohi_d;
   double **funvalrelolohi_d;
   double **funvalrehihilo_d;
   double **funvalrelohilo_d;
   double **funvalrehilolo_d;
   double **funvalrelololo_d;
   double **funvalimhihihi_d;
   double **funvalimlohihi_d;
   double **funvalimhilohi_d;
   double **funvalimlolohi_d;
   double **funvalimhihilo_d;
   double **funvalimlohilo_d;
   double **funvalimhilolo_d;
   double **funvalimlololo_d;

   funvalrehihihi_h = new double*[dim];
   funvalrelohihi_h = new double*[dim];
   funvalrehilohi_h = new double*[dim];
   funvalrelolohi_h = new double*[dim];
   funvalrehihilo_h = new double*[dim];
   funvalrelohilo_h = new double*[dim];
   funvalrehilolo_h = new double*[dim];
   funvalrelololo_h = new double*[dim];
   funvalimhihihi_h = new double*[dim];
   funvalimlohihi_h = new double*[dim];
   funvalimhilohi_h = new double*[dim];
   funvalimlolohi_h = new double*[dim];
   funvalimhihilo_h = new double*[dim];
   funvalimlohilo_h = new double*[dim];
   funvalimhilolo_h = new double*[dim];
   funvalimlololo_h = new double*[dim];
   funvalrehihihi_d = new double*[dim];
   funvalrelohihi_d = new double*[dim];
   funvalrehilohi_d = new double*[dim];
   funvalrelolohi_d = new double*[dim];
   funvalrehihilo_d = new double*[dim];
   funvalrelohilo_d = new double*[dim];
   funvalrehilolo_d = new double*[dim];
   funvalrelololo_d = new double*[dim];
   funvalimhihihi_d = new double*[dim];
   funvalimlohihi_d = new double*[dim];
   funvalimhilohi_d = new double*[dim];
   funvalimlolohi_d = new double*[dim];
   funvalimhihilo_d = new double*[dim];
   funvalimlohilo_d = new double*[dim];
   funvalimhilolo_d = new double*[dim];
   funvalimlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalrehihihi_h[i] = new double[degp1];
      funvalrelohihi_h[i] = new double[degp1];
      funvalrehilohi_h[i] = new double[degp1];
      funvalrelolohi_h[i] = new double[degp1];
      funvalrehihilo_h[i] = new double[degp1];
      funvalrelohilo_h[i] = new double[degp1];
      funvalrehilolo_h[i] = new double[degp1];
      funvalrelololo_h[i] = new double[degp1];
      funvalimhihihi_h[i] = new double[degp1];
      funvalimlohihi_h[i] = new double[degp1];
      funvalimhilohi_h[i] = new double[degp1];
      funvalimlolohi_h[i] = new double[degp1];
      funvalimhihilo_h[i] = new double[degp1];
      funvalimlohilo_h[i] = new double[degp1];
      funvalimhilolo_h[i] = new double[degp1];
      funvalimlololo_h[i] = new double[degp1];
      funvalrehihihi_d[i] = new double[degp1];
      funvalrelohihi_d[i] = new double[degp1];
      funvalrehilohi_d[i] = new double[degp1];
      funvalrelolohi_d[i] = new double[degp1];
      funvalrehihilo_d[i] = new double[degp1];
      funvalrelohilo_d[i] = new double[degp1];
      funvalrehilolo_d[i] = new double[degp1];
      funvalrelololo_d[i] = new double[degp1];
      funvalimhihihi_d[i] = new double[degp1];
      funvalimlohihi_d[i] = new double[degp1];
      funvalimhilohi_d[i] = new double[degp1];
      funvalimlolohi_d[i] = new double[degp1];
      funvalimhihilo_d[i] = new double[degp1];
      funvalimlohilo_d[i] = new double[degp1];
      funvalimhilolo_d[i] = new double[degp1];
      funvalimlololo_d[i] = new double[degp1];

      for(int j=0; j<degp1; j++)
      {
         funvalrehihihi_h[i][j] = 0.0; funvalrelohihi_h[i][j] = 0.0;
         funvalrehilohi_h[i][j] = 0.0; funvalrelolohi_h[i][j] = 0.0;
         funvalrehihilo_h[i][j] = 0.0; funvalrelohilo_h[i][j] = 0.0;
         funvalrehilolo_h[i][j] = 0.0; funvalrelololo_h[i][j] = 0.0;
         funvalimhihihi_h[i][j] = 0.0; funvalimlohihi_h[i][j] = 0.0;
         funvalimhilohi_h[i][j] = 0.0; funvalimlolohi_h[i][j] = 0.0;
         funvalimhihilo_h[i][j] = 0.0; funvalimlohilo_h[i][j] = 0.0;
         funvalimhilolo_h[i][j] = 0.0; funvalimlololo_h[i][j] = 0.0;
         funvalrehihihi_d[i][j] = 0.0; funvalrelohihi_d[i][j] = 0.0;
         funvalrehilohi_d[i][j] = 0.0; funvalrelolohi_d[i][j] = 0.0;
         funvalrehihilo_d[i][j] = 0.0; funvalrelohilo_d[i][j] = 0.0;
         funvalrehilolo_d[i][j] = 0.0; funvalrelololo_d[i][j] = 0.0;
         funvalimhihihi_d[i][j] = 0.0; funvalimlohihi_d[i][j] = 0.0;
         funvalimhilohi_d[i][j] = 0.0; funvalimlolohi_d[i][j] = 0.0;
         funvalimhihilo_d[i][j] = 0.0; funvalimlohilo_d[i][j] = 0.0;
         funvalimhilolo_d[i][j] = 0.0; funvalimlololo_d[i][j] = 0.0;
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalrehihihi_h;
   double ***jacvalrelohihi_h;
   double ***jacvalrehilohi_h;
   double ***jacvalrelolohi_h;
   double ***jacvalrehihilo_h;
   double ***jacvalrelohilo_h;
   double ***jacvalrehilolo_h;
   double ***jacvalrelololo_h;
   double ***jacvalimhihihi_h;
   double ***jacvalimlohihi_h;
   double ***jacvalimhilohi_h;
   double ***jacvalimlolohi_h;
   double ***jacvalimhihilo_h;
   double ***jacvalimlohilo_h;
   double ***jacvalimhilolo_h;
   double ***jacvalimlololo_h;
   double ***jacvalrehihihi_d;
   double ***jacvalrelohihi_d;
   double ***jacvalrehilohi_d;
   double ***jacvalrelolohi_d;
   double ***jacvalrehihilo_d;
   double ***jacvalrelohilo_d;
   double ***jacvalrehilolo_d;
   double ***jacvalrelololo_d;
   double ***jacvalimhihihi_d;
   double ***jacvalimlohihi_d;
   double ***jacvalimhilohi_d;
   double ***jacvalimlolohi_d;
   double ***jacvalimhihilo_d;
   double ***jacvalimlohilo_d;
   double ***jacvalimhilolo_d;
   double ***jacvalimlololo_d;

   jacvalrehihihi_h = new double**[degp1];
   jacvalrelohihi_h = new double**[degp1];
   jacvalrehilohi_h = new double**[degp1];
   jacvalrelolohi_h = new double**[degp1];
   jacvalrehihilo_h = new double**[degp1];
   jacvalrelohilo_h = new double**[degp1];
   jacvalrehilolo_h = new double**[degp1];
   jacvalrelololo_h = new double**[degp1];
   jacvalimhihihi_h = new double**[degp1];
   jacvalimlohihi_h = new double**[degp1];
   jacvalimhilohi_h = new double**[degp1];
   jacvalimlolohi_h = new double**[degp1];
   jacvalimhihilo_h = new double**[degp1];
   jacvalimlohilo_h = new double**[degp1];
   jacvalimhilolo_h = new double**[degp1];
   jacvalimlololo_h = new double**[degp1];
   jacvalrehihihi_d = new double**[degp1];
   jacvalrelohihi_d = new double**[degp1];
   jacvalrehilohi_d = new double**[degp1];
   jacvalrelolohi_d = new double**[degp1];
   jacvalrehihilo_d = new double**[degp1];
   jacvalrelohilo_d = new double**[degp1];
   jacvalrehilolo_d = new double**[degp1];
   jacvalrelololo_d = new double**[degp1];
   jacvalimhihihi_d = new double**[degp1];
   jacvalimlohihi_d = new double**[degp1];
   jacvalimhilohi_d = new double**[degp1];
   jacvalimlolohi_d = new double**[degp1];
   jacvalimhihilo_d = new double**[degp1];
   jacvalimlohilo_d = new double**[degp1];
   jacvalimhilolo_d = new double**[degp1];
   jacvalimlololo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalrehihihi_h[i] = new double*[dim];
      jacvalrelohihi_h[i] = new double*[dim];
      jacvalrehilohi_h[i] = new double*[dim];
      jacvalrelolohi_h[i] = new double*[dim];
      jacvalrehihilo_h[i] = new double*[dim];
      jacvalrelohilo_h[i] = new double*[dim];
      jacvalrehilolo_h[i] = new double*[dim];
      jacvalrelololo_h[i] = new double*[dim];
      jacvalimhihihi_h[i] = new double*[dim];
      jacvalimlohihi_h[i] = new double*[dim];
      jacvalimhilohi_h[i] = new double*[dim];
      jacvalimlolohi_h[i] = new double*[dim];
      jacvalimhihilo_h[i] = new double*[dim];
      jacvalimlohilo_h[i] = new double*[dim];
      jacvalimhilolo_h[i] = new double*[dim];
      jacvalimlololo_h[i] = new double*[dim];
      jacvalrehihihi_d[i] = new double*[dim];
      jacvalrelohihi_d[i] = new double*[dim];
      jacvalrehilohi_d[i] = new double*[dim];
      jacvalrelolohi_d[i] = new double*[dim];
      jacvalrehihilo_d[i] = new double*[dim];
      jacvalrelohilo_d[i] = new double*[dim];
      jacvalrehilolo_d[i] = new double*[dim];
      jacvalrelololo_d[i] = new double*[dim];
      jacvalimhihihi_d[i] = new double*[dim];
      jacvalimlohihi_d[i] = new double*[dim];
      jacvalimhilohi_d[i] = new double*[dim];
      jacvalimlolohi_d[i] = new double*[dim];
      jacvalimhihilo_d[i] = new double*[dim];
      jacvalimlohilo_d[i] = new double*[dim];
      jacvalimhilolo_d[i] = new double*[dim];
      jacvalimlololo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalrehihihi_h[i][j] = new double[dim];
         jacvalrelohihi_h[i][j] = new double[dim];
         jacvalrehilohi_h[i][j] = new double[dim];
         jacvalrelolohi_h[i][j] = new double[dim];
         jacvalrehihilo_h[i][j] = new double[dim];
         jacvalrelohilo_h[i][j] = new double[dim];
         jacvalrehilolo_h[i][j] = new double[dim];
         jacvalrelololo_h[i][j] = new double[dim];
         jacvalimhihihi_h[i][j] = new double[dim];
         jacvalimlohihi_h[i][j] = new double[dim];
         jacvalimhilohi_h[i][j] = new double[dim];
         jacvalimlolohi_h[i][j] = new double[dim];
         jacvalimhihilo_h[i][j] = new double[dim];
         jacvalimlohilo_h[i][j] = new double[dim];
         jacvalimhilolo_h[i][j] = new double[dim];
         jacvalimlololo_h[i][j] = new double[dim];
         jacvalrehihihi_d[i][j] = new double[dim];
         jacvalrelohihi_d[i][j] = new double[dim];
         jacvalrehilohi_d[i][j] = new double[dim];
         jacvalrelolohi_d[i][j] = new double[dim];
         jacvalrehihilo_d[i][j] = new double[dim];
         jacvalrelohilo_d[i][j] = new double[dim];
         jacvalrehilolo_d[i][j] = new double[dim];
         jacvalrelololo_d[i][j] = new double[dim];
         jacvalimhihihi_d[i][j] = new double[dim];
         jacvalimlohihi_d[i][j] = new double[dim];
         jacvalimhilohi_d[i][j] = new double[dim];
         jacvalimlolohi_d[i][j] = new double[dim];
         jacvalimhihilo_d[i][j] = new double[dim];
         jacvalimlohilo_d[i][j] = new double[dim];
         jacvalimhilolo_d[i][j] = new double[dim];
         jacvalimlololo_d[i][j] = new double[dim];

         for(int k=0; k<dim; k++)
         {
            jacvalrehihihi_h[i][j][k] = 0.0; jacvalrelohihi_h[i][j][k] = 0.0;
            jacvalrehilohi_h[i][j][k] = 0.0; jacvalrelolohi_h[i][j][k] = 0.0;
            jacvalrehihilo_h[i][j][k] = 0.0; jacvalrelohilo_h[i][j][k] = 0.0;
            jacvalrehilolo_h[i][j][k] = 0.0; jacvalrelololo_h[i][j][k] = 0.0;
            jacvalimhihihi_h[i][j][k] = 0.0; jacvalimlohihi_h[i][j][k] = 0.0;
            jacvalimhilohi_h[i][j][k] = 0.0; jacvalimlolohi_h[i][j][k] = 0.0;
            jacvalimhihilo_h[i][j][k] = 0.0; jacvalimlohilo_h[i][j][k] = 0.0;
            jacvalimhilolo_h[i][j][k] = 0.0; jacvalimlololo_h[i][j][k] = 0.0;
            jacvalrehihihi_d[i][j][k] = 0.0; jacvalrelohihi_d[i][j][k] = 0.0;
            jacvalrehilohi_d[i][j][k] = 0.0; jacvalrelolohi_d[i][j][k] = 0.0;
            jacvalrehihilo_d[i][j][k] = 0.0; jacvalrelohilo_d[i][j][k] = 0.0;
            jacvalrehilolo_d[i][j][k] = 0.0; jacvalrelololo_d[i][j][k] = 0.0;
            jacvalimhihihi_d[i][j][k] = 0.0; jacvalimlohihi_d[i][j][k] = 0.0;
            jacvalimhilohi_d[i][j][k] = 0.0; jacvalimlolohi_d[i][j][k] = 0.0;
            jacvalimhihilo_d[i][j][k] = 0.0; jacvalimlohilo_d[i][j][k] = 0.0;
            jacvalimhilolo_d[i][j][k] = 0.0; jacvalimlololo_d[i][j][k] = 0.0;
         }
      }
   }
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **accrehihihi = new double*[dim+1]; // to accumulate series
   double **accrelohihi = new double*[dim+1];
   double **accrehilohi = new double*[dim+1];
   double **accrelolohi = new double*[dim+1];
   double **accrehihilo = new double*[dim+1];
   double **accrelohilo = new double*[dim+1];
   double **accrehilolo = new double*[dim+1];
   double **accrelololo = new double*[dim+1];
   double **accimhihihi = new double*[dim+1];
   double **accimlohihi = new double*[dim+1];
   double **accimhilohi = new double*[dim+1];
   double **accimlolohi = new double*[dim+1];
   double **accimhihilo = new double*[dim+1];
   double **accimlohilo = new double*[dim+1];
   double **accimhilolo = new double*[dim+1];
   double **accimlololo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehihihi[i] = new double[degp1]; accrelohihi[i] = new double[degp1];
      accrehilohi[i] = new double[degp1]; accrelolohi[i] = new double[degp1];
      accrehihilo[i] = new double[degp1]; accrelohilo[i] = new double[degp1];
      accrehilolo[i] = new double[degp1]; accrelololo[i] = new double[degp1];
      accimhihihi[i] = new double[degp1]; accimlohihi[i] = new double[degp1];
      accimhilohi[i] = new double[degp1]; accimlolohi[i] = new double[degp1];
      accimhihilo[i] = new double[degp1]; accimlohilo[i] = new double[degp1];
      accimhilolo[i] = new double[degp1]; accimlololo[i] = new double[degp1];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputrehihihi_h[i][j] = solrehihihi[i][j];
         inputrelohihi_h[i][j] = solrelohihi[i][j];
         inputrehilohi_h[i][j] = solrehilohi[i][j];
         inputrelolohi_h[i][j] = solrelolohi[i][j];
         inputrehihilo_h[i][j] = solrehihilo[i][j];
         inputrelohilo_h[i][j] = solrelohilo[i][j];
         inputrehilolo_h[i][j] = solrehilolo[i][j];
         inputrelololo_h[i][j] = solrelololo[i][j];
         inputimhihihi_h[i][j] = solimhihihi[i][j];
         inputimlohihi_h[i][j] = solimlohihi[i][j];
         inputimhilohi_h[i][j] = solimhilohi[i][j];
         inputimlolohi_h[i][j] = solimlolohi[i][j];
         inputimhihilo_h[i][j] = solimhihilo[i][j];
         inputimlohilo_h[i][j] = solimlohilo[i][j];
         inputimhilolo_h[i][j] = solimhilolo[i][j];
         inputimlololo_h[i][j] = solimlololo[i][j];
         inputrehihihi_d[i][j] = inputrehihihi_h[i][j];
         inputrelohihi_d[i][j] = inputrelohihi_h[i][j];
         inputrehilohi_d[i][j] = inputrehilohi_h[i][j];
         inputrelolohi_d[i][j] = inputrelolohi_h[i][j];
         inputrehihilo_d[i][j] = inputrehihilo_h[i][j];
         inputrelohilo_d[i][j] = inputrelohilo_h[i][j];
         inputrehilolo_d[i][j] = inputrehilolo_h[i][j];
         inputrelololo_d[i][j] = inputrelololo_h[i][j];
         inputimhihihi_d[i][j] = inputimhihihi_h[i][j];
         inputimlohihi_d[i][j] = inputimlohihi_h[i][j];
         inputimhilohi_d[i][j] = inputimhilohi_h[i][j];
         inputimlolohi_d[i][j] = inputimlolohi_h[i][j];
         inputimhihilo_d[i][j] = inputimhihilo_h[i][j];
         inputimlohilo_d[i][j] = inputimlohilo_h[i][j];
         inputimhilolo_d[i][j] = inputimhilolo_h[i][j];
         inputimlololo_d[i][j] = inputimlololo_h[i][j];
      }

   CPU_cmplx8_evaluate_columns
      (dim,deg,dim,nvr,idx,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       accrehihihi,accrelohihi,accrehilohi,accrelolohi,
       accrehihilo,accrelohilo,accrehilolo,accrelololo,
       accimhihihi,accimlohihi,accimhilohi,accimlolohi,
       accimhihilo,accimlohilo,accimhilolo,accimlololo,
       inputrehihihi_h,inputrelohihi_h,inputrehilohi_h,inputrelolohi_h,
       inputrehihilo_h,inputrelohilo_h,inputrehilolo_h,inputrelololo_h,
       inputimhihihi_h,inputimlohihi_h,inputimhilohi_h,inputimlolohi_h,
       inputimhihilo_h,inputimlohilo_h,inputimhilolo_h,inputimlololo_h,
       funvalrehihihi_h,funvalrelohihi_h,funvalrehilohi_h,funvalrelolohi_h,
       funvalrehihilo_h,funvalrelohilo_h,funvalrehilolo_h,funvalrelololo_h,
       funvalimhihihi_h,funvalimlohihi_h,funvalimhilohi_h,funvalimlolohi_h,
       funvalimhihilo_h,funvalimlohilo_h,funvalimhilolo_h,funvalimlololo_h,
       jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
       jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
       jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
       jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,0);

   double totcnvlapsedms;

   GPU_cmplx8_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,
       cffrehihihi,cffrelohihi,cffrehilohi,cffrelolohi,
       cffrehihilo,cffrelohilo,cffrehilolo,cffrelololo,
       cffimhihihi,cffimlohihi,cffimhilohi,cffimlolohi,
       cffimhihilo,cffimlohilo,cffimhilolo,cffimlololo,
       inputrehihihi_d,inputrelohihi_d,inputrehilohi_d,inputrelolohi_d,
       inputrehihilo_d,inputrelohilo_d,inputrehilolo_d,inputrelololo_d,
       inputimhihihi_d,inputimlohihi_d,inputimhilohi_d,inputimlolohi_d,
       inputimhihilo_d,inputimlohilo_d,inputimhilolo_d,inputimlololo_d,
       outputrehihihi_d,outputrelohihi_d,outputrehilohi_d,outputrelolohi_d,
       outputrehihilo_d,outputrelohilo_d,outputrehilolo_d,outputrelololo_d,
       outputimhihihi_d,outputimlohihi_d,outputimhilohi_d,outputimlolohi_d,
       outputimhihilo_d,outputimlohilo_d,outputimhilolo_d,outputimlololo_d,
       funvalrehihihi_d,funvalrelohihi_d,funvalrehilohi_d,funvalrelolohi_d,
       funvalrehihilo_d,funvalrelohilo_d,funvalrehilolo_d,funvalrelololo_d,
       funvalimhihihi_d,funvalimlohihi_d,funvalimhilohi_d,funvalimlolohi_d,
       funvalimhihilo_d,funvalimlohilo_d,funvalimhilolo_d,funvalimlololo_d,
       jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
       jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
       jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
       jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
       &totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = cmplx8_error2sum(dim,degp1,
      funvalrehihihi_h,funvalrelohihi_h,funvalrehilohi_h,funvalrelolohi_h,
      funvalrehihilo_h,funvalrelohilo_h,funvalrehilolo_h,funvalrelololo_h,
      funvalimhihihi_h,funvalimlohihi_h,funvalimhilohi_h,funvalimlolohi_h,
      funvalimhihilo_h,funvalimlohilo_h,funvalimhilolo_h,funvalimlololo_h,
      funvalrehihihi_d,funvalrelohihi_d,funvalrehilohi_d,funvalrelolohi_d,
      funvalrehihilo_d,funvalrelohilo_d,funvalrehilolo_d,funvalrelololo_d,
      funvalimhihihi_d,funvalimlohihi_d,funvalimhilohi_d,funvalimlolohi_d,
      funvalimhihilo_d,funvalimlohilo_d,funvalimhilolo_d,funvalimlololo_d,
      "funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = cmplx8_error3sum(degp1,dim,dim,
      jacvalrehihihi_h,jacvalrelohihi_h,jacvalrehilohi_h,jacvalrelolohi_h,
      jacvalrehihilo_h,jacvalrelohilo_h,jacvalrehilolo_h,jacvalrelololo_h,
      jacvalimhihihi_h,jacvalimlohihi_h,jacvalimhilohi_h,jacvalimlolohi_h,
      jacvalimhihilo_h,jacvalimlohilo_h,jacvalimhilolo_h,jacvalimlololo_h,
      jacvalrehihihi_d,jacvalrelohihi_d,jacvalrehilohi_d,jacvalrelolohi_d,
      jacvalrehihilo_d,jacvalrelohilo_d,jacvalrehilolo_d,jacvalrelololo_d,
      jacvalimhihihi_d,jacvalimlohihi_d,jacvalimhilohi_d,jacvalimlolohi_d,
      jacvalimhihilo_d,jacvalimlohilo_d,jacvalimhilolo_d,jacvalimlololo_d,
      "jacval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-12);
}
