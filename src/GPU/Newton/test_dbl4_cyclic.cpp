/* Tests evaluation and differentiation of the column representation of 
 * the cyclic n-roots system in quad double precision.
 * In this column representation, the polynomial system is a sequence
 * of monomial systems, each column defines one monomial system.
 * Each column can be evaluated and differentiated independently,
 * just as each monomial in every column. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "random_numbers.h"
#include "dbl4_systems_host.h"
#include "dbl4_monomial_systems.h"
#include "dbl4_systems_kernels.h"
#include "dbl4_newton_testers.h"
#include "cyclic_columns.h"

using namespace std;

int dbl4_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns
 *   on real data in quad double precision, both on CPU and GPU.
 *
 * ON ENTRY :
 *   nbt      number of tiles;
 *   szt      size of each tile;
 *   deg      degree of the series;
 *   dim      dimension of the problem;
 *   nvr      number of variables for each monomial;
 *   idx      index representation of the monomials;
 *   vrblvl   is the verbose level. */

int dbl4_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl );
/*
 * DESCRIPTION :
 *   Generates a random series, truncated at degree deg,
 *   evaluates and differentiation the cyclic dim-roots columns,
 *   on complex data in quad double precision, both on CPU and GPU.
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

   // fail += dbl4_real_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);
   fail += dbl4_complex_evaltest(nbt,szt,deg,dim,nvr,idx,vrblvl);

   if(fail == 0)
      cout << "-> Test passed." << endl;
   else
      cout << "-> " << fail << " tests failed." << endl;

   return 0;
}

int dbl4_real_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> setting up the test system ..." << endl;

   double ***cffhihi = new double**[dim];
   double ***cfflohi = new double**[dim];
   double ***cffhilo = new double**[dim];
   double ***cfflolo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihi[i] = new double*[dim]; // the coefficients of monomials
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
   double **solhihi = new double*[dim];
   double **sollohi = new double*[dim];
   double **solhilo = new double*[dim];
   double **sollolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solhihi[i] = new double[degp1];
      sollohi[i] = new double[degp1];
      solhilo[i] = new double[degp1];
      sollolo[i] = new double[degp1];
   }
   make_real4_exponentials(dim,deg,solhihi,sollohi,solhilo,sollolo);
   make_real4_coefficients(dim,dim,cffhihi,cfflohi,cffhilo,cfflolo);
   for(int i=0; i<dim; i++)
      for(int j=0; j<dim; j++)
         for(int k=1; k<degp1; k++)
         {
            cffhihi[i][j][k] = 0.0; cfflohi[i][j][k] = 0.0;
            cffhilo[i][j][k] = 0.0; cfflolo[i][j][k] = 0.0;
         }

   if(vrblvl > 1)
   {
      cout << scientific << setprecision(16);
      cout << "The coefficients of the solution series :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<degp1; j++)
            cout << "sol[" << i << "][" << j << "] : "
                 << solhihi[i][j] << "  "
                 << sollohi[i][j] << endl
                 << solhilo[i][j] << "  "
                 << sollolo[i][j] << endl;
      }
      cout << "The coefficients of the system :" << endl;
      for(int i=0; i<dim; i++)
      {
         for(int j=0; j<dim; j++)
            for(int k=0; k<degp1; k++)
               cout << "cff[" << i << "][" << j << "][" << k << "] : "
                    << cffhihi[i][j][k] << "  "
                    << cfflohi[i][j][k] << endl
                    << cffhilo[i][j][k] << "  "
                    << cfflolo[i][j][k] << endl;
      }
   }

   cout << "-> allocating space for input and output ..." << endl;

   double **inputhihi_h = new double*[dim];
   double **inputlohi_h = new double*[dim];
   double **inputhilo_h = new double*[dim];
   double **inputlolo_h = new double*[dim];
   double **inputhihi_d = new double*[dim];
   double **inputlohi_d = new double*[dim];
   double **inputhilo_d = new double*[dim];
   double **inputlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
       inputhihi_h[i] = new double[degp1];
       inputlohi_h[i] = new double[degp1];
       inputhilo_h[i] = new double[degp1];
       inputlolo_h[i] = new double[degp1];
       inputhihi_d[i] = new double[degp1];
       inputlohi_d[i] = new double[degp1];
       inputhilo_d[i] = new double[degp1];
       inputlolo_d[i] = new double[degp1];
   }
   double ***outputhihi_h;
   double ***outputlohi_h;
   double ***outputhilo_h;
   double ***outputlolo_h;
   double ***outputhihi_d;
   double ***outputlohi_d;
   double ***outputhilo_d;
   double ***outputlolo_d;

   outputhihi_h = new double**[dim];
   outputlohi_h = new double**[dim];
   outputhilo_h = new double**[dim];
   outputlolo_h = new double**[dim];
   outputhihi_d = new double**[dim];
   outputlohi_d = new double**[dim];
   outputhilo_d = new double**[dim];
   outputlolo_d = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      outputhihi_h[i] = new double*[dim+1];
      outputlohi_h[i] = new double*[dim+1];
      outputhilo_h[i] = new double*[dim+1];
      outputlolo_h[i] = new double*[dim+1];
      outputhihi_d[i] = new double*[dim+1];
      outputlohi_d[i] = new double*[dim+1];
      outputhilo_d[i] = new double*[dim+1];
      outputlolo_d[i] = new double*[dim+1];

      for(int j=0; j<=dim; j++)
      {
         outputhihi_h[i][j] = new double[degp1];
         outputlohi_h[i][j] = new double[degp1];
         outputhilo_h[i][j] = new double[degp1];
         outputlolo_h[i][j] = new double[degp1];
         outputhihi_d[i][j] = new double[degp1];
         outputlohi_d[i][j] = new double[degp1];
         outputhilo_d[i][j] = new double[degp1];
         outputlolo_d[i][j] = new double[degp1];
      }
   }
   int **rowsA = new int*[dim]; // exponents in the rows
   for(int i=0; i<dim; i++) rowsA[i] = new int[dim]; // work space

   double **rhshihi; // result of evaluate columns
   double **rhslohi;
   double **rhshilo;
   double **rhslolo;

   rhshihi = new double*[dim];
   rhslohi = new double*[dim];
   rhshilo = new double*[dim];
   rhslolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhshihi[i] = new double[degp1];
      rhslohi[i] = new double[degp1];
      rhshilo[i] = new double[degp1];
      rhslolo[i] = new double[degp1];
   }
   double **funvalhihi_h;  // function values on host
   double **funvallohi_h;
   double **funvalhilo_h;
   double **funvallolo_h;
   double **funvalhihi_d;  // function values on device
   double **funvallohi_d;
   double **funvalhilo_d;
   double **funvallolo_d;

   funvalhihi_h = new double*[dim];
   funvallohi_h = new double*[dim];
   funvalhilo_h = new double*[dim];
   funvallolo_h = new double*[dim];
   funvalhihi_d = new double*[dim];
   funvallohi_d = new double*[dim];
   funvalhilo_d = new double*[dim];
   funvallolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      funvalhihi_h[i] = new double[degp1];
      funvallohi_h[i] = new double[degp1];
      funvalhilo_h[i] = new double[degp1];
      funvallolo_h[i] = new double[degp1];
      funvalhihi_d[i] = new double[degp1];
      funvallohi_d[i] = new double[degp1];
      funvalhilo_d[i] = new double[degp1];
      funvallolo_d[i] = new double[degp1];
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
   double ***jacvalhihi_h;
   double ***jacvallohi_h;
   double ***jacvalhilo_h;
   double ***jacvallolo_h;
   double ***jacvalhihi_d;
   double ***jacvallohi_d;
   double ***jacvalhilo_d;
   double ***jacvallolo_d;

   jacvalhihi_h = new double**[degp1];
   jacvallohi_h = new double**[degp1];
   jacvalhilo_h = new double**[degp1];
   jacvallolo_h = new double**[degp1];
   jacvalhihi_d = new double**[degp1];
   jacvallohi_d = new double**[degp1];
   jacvalhilo_d = new double**[degp1];
   jacvallolo_d = new double**[degp1];

   for(int i=0; i<degp1; i++) // jacval[i] is matrix of dimension dim
   {
      jacvalhihi_h[i] = new double*[dim];
      jacvallohi_h[i] = new double*[dim];
      jacvalhilo_h[i] = new double*[dim];
      jacvallolo_h[i] = new double*[dim];
      jacvalhihi_d[i] = new double*[dim];
      jacvallohi_d[i] = new double*[dim];
      jacvalhilo_d[i] = new double*[dim];
      jacvallolo_d[i] = new double*[dim];

      for(int j=0; j<dim; j++)
      {
         jacvalhihi_h[i][j] = new double[dim];
         jacvallohi_h[i][j] = new double[dim];
         jacvalhilo_h[i][j] = new double[dim];
         jacvallolo_h[i][j] = new double[dim];
         jacvalhihi_d[i][j] = new double[dim];
         jacvallohi_d[i][j] = new double[dim];
         jacvalhilo_d[i][j] = new double[dim];
         jacvallolo_d[i][j] = new double[dim];
      }
   }
   cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

   double **acchihi = new double*[dim+1]; // accumulate series in one column
   double **acclohi = new double*[dim+1];
   double **acchilo = new double*[dim+1];
   double **acclolo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      acchihi[i] = new double[degp1];
      acclohi[i] = new double[degp1];
      acchilo[i] = new double[degp1];
      acclolo[i] = new double[degp1];
   }
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
      {
         inputhihi_h[i][j] = solhihi[i][j];
         inputlohi_h[i][j] = sollohi[i][j];
         inputhilo_h[i][j] = solhilo[i][j];
         inputlolo_h[i][j] = sollolo[i][j];
         inputhihi_d[i][j] = inputhihi_h[i][j];
         inputlohi_d[i][j] = inputlohi_h[i][j];
         inputhilo_d[i][j] = inputhilo_h[i][j];
         inputlolo_d[i][j] = inputlolo_h[i][j];
      }

   evaluate_real4_columns
      (dim,deg,dim,nvr,idx,rowsA,
       cffhihi,cfflohi,cffhilo,cfflolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       rhshihi,rhslohi,rhshilo,rhslolo,0);

   CPU_dbl4_evaluate_columns
      (dim,deg,dim,nvr,idx,
       cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
       inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h,
       funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
       jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing column evaluations ... " << endl;
   double errsum0 = 0.0;

   errsum0 = dbl4_error2sum(dim,degp1,
                rhshihi,rhslohi,rhshilo,rhslolo,
                funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
                "rhsfun",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum0 << endl;

   double totcnvlapsedms;

   GPU_dbl4_evaluate_columns
      (dim,deg,dim,szt,nbt,nvr,idx,
       cffhihi,cfflohi,cffhilo,cfflolo,
       inputhihi_d,inputlohi_d,inputhilo_d,inputlolo_d,
       outputhihi_d,outputlohi_d,outputhilo_d,outputlolo_d,
       funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
       jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
       &totcnvlapsedms,0);

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU evaluations ... " << endl;
   double errsum1 = 0.0;

   errsum1 = dbl4_error2sum(dim,degp1,
                funvalhihi_h,funvallohi_h,funvalhilo_h,funvallolo_h,
                funvalhihi_d,funvallohi_d,funvalhilo_d,funvallolo_d,
                "funval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum1 << endl;

   cout << scientific << setprecision(16);
   cout << "-> comparing CPU with GPU Jacobians ... " << endl;
   double errsum2 = 0.0;

   errsum2 = dbl4_error3sum(degp1,dim,dim,
                jacvalhihi_h,jacvallohi_h,jacvalhilo_h,jacvallolo_h,
                jacvalhihi_d,jacvallohi_d,jacvalhilo_d,jacvallolo_d,
                "jacval",vrblvl);

   cout << scientific << setprecision(3);
   cout << "sum of errors : " << errsum2 << endl;

   double errsum = errsum0 + errsum1 + errsum2;
   cout << "total sum of errors : " << errsum << endl;

   return (errsum > 1.0e-48);
}

int dbl4_complex_evaltest
 ( int nbt, int szt, int deg, int dim, int **nvr, int ***idx, int vrblvl )
{
   const int degp1 = deg+1;

   cout << "-> allocating coefficients and solution series ..." << endl;

   double ***cffrehihi = new double**[dim];
   double ***cffrelohi = new double**[dim];
   double ***cffrehilo = new double**[dim];
   double ***cffrelolo = new double**[dim];
   double ***cffimhihi = new double**[dim];
   double ***cffimlohi = new double**[dim];
   double ***cffimhilo = new double**[dim];
   double ***cffimlolo = new double**[dim];

   for(int i=0; i<dim; i++)
   {
      cffrehihi[i] = new double*[dim]; // the coefficients of monomials
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
   double **solrehihi = new double*[dim];
   double **solrelohi = new double*[dim];
   double **solrehilo = new double*[dim];
   double **solrelolo = new double*[dim];
   double **solimhihi = new double*[dim];
   double **solimlohi = new double*[dim];
   double **solimhilo = new double*[dim];
   double **solimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      solrehihi[i] = new double[degp1];
      solrelohi[i] = new double[degp1];
      solrehilo[i] = new double[degp1];
      solrelolo[i] = new double[degp1];
      solimhihi[i] = new double[degp1];
      solimlohi[i] = new double[degp1];
      solimhilo[i] = new double[degp1];
      solimlolo[i] = new double[degp1];
   }
   cout << "-> allocating space for input and output ..." << endl;

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

   outputrehihi_h = new double**[dim];
   outputrelohi_h = new double**[dim];
   outputrehilo_h = new double**[dim];
   outputrelolo_h = new double**[dim];
   outputimhihi_h = new double**[dim];
   outputimlohi_h = new double**[dim];
   outputimhilo_h = new double**[dim];
   outputimlolo_h = new double**[dim];
   outputrehihi_d = new double**[dim];
   outputrelohi_d = new double**[dim];
   outputrehilo_d = new double**[dim];
   outputrelolo_d = new double**[dim];
   outputimhihi_d = new double**[dim];
   outputimlohi_d = new double**[dim];
   outputimhilo_d = new double**[dim];
   outputimlolo_d = new double**[dim];

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
   int **rowsA = new int*[dim]; // exponents in the rows
   for(int i=0; i<dim; i++) rowsA[i] = new int[dim]; // work space

   double **rhsrehihi; // result of evaluate columns
   double **rhsrelohi;
   double **rhsrehilo;
   double **rhsrelolo;
   double **rhsimhihi;
   double **rhsimlohi;
   double **rhsimhilo;
   double **rhsimlolo;

   rhsrehihi = new double*[dim];
   rhsrelohi = new double*[dim];
   rhsrehilo = new double*[dim];
   rhsrelolo = new double*[dim];
   rhsimhihi = new double*[dim];
   rhsimlohi = new double*[dim];
   rhsimhilo = new double*[dim];
   rhsimlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhsrehihi[i] = new double[degp1];
      rhsrelohi[i] = new double[degp1];
      rhsrehilo[i] = new double[degp1];
      rhsrelolo[i] = new double[degp1];
      rhsimhihi[i] = new double[degp1];
      rhsimlohi[i] = new double[degp1];
      rhsimhilo[i] = new double[degp1];
      rhsimlolo[i] = new double[degp1];
   }

   double **funvalrehihi_h;  // function values on host
   double **funvalrelohi_h;
   double **funvalrehilo_h;
   double **funvalrelolo_h;
   double **funvalimhihi_h;
   double **funvalimlohi_h;
   double **funvalimhilo_h;
   double **funvalimlolo_h;
   double **funvalrehihi_d;  // function values on device
   double **funvalrelohi_d;
   double **funvalrehilo_d;
   double **funvalrelolo_d;
   double **funvalimhihi_d;
   double **funvalimlohi_d;
   double **funvalimhilo_d;
   double **funvalimlolo_d;

   funvalrehihi_h = new double*[dim];
   funvalrelohi_h = new double*[dim];
   funvalrehilo_h = new double*[dim];
   funvalrelolo_h = new double*[dim];
   funvalimhihi_h = new double*[dim];
   funvalimlohi_h = new double*[dim];
   funvalimhilo_h = new double*[dim];
   funvalimlolo_h = new double*[dim];
   funvalrehihi_d = new double*[dim];
   funvalrelohi_d = new double*[dim];
   funvalrehilo_d = new double*[dim];
   funvalrelolo_d = new double*[dim];
   funvalimhihi_d = new double*[dim];
   funvalimlohi_d = new double*[dim];
   funvalimhilo_d = new double*[dim];
   funvalimlolo_d = new double*[dim];

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

      for(int j=0; j<degp1; j++)
      {
         funvalrehihi_h[i][j] = 0.0; funvalrelohi_h[i][j] = 0.0;
         funvalrehilo_h[i][j] = 0.0; funvalrelolo_h[i][j] = 0.0;
         funvalimhihi_h[i][j] = 0.0; funvalimlohi_h[i][j] = 0.0;
         funvalimhilo_h[i][j] = 0.0; funvalimlolo_h[i][j] = 0.0;
         funvalrehihi_d[i][j] = 0.0; funvalrelohi_d[i][j] = 0.0;
         funvalrehilo_d[i][j] = 0.0; funvalrelolo_d[i][j] = 0.0;
         funvalimhihi_d[i][j] = 0.0; funvalimlohi_d[i][j] = 0.0;
         funvalimhilo_d[i][j] = 0.0; funvalimlolo_d[i][j] = 0.0;
      }
   }
   // The derivatives in the output are a series truncated at degree deg.
   // The coefficients of the series are matrices of dimension dim.
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

   jacvalrehihi_h = new double**[degp1];
   jacvalrelohi_h = new double**[degp1];
   jacvalrehilo_h = new double**[degp1];
   jacvalrelolo_h = new double**[degp1];
   jacvalimhihi_h = new double**[degp1];
   jacvalimlohi_h = new double**[degp1];
   jacvalimhilo_h = new double**[degp1];
   jacvalimlolo_h = new double**[degp1];
   jacvalrehihi_d = new double**[degp1];
   jacvalrelohi_d = new double**[degp1];
   jacvalrehilo_d = new double**[degp1];
   jacvalrelolo_d = new double**[degp1];
   jacvalimhihi_d = new double**[degp1];
   jacvalimlohi_d = new double**[degp1];
   jacvalimhilo_d = new double**[degp1];
   jacvalimlolo_d = new double**[degp1];

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

         for(int k=0; k<dim; k++)
         {
            jacvalrehihi_h[i][j][k] = 0.0; jacvalrelohi_h[i][j][k] = 0.0;
            jacvalrehilo_h[i][j][k] = 0.0; jacvalrelolo_h[i][j][k] = 0.0;
            jacvalimhihi_h[i][j][k] = 0.0; jacvalimlohi_h[i][j][k] = 0.0;
            jacvalimhilo_h[i][j][k] = 0.0; jacvalimlolo_h[i][j][k] = 0.0;
            jacvalrehihi_d[i][j][k] = 0.0; jacvalrelohi_d[i][j][k] = 0.0;
            jacvalrehilo_d[i][j][k] = 0.0; jacvalrelolo_d[i][j][k] = 0.0;
            jacvalimhihi_d[i][j][k] = 0.0; jacvalimlohi_d[i][j][k] = 0.0;
            jacvalimhilo_d[i][j][k] = 0.0; jacvalimlolo_d[i][j][k] = 0.0;
         }
      }
   }
   double **accrehihi = new double*[dim+1]; // to accumulate series
   double **accrelohi = new double*[dim+1];
   double **accrehilo = new double*[dim+1];
   double **accrelolo = new double*[dim+1];
   double **accimhihi = new double*[dim+1];
   double **accimlohi = new double*[dim+1];
   double **accimhilo = new double*[dim+1];
   double **accimlolo = new double*[dim+1];

   for(int i=0; i<=dim; i++)
   {
      accrehihi[i] = new double[degp1]; accrelohi[i] = new double[degp1];
      accrehilo[i] = new double[degp1]; accrelolo[i] = new double[degp1];
      accimhihi[i] = new double[degp1]; accimlohi[i] = new double[degp1];
      accimhilo[i] = new double[degp1]; accimlolo[i] = new double[degp1];
   }
   double errsum = 0.0;

   for(int run=0; run<2; run++)
   {
      cout << "-> defining input and coefficients ..." << endl;

      make_complex4_exponentials
         (dim,deg,solrehihi,solrelohi,solrehilo,solrelolo,
                  solimhihi,solimlohi,solimhilo,solimlolo);

/*
   for(int i=0; i<dim; i++)  // set leading coefficients to random doubles
   { 
      solrehihi[i][0] = random_double();
      solimhihi[i][0] = random_double();
   }

 */     make_complex4_coefficients
         (dim,dim,cffrehihi,cffrelohi,cffrehilo,cffrelolo,
                  cffimhihi,cffimlohi,cffimhilo,cffimlolo);

      for(int i=0; i<dim; i++)
         for(int j=0; j<dim; j++)
         {
            // reset coefficient to one
            if(run == 0)
            {
               cffrehihi[i][j][0] = 1.0; cffrelohi[i][j][0] = 0.0;
               cffrehilo[i][j][0] = 0.0; cffrelolo[i][j][0] = 0.0;
               cffimhihi[i][j][0] = 0.0; cffimlohi[i][j][0] = 0.0;
               cffimhilo[i][j][0] = 0.0; cffimlolo[i][j][0] = 0.0;
            }
            for(int k=1; k<degp1; k++)
            {
               cffrehihi[i][j][k] = 0.0; cffrelohi[i][j][k] = 0.0;
               cffrehilo[i][j][k] = 0.0; cffrelolo[i][j][k] = 0.0;
               cffimhihi[i][j][k] = 0.0; cffimlohi[i][j][k] = 0.0;
               cffimhilo[i][j][k] = 0.0; cffimlolo[i][j][k] = 0.0;
            }
         }

      if(vrblvl > 1)
      {
         cout << scientific << setprecision(16);
         cout << "The coefficients of the solution series :" << endl;
         for(int i=0; i<dim; i++)
         {
            for(int j=0; j<degp1; j++)
               cout << "sol[" << i << "][" << j << "] : "
                    << solrehihi[i][j] << "  "
                    << solrelohi[i][j] << endl
                    << solrehilo[i][j] << "  "
                    << solrelolo[i][j] << endl
                    << solimhihi[i][j] << "  "
                    << solimlohi[i][j] << endl
                    << solimhilo[i][j] << "  "
                    << solimlolo[i][j] << endl;
         }
         cout << "The coefficients of the system :" << endl;
         for(int i=0; i<dim; i++)
         {
            for(int j=0; j<dim; j++)
               for(int k=0; k<degp1; k++)
                  cout << "cff[" << i << "][" << j << "][" << k << "] : "
                       << cffrehihi[i][j][k] << "  "
                       << cffrelohi[i][j][k] << endl
                       << cffrehilo[i][j][k] << "  "
                       << cffrelolo[i][j][k] << endl
                       << cffimhihi[i][j][k] << "  "
                       << cffimlohi[i][j][k] << endl
                       << cffimhilo[i][j][k] << "  "
                       << cffimlolo[i][j][k] << endl;
         }
      }
      cout << "-> calling CPU_dbl_evaluate_monomials ..." << endl;

      for(int i=0; i<dim; i++)
         for(int j=0; j<degp1; j++)
         {
            inputrehihi_h[i][j] = solrehihi[i][j];
            inputrelohi_h[i][j] = solrelohi[i][j];
            inputrehilo_h[i][j] = solrehilo[i][j];
            inputrelolo_h[i][j] = solrelolo[i][j];
            inputimhihi_h[i][j] = solimhihi[i][j];
            inputimlohi_h[i][j] = solimlohi[i][j];
            inputimhilo_h[i][j] = solimhilo[i][j];
            inputimlolo_h[i][j] = solimlolo[i][j];
            inputrehihi_d[i][j] = inputrehihi_h[i][j];
            inputrelohi_d[i][j] = inputrelohi_h[i][j];
            inputrehilo_d[i][j] = inputrehilo_h[i][j];
            inputrelolo_d[i][j] = inputrelolo_h[i][j];
            inputimhihi_d[i][j] = inputimhihi_h[i][j];
            inputimlohi_d[i][j] = inputimlohi_h[i][j];
            inputimhilo_d[i][j] = inputimhilo_h[i][j];
            inputimlolo_d[i][j] = inputimlolo_h[i][j];
         }

      const int nbrcol = 1;

      evaluate_complex4_columns
         (dim,deg,nbrcol,nvr,idx,rowsA,
          cffrehihi,cffrelohi,cffrehilo,cffrelolo,
          cffimhihi,cffimlohi,cffimhilo,cffimlolo,
          inputrehihi_h,inputrelohi_h,inputrehilo_h,inputrelolo_h,
          inputimhihi_h,inputimlohi_h,inputimhilo_h,inputimlolo_h,
          rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
          rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,0);

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
          jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,0);

      cout << "*** run " << run << " ***" << endl;

      cout << scientific << setprecision(16);
      cout << "-> comparing column evaluations ... " << endl;
      double errsum0 = 0.0;

      errsum0 = cmplx4_error2sum(dim,degp1,
                   rhsrehihi,rhsrelohi,rhsrehilo,rhsrelolo,
                   rhsimhihi,rhsimlohi,rhsimhilo,rhsimlolo,
                   funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
                   funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
                   "rhsfun",vrblvl);

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum0 << endl;

      double totcnvlapsedms;

      if(vrblvl > 0) cout << scientific << setprecision(16);

      GPU_cmplx4vectorized_evaluate_columns
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
          &totcnvlapsedms,vrblvl);

      cout << scientific << setprecision(16);
      cout << "-> comparing CPU with GPU evaluations ... " << endl;
      double errsum1 = 0.0;

      errsum1 = cmplx4_error2sum (dim,degp1,
         funvalrehihi_h,funvalrelohi_h,funvalrehilo_h,funvalrelolo_h,
         funvalimhihi_h,funvalimlohi_h,funvalimhilo_h,funvalimlolo_h,
         funvalrehihi_d,funvalrelohi_d,funvalrehilo_d,funvalrelolo_d,
         funvalimhihi_d,funvalimlohi_d,funvalimhilo_d,funvalimlolo_d,
         "funval",vrblvl);

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum1 << endl;

      cout << scientific << setprecision(16);
      cout << "-> comparing CPU with GPU Jacobians ... " << endl;
      double errsum2 = 0.0;

      errsum2 = cmplx4_error3sum(degp1,dim,dim,
         jacvalrehihi_h,jacvalrelohi_h,jacvalrehilo_h,jacvalrelolo_h,
         jacvalimhihi_h,jacvalimlohi_h,jacvalimhilo_h,jacvalimlolo_h,
         jacvalrehihi_d,jacvalrelohi_d,jacvalrehilo_d,jacvalrelolo_d,
         jacvalimhihi_d,jacvalimlohi_d,jacvalimhilo_d,jacvalimlolo_d,
         "jacval",vrblvl);

      cout << scientific << setprecision(3);
      cout << "sum of errors : " << errsum2 << endl;

      errsum = errsum0 + errsum1 + errsum2;
      cout << "total sum of errors : " << errsum << endl;
   }
   return (errsum > 1.0e-48);
}
