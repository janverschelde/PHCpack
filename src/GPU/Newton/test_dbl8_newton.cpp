/* Tests operations on Newton for series in octo double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl8_newton_testers.h"

using namespace std;

int main ( void )
{
   cout << "testing Newton ..." << endl;

   int dim,deg,size,posvals,vrblvl,nbritr,nbsteps,szt,nbt,mode;
   prompt_dimensions(&dim,&deg,&size,&posvals,&vrblvl,&nbritr,&nbsteps);

   cout << "-> enter 0 (GPU only), 1 (CPU only), or 2 (GPU+CPU) : ";
   cin >> mode;

   if(mode != 1)
   {
      cout << "-> give the number of tiles : "; cin >> nbt;
      cout << "-> give the size of each tile : "; cin >> szt;
      int p = szt*nbt;

      while(p != dim)
      {
          cout << "Dimension = " << dim << " != " << szt << " * " << nbt
               << ", retry." << endl;
          cout << "-> give the size of each tile : "; cin >> szt;
          cout << "-> give the number of tiles : "; cin >> nbt;
          p = szt*nbt;
      }
   }

   cout << "generating an integer matrix of dimension " << dim
        << " ..." << endl;
/*
 * 1. defining the unimodular matrix of monomials
 */
   int **rowsA = new int*[dim];  // exponents in the rows
   int *nvr = new int[dim];      // number of variables in each monomial
   int **idx = new int*[dim];    // indexes of variables in each monomial
   int **exp = new int*[dim];    // exponents of the variables
   int *nbrfac = new int[dim];   // number of exponents > 1 in each monomial
   int **expfac = new int*[dim]; // exponents of the common factors

   make_monomial_system
      (dim,size,posvals,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);
/*
 * 2. allocating input and output space for evaluation and differentiation
 */

   const int degp1 = deg+1;
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
   // allocate memory for coefficients and the output
   double *acchihihi = new double[degp1]; // accumulated power series
   double *acclohihi = new double[degp1];
   double *acchilohi = new double[degp1];
   double *acclolohi = new double[degp1];
   double *acchihilo = new double[degp1];
   double *acclohilo = new double[degp1];
   double *acchilolo = new double[degp1];
   double *acclololo = new double[degp1];
   double **cffhihihi = new double*[dim]; // the coefficients of monomials
   double **cfflohihi = new double*[dim];
   double **cffhilohi = new double*[dim];
   double **cfflolohi = new double*[dim];
   double **cffhihilo = new double*[dim];
   double **cfflohilo = new double*[dim];
   double **cffhilolo = new double*[dim];
   double **cfflololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihihi[i] = new double[degp1];
      cfflohihi[i] = new double[degp1];
      cffhilohi[i] = new double[degp1];
      cfflolohi[i] = new double[degp1];
      cffhihilo[i] = new double[degp1];
      cfflohilo[i] = new double[degp1];
      cffhilolo[i] = new double[degp1];
      cfflololo[i] = new double[degp1];
   }
   double ***outputhihihi_h = new double**[dim];
   double ***outputlohihi_h = new double**[dim];
   double ***outputhilohi_h = new double**[dim];
   double ***outputlolohi_h = new double**[dim];
   double ***outputhihilo_h = new double**[dim];
   double ***outputlohilo_h = new double**[dim];
   double ***outputhilolo_h = new double**[dim];
   double ***outputlololo_h = new double**[dim];
   double ***outputhihihi_d = new double**[dim];
   double ***outputlohihi_d = new double**[dim];
   double ***outputhilohi_d = new double**[dim];
   double ***outputlolohi_d = new double**[dim];
   double ***outputhihilo_d = new double**[dim];
   double ***outputlohilo_d = new double**[dim];
   double ***outputhilolo_d = new double**[dim];
   double ***outputlololo_d = new double**[dim];

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
   // The function values are power series truncated at degree deg.
   double **funvalhihihi_h = new double*[dim];
   double **funvallohihi_h = new double*[dim];
   double **funvalhilohi_h = new double*[dim];
   double **funvallolohi_h = new double*[dim];
   double **funvalhihilo_h = new double*[dim];
   double **funvallohilo_h = new double*[dim];
   double **funvalhilolo_h = new double*[dim];
   double **funvallololo_h = new double*[dim];
   double **funvalhihihi_d = new double*[dim];
   double **funvallohihi_d = new double*[dim];
   double **funvalhilohi_d = new double*[dim];
   double **funvallolohi_d = new double*[dim];
   double **funvalhihilo_d = new double*[dim];
   double **funvallohilo_d = new double*[dim];
   double **funvalhilolo_d = new double*[dim];
   double **funvallololo_d = new double*[dim];

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
   double ***jacvalhihihi_h = new double**[degp1];
   double ***jacvallohihi_h = new double**[degp1];
   double ***jacvalhilohi_h = new double**[degp1];
   double ***jacvallolohi_h = new double**[degp1];
   double ***jacvalhihilo_h = new double**[degp1];
   double ***jacvallohilo_h = new double**[degp1];
   double ***jacvalhilolo_h = new double**[degp1];
   double ***jacvallololo_h = new double**[degp1];
   double ***jacvalhihihi_d = new double**[degp1];
   double ***jacvallohihi_d = new double**[degp1];
   double ***jacvalhilohi_d = new double**[degp1];
   double ***jacvallolohi_d = new double**[degp1];
   double ***jacvalhihilo_d = new double**[degp1];
   double ***jacvallohilo_d = new double**[degp1];
   double ***jacvalhilolo_d = new double**[degp1];
   double ***jacvallololo_d = new double**[degp1];

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
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihihi_h = new double*[degp1];
   double **sollohihi_h = new double*[degp1];
   double **solhilohi_h = new double*[degp1];
   double **sollolohi_h = new double*[degp1];
   double **solhihilo_h = new double*[degp1];
   double **sollohilo_h = new double*[degp1];
   double **solhilolo_h = new double*[degp1];
   double **sollololo_h = new double*[degp1];
   double **solhihihi_d = new double*[degp1];
   double **sollohihi_d = new double*[degp1];
   double **solhilohi_d = new double*[degp1];
   double **sollolohi_d = new double*[degp1];
   double **solhihilo_d = new double*[degp1];
   double **sollohilo_d = new double*[degp1];
   double **solhilolo_d = new double*[degp1];
   double **sollololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihihi_h[i] = new double[dim];
      sollohihi_h[i] = new double[dim];
      solhilohi_h[i] = new double[dim];
      sollolohi_h[i] = new double[dim];
      solhihilo_h[i] = new double[dim];
      sollohilo_h[i] = new double[dim];
      solhilolo_h[i] = new double[dim];
      sollololo_h[i] = new double[dim];
      solhihihi_d[i] = new double[dim];
      sollohihi_d[i] = new double[dim];
      solhilohi_d[i] = new double[dim];
      sollolohi_d[i] = new double[dim];
      solhihilo_d[i] = new double[dim];
      sollohilo_d[i] = new double[dim];
      solhilolo_d[i] = new double[dim];
      sollololo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihihi_h = new double*[degp1];
   double **rhslohihi_h = new double*[degp1];
   double **rhshilohi_h = new double*[degp1];
   double **rhslolohi_h = new double*[degp1];
   double **rhshihilo_h = new double*[degp1];
   double **rhslohilo_h = new double*[degp1];
   double **rhshilolo_h = new double*[degp1];
   double **rhslololo_h = new double*[degp1];
   double **rhshihihi_d = new double*[degp1];
   double **rhslohihi_d = new double*[degp1];
   double **rhshilohi_d = new double*[degp1];
   double **rhslolohi_d = new double*[degp1];
   double **rhshihilo_d = new double*[degp1];
   double **rhslohilo_d = new double*[degp1];
   double **rhshilolo_d = new double*[degp1];
   double **rhslololo_d = new double*[degp1];

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
      rhshihihi_d[i] = new double[dim];
      rhslohihi_d[i] = new double[dim];
      rhshilohi_d[i] = new double[dim];
      rhslolohi_d[i] = new double[dim];
      rhshihilo_d[i] = new double[dim];
      rhslohilo_d[i] = new double[dim];
      rhshilolo_d[i] = new double[dim];
      rhslololo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihihi_h = new double*[degp1];
   double **urhslohihi_h = new double*[degp1];
   double **urhshilohi_h = new double*[degp1];
   double **urhslolohi_h = new double*[degp1];
   double **urhshihilo_h = new double*[degp1];
   double **urhslohilo_h = new double*[degp1];
   double **urhshilolo_h = new double*[degp1];
   double **urhslololo_h = new double*[degp1];
   double **urhshihihi_d = new double*[degp1];
   double **urhslohihi_d = new double*[degp1];
   double **urhshilohi_d = new double*[degp1];
   double **urhslolohi_d = new double*[degp1];
   double **urhshihilo_d = new double*[degp1];
   double **urhslohilo_d = new double*[degp1];
   double **urhshilolo_d = new double*[degp1];
   double **urhslololo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhshihihi_h[i] = new double[dim];
      urhslohihi_h[i] = new double[dim];
      urhshilohi_h[i] = new double[dim];
      urhslolohi_h[i] = new double[dim];
      urhshihilo_h[i] = new double[dim];
      urhslohilo_h[i] = new double[dim];
      urhshilolo_h[i] = new double[dim];
      urhslololo_h[i] = new double[dim];
      urhshihihi_d[i] = new double[dim];
      urhslohihi_d[i] = new double[dim];
      urhshilohi_d[i] = new double[dim];
      urhslolohi_d[i] = new double[dim];
      urhshihilo_d[i] = new double[dim];
      urhslohilo_d[i] = new double[dim];
      urhshilolo_d[i] = new double[dim];
      urhslololo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihihi = new double*[dim];
   double **workmatlohihi = new double*[dim];
   double **workmathilohi = new double*[dim];
   double **workmatlolohi = new double*[dim];
   double **workmathihilo = new double*[dim];
   double **workmatlohilo = new double*[dim];
   double **workmathilolo = new double*[dim];
   double **workmatlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihihi[i] = new double[dim];
      workmatlohihi[i] = new double[dim];
      workmathilohi[i] = new double[dim];
      workmatlolohi[i] = new double[dim];
      workmathihilo[i] = new double[dim];
      workmatlohilo[i] = new double[dim];
      workmathilolo[i] = new double[dim];
      workmatlololo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
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

   double **Qhihihi_h = new double*[dim];
   double **Qlohihi_h = new double*[dim];
   double **Qhilohi_h = new double*[dim];
   double **Qlolohi_h = new double*[dim];
   double **Qhihilo_h = new double*[dim];
   double **Qlohilo_h = new double*[dim];
   double **Qhilolo_h = new double*[dim];
   double **Qlololo_h = new double*[dim];
   double **Qhihihi_d = new double*[dim];
   double **Qlohihi_d = new double*[dim];
   double **Qhilohi_d = new double*[dim];
   double **Qlolohi_d = new double*[dim];
   double **Qhihilo_d = new double*[dim];
   double **Qlohilo_d = new double*[dim];
   double **Qhilolo_d = new double*[dim];
   double **Qlololo_d = new double*[dim];

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
      Qhihihi_d[i] = new double[dim];
      Qlohihi_d[i] = new double[dim];
      Qhilohi_d[i] = new double[dim];
      Qlolohi_d[i] = new double[dim];
      Qhihilo_d[i] = new double[dim];
      Qlohilo_d[i] = new double[dim];
      Qhilolo_d[i] = new double[dim];
      Qlololo_d[i] = new double[dim];
   }
   double **Rhihihi_h = new double*[dim];
   double **Rlohihi_h = new double*[dim];
   double **Rhilohi_h = new double*[dim];
   double **Rlolohi_h = new double*[dim];
   double **Rhihilo_h = new double*[dim];
   double **Rlohilo_h = new double*[dim];
   double **Rhilolo_h = new double*[dim];
   double **Rlololo_h = new double*[dim];
   double **Rhihihi_d = new double*[dim];
   double **Rlohihi_d = new double*[dim];
   double **Rhilohi_d = new double*[dim];
   double **Rlolohi_d = new double*[dim];
   double **Rhihilo_d = new double*[dim];
   double **Rlohilo_d = new double*[dim];
   double **Rhilolo_d = new double*[dim];
   double **Rlololo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhihihi_h[i] = new double[dim];
      Rlohihi_h[i] = new double[dim];
      Rhilohi_h[i] = new double[dim];
      Rlolohi_h[i] = new double[dim];
      Rhihilo_h[i] = new double[dim];
      Rlohilo_h[i] = new double[dim];
      Rhilolo_h[i] = new double[dim];
      Rlololo_h[i] = new double[dim];
      Rhihihi_d[i] = new double[dim];
      Rlohihi_d[i] = new double[dim];
      Rhilohi_d[i] = new double[dim];
      Rlolohi_d[i] = new double[dim];
      Rhihilo_d[i] = new double[dim];
      Rlohilo_d[i] = new double[dim];
      Rhilolo_d[i] = new double[dim];
      Rlololo_d[i] = new double[dim];
   }
/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl8_unit_series_vector
      (dim,deg,inputhihihi_h,inputlohihi_h,inputhilohi_h,inputlolohi_h,
               inputhihilo_h,inputlohilo_h,inputhilolo_h,inputlololo_h);
   for(int i=0; i<dim; i++)
      for(int j=0; j<degp1; j++)
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

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
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
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;
/*
      dbl8_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihihi,cfflohihi,cffhilohi,cfflolohi,
          cffhihilo,cfflohilo,cffhilolo,cfflololo,
          acchihihi,acclohihi,acchilohi,acclolohi,
          acchihilo,acclohilo,acchilolo,acclololo,
          inputhihihi, inputlohihi, inputhilohi, inputlolohi,
          inputhihilo, inputlohilo, inputhilolo, inputlololo,
          outputhihihi,outputlohihi,outputhilohi,outputlolohi,
          outputhihilo,outputlohilo,outputhilolo,outputlololo,
          funvalhihihi,funvallohihi,funvalhilohi,funvallolohi,
          funvalhihilo,funvallohilo,funvalhilolo,funvallololo,
          jacvalhihihi,jacvallohihi,jacvalhilohi,jacvallolohi,
          jacvalhihilo,jacvallohilo,jacvalhilolo,jacvallololo,
          rhshihihi,rhslohihi,rhshilohi,rhslolohi,
          rhshihilo,rhslohilo,rhshilolo,rhslololo,
          solhihihi,sollohihi,solhilohi,sollolohi,
          solhihilo,sollohilo,solhilolo,sollololo,
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          workrhshihihi,workrhslohihi,workrhshilohi,workrhslolohi,
          workrhshihilo,workrhslohilo,workrhshilolo,workrhslololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          ipvt,vrblvl);
 */
      dbl8_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
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
          workmathihihi,workmatlohihi,workmathilohi,workmatlolohi,
          workmathihilo,workmatlohilo,workmathilolo,workmatlololo,
          workvechihihi,workveclohihi,workvechilohi,workveclolohi,
          workvechihilo,workveclohilo,workvechilolo,workveclololo,
          resvechihihi,resveclohihi,resvechilohi,resveclolohi,
          resvechihilo,resveclohilo,resvechilolo,resveclololo,
          &resmaxhihihi,&resmaxlohihi,&resmaxhilohi,&resmaxlolohi,
          &resmaxhihilo,&resmaxlohilo,&resmaxhilolo,&resmaxlololo,
          vrblvl,mode);
   }
   return 0;
}
