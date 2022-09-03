/* Tests operations on Newton for series in quad double precision */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "unimodular_matrices.h"
#include "dbl4_newton_testers.h"

using namespace std;

int main ( void )
{
   srand(time(NULL));

   cout << "testing Newton ..." << endl;

   int dim,deg,size,posvals,vrblvl,nbritr,nbsteps,szt,nbt,mode;
   prompt_dimensions(&dim,&deg,&size,&posvals,&vrblvl,&nbritr,&nbsteps);

   szt = dim; nbt = 1; mode = 1; // stub values

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
   // allocate memory for coefficients and the output
   double *acchihi = new double[degp1]; // accumulated power series
   double *acclohi = new double[degp1];
   double *acchilo = new double[degp1];
   double *acclolo = new double[degp1];
   double **cffhihi = new double*[dim]; // the coefficients of monomials
   double **cfflohi = new double*[dim];
   double **cffhilo = new double*[dim];
   double **cfflolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      cffhihi[i] = new double[degp1];
      cfflohi[i] = new double[degp1];
      cffhilo[i] = new double[degp1];
      cfflolo[i] = new double[degp1];
   }
   double ***outputhihi_h = new double**[dim];
   double ***outputlohi_h = new double**[dim];
   double ***outputhilo_h = new double**[dim];
   double ***outputlolo_h = new double**[dim];
   double ***outputhihi_d = new double**[dim];
   double ***outputlohi_d = new double**[dim];
   double ***outputhilo_d = new double**[dim];
   double ***outputlolo_d = new double**[dim];

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
   // The function values are power series truncated at degree deg.
   double **funvalhihi_h = new double*[dim];
   double **funvallohi_h = new double*[dim];
   double **funvalhilo_h = new double*[dim];
   double **funvallolo_h = new double*[dim];
   double **funvalhihi_d = new double*[dim];
   double **funvallohi_d = new double*[dim];
   double **funvalhilo_d = new double*[dim];
   double **funvallolo_d = new double*[dim];

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
   double ***jacvalhihi_h = new double**[degp1];
   double ***jacvallohi_h = new double**[degp1];
   double ***jacvalhilo_h = new double**[degp1];
   double ***jacvallolo_h = new double**[degp1];
   double ***jacvalhihi_d = new double**[degp1];
   double ***jacvallohi_d = new double**[degp1];
   double ***jacvalhilo_d = new double**[degp1];
   double ***jacvallolo_d = new double**[degp1];

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
/*
 * 3. allocate space to solve the linearized power series system
 */
   // The solution x(t) to jacval(t)*x(t) = -funval(t) in linearized
   // format is a series truncated at degree deg, with as coefficients
   // arrays of dimension dim.
   double **solhihi_h = new double*[degp1];
   double **sollohi_h = new double*[degp1];
   double **solhilo_h = new double*[degp1];
   double **sollolo_h = new double*[degp1];
   double **solhihi_d = new double*[degp1];
   double **sollohi_d = new double*[degp1];
   double **solhilo_d = new double*[degp1];
   double **sollolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      solhihi_h[i] = new double[dim];
      sollohi_h[i] = new double[dim];
      solhilo_h[i] = new double[dim];
      sollolo_h[i] = new double[dim];
      solhihi_d[i] = new double[dim];
      sollohi_d[i] = new double[dim];
      solhilo_d[i] = new double[dim];
      sollolo_d[i] = new double[dim];
   }
   // The right hand side -funval(t) in linearized format is a series
   // truncated at degree deg, with arrays of dimension dim as coefficients.
   double **rhshihi_h = new double*[degp1];
   double **rhslohi_h = new double*[degp1];
   double **rhshilo_h = new double*[degp1];
   double **rhslolo_h = new double*[degp1];
   double **rhshihi_d = new double*[degp1];
   double **rhslohi_d = new double*[degp1];
   double **rhshilo_d = new double*[degp1];
   double **rhslolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      rhshihi_h[i] = new double[dim];
      rhslohi_h[i] = new double[dim];
      rhshilo_h[i] = new double[dim];
      rhslolo_h[i] = new double[dim];
      rhshihi_d[i] = new double[dim];
      rhslohi_d[i] = new double[dim];
      rhshilo_d[i] = new double[dim];
      rhslolo_d[i] = new double[dim];
   }
   // Copy the rhs vector into work space for inplace solver.
   double **urhshihi_h = new double*[degp1];
   double **urhslohi_h = new double*[degp1];
   double **urhshilo_h = new double*[degp1];
   double **urhslolo_h = new double*[degp1];
   double **urhshihi_d = new double*[degp1];
   double **urhslohi_d = new double*[degp1];
   double **urhshilo_d = new double*[degp1];
   double **urhslolo_d = new double*[degp1];

   for(int i=0; i<degp1; i++)
   {
      urhshihi_h[i] = new double[dim];
      urhslohi_h[i] = new double[dim];
      urhshilo_h[i] = new double[dim];
      urhslolo_h[i] = new double[dim];
      urhshihi_d[i] = new double[dim];
      urhslohi_d[i] = new double[dim];
      urhshilo_d[i] = new double[dim];
      urhslolo_d[i] = new double[dim];
   }
   // Allocate work space for the inplace LU solver.
   double **workmathihi = new double*[dim];
   double **workmatlohi = new double*[dim];
   double **workmathilo = new double*[dim];
   double **workmatlolo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      workmathihi[i] = new double[dim];
      workmatlohi[i] = new double[dim];
      workmathilo[i] = new double[dim];
      workmatlolo[i] = new double[dim];
   }
   int *ipvt = new int[dim];
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
   double **Qhihi_h = new double*[dim];
   double **Qlohi_h = new double*[dim];
   double **Qhilo_h = new double*[dim];
   double **Qlolo_h = new double*[dim];
   double **Qhihi_d = new double*[dim];
   double **Qlohi_d = new double*[dim];
   double **Qhilo_d = new double*[dim];
   double **Qlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Qhihi_h[i] = new double[dim];
      Qlohi_h[i] = new double[dim];
      Qhilo_h[i] = new double[dim];
      Qlolo_h[i] = new double[dim];
      Qhihi_d[i] = new double[dim];
      Qlohi_d[i] = new double[dim];
      Qhilo_d[i] = new double[dim];
      Qlolo_d[i] = new double[dim];
   }
   double **Rhihi_h = new double*[dim];
   double **Rlohi_h = new double*[dim];
   double **Rhilo_h = new double*[dim];
   double **Rlolo_h = new double*[dim];
   double **Rhihi_d = new double*[dim];
   double **Rlohi_d = new double*[dim];
   double **Rhilo_d = new double*[dim];
   double **Rlolo_d = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      Rhihi_h[i] = new double[dim];
      Rlohi_h[i] = new double[dim];
      Rhilo_h[i] = new double[dim];
      Rlolo_h[i] = new double[dim];
      Rhihi_d[i] = new double[dim];
      Rlohi_d[i] = new double[dim];
      Rhilo_d[i] = new double[dim];
      Rlolo_d[i] = new double[dim];
   }
/*
 * 4. initialize input, coefficient, evaluate, differentiate, and solve
 */
   // Define the initial input, a vector of ones.
   dbl4_unit_series_vector
      (dim,deg,inputhihi_h,inputlohi_h,inputhilo_h,inputlolo_h);

   if(vrblvl > 0)
   {
      cout << scientific << setprecision(16);
      cout << "The leading coefficients of the input series :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << i << " : " << inputhihi_h[i][0] << "  "
                            << inputlohi_h[i][0] << endl;
         cout << "     " << inputhilo_h[i][0] << "  "
                         << inputlolo_h[i][0] << endl;
      }
   }
   for(int step=0; step<nbsteps; step++)
   {
      cout << "step " << step << " ..." << endl;
/*
      dbl4_newton_lustep
         (dim,deg,nvr,idx,exp,nbrfac,expfac,
          cffhihi,cfflohi,cffhilo,cfflolo,acchihi,acclohi,acchilo,acclolo,
           inputhihi, inputlohi, inputhilo, inputlolo,
          outputhihi,outputlohi,outputhilo,outputlolo,
          funvalhihi,funvallohi,funvalhilo,funvallolo,
          jacvalhihi,jacvallohi,jacvalhilo,jacvallolo,
          rhshihi,rhslohi,rhshilo,rhslolo,solhihi,sollohi,solhilo,sollolo,
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          workvechihi,workveclohi,workvechilo,workveclolo,
          workrhshihi,workrhslohi,workrhshilo,workrhslolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,ipvt,vrblvl);
 */
      dbl4_newton_qrstep
         (szt,nbt,dim,deg,nvr,idx,exp,nbrfac,expfac,
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
          workmathihi,workmatlohi,workmathilo,workmatlolo,
          workvechihi,workveclohi,workvechilo,workveclolo,
          resvechihi,resveclohi,resvechilo,resveclolo,
          &resmaxhihi,&resmaxlohi,&resmaxhilo,&resmaxlolo,vrblvl,mode);
   }
   return 0;
}
