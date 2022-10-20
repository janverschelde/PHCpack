/* Tests the making of a monomial systems in octo double precision. */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "random8_series.h"
#include "unimodular_matrices.h"
#include "dbl8_monomial_systems.h"

using namespace std;

int make_real_system ( int dim, int deg, int **rowsA );
/*
 * DESCRIPTION :
 *   Sets up a monomial system of dim equations,
 *   of series with real coefficients truncated to degree deg,
 *   using the exponents in the rows of rowsA. */

int make_complex_system ( int dim, int deg, int **rowsA );
/*
 * DESCRIPTION :
 *   Sets up a monomial system of dim equations,
 *   of series with complex coefficients truncated to degree deg,
 *   using the exponents in the rows of rowsA. */

int main ( void )
{
   cout << "testing the making of a monomial system ..." << endl;

   const int vrblvl = 2;
   int seed,dim,deg,size,nbritr,cdata;

   cout << "-> give the seed (0 for time) : "; cin >> seed;
   cout << "-> give the number of series : "; cin >> dim;
   while(true)
   {
      cout << "-> give the truncation degree : "; cin >> deg;
      if(deg > 0) break;
      cout << "The degree must be one or larger.  Retry" << endl;
   }
   cout << "-> give the size of the numbers : "; cin >> size;
   cout << "-> give the number of iterations : "; cin >> nbritr;
   cout << "-> on complex data (1 is yes, 0 is no) : "; cin >> cdata;

   if(seed == 0)
      srand(time(NULL));
   else
      srand(seed);

// make the unimodular matrix of the system

   int **rowsA = new int*[dim];  // exponents in the rows
   int *nvr = new int[dim];      // number of variables in each monomial
   int **idx = new int*[dim];    // indexes of variables in each monomial
   int **exp = new int*[dim];    // exponents of the variables
   int *nbrfac = new int[dim];   // number of exponents > 1 in each monomial
   int **expfac = new int*[dim]; // exponents of the common factors

   make_monomial_system
      (dim,size,1,nbritr,nvr,idx,exp,nbrfac,expfac,rowsA,vrblvl);

   int *expsol = new int[dim];
   int sing = exponents_check(dim,rowsA,expsol,vrblvl);

// generate the solution series

   if(cdata == 1)
      return make_complex_system(dim,deg,rowsA);
   else
      return make_real_system(dim,deg,rowsA);
}

int make_real_system ( int dim, int deg, int **rowsA )
{
   double **solhihihi = new double*[dim];
   double **sollohihi = new double*[dim];
   double **solhilohi = new double*[dim];
   double **sollolohi = new double*[dim];
   double **solhihilo = new double*[dim];
   double **sollohilo = new double*[dim];
   double **solhilolo = new double*[dim];
   double **sollololo = new double*[dim];

   const int degp1 = deg+1;

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

   cout << scientific << setprecision(16);

   for(int j=0; j<degp1; j++)
   {
      cout << "coefficients of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "c[" << i << "][" << j << "] : ";     
         cout << solhihihi[i][j] << "  " << sollohihi[i][j] << endl << "  "
              << solhilohi[i][j] << "  " << sollolohi[i][j] << endl << "  "
              << solhihilo[i][j] << "  " << sollohilo[i][j] << endl << "  "
              << solhilolo[i][j] << "  " << sollololo[i][j] << endl;
      }
   }

// compute the right hand sides via evaluation

   double **rhshihihi = new double*[dim];
   double **rhslohihi = new double*[dim];
   double **rhshilohi = new double*[dim];
   double **rhslolohi = new double*[dim];
   double **rhshihilo = new double*[dim];
   double **rhslohilo = new double*[dim];
   double **rhshilolo = new double*[dim];
   double **rhslololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhshihihi[i] = new double[degp1];
      rhslohihi[i] = new double[degp1];
      rhshilohi[i] = new double[degp1];
      rhslolohi[i] = new double[degp1];
      rhshihilo[i] = new double[degp1];
      rhslohilo[i] = new double[degp1];
      rhshilolo[i] = new double[degp1];
      rhslololo[i] = new double[degp1];

      rhshihihi[i][0] = 1.0;     // initialize product to one
      rhslohihi[i][0] = 0.0;
      rhshilohi[i][0] = 0.0; 
      rhslolohi[i][0] = 0.0;
      rhshihilo[i][0] = 0.0; 
      rhslohilo[i][0] = 0.0;
      rhshilolo[i][0] = 0.0; 
      rhslololo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         rhshihihi[i][k] = 0.0; rhslohihi[i][k] = 0.0;
         rhshilohi[i][k] = 0.0; rhslolohi[i][k] = 0.0;
         rhshihilo[i][k] = 0.0; rhslohilo[i][k] = 0.0;
         rhshilolo[i][k] = 0.0; rhslololo[i][k] = 0.0;
      }
   }
   evaluate_real8_monomials
      (dim,deg,rowsA,solhihihi,sollohihi,solhilohi,sollolohi,
                     solhihilo,sollohilo,solhilolo,sollololo,
                     rhshihihi,rhslohihi,rhshilohi,rhslolohi,
                     rhshihilo,rhslohilo,rhshilolo,rhslololo);

   cout << "the evaluated right hand sides :" << endl;

   for(int j=0; j<degp1; j++)
   {
      cout << "coefficients of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "r[" << i << "][" << j << "] : ";     
         cout << rhshihihi[i][j] << "  " << rhslohihi[i][j] << endl << "  "
              << rhshilohi[i][j] << "  " << rhslolohi[i][j] << endl << "  "
              << rhshihilo[i][j] << "  " << rhslohilo[i][j] << endl << "  "
              << rhshilolo[i][j] << "  " << rhslololo[i][j] << endl;
      }
   }
   return 0;
}

int make_complex_system ( int dim, int deg, int **rowsA )
{
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

   const int degp1 = deg+1;

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

   cout << scientific << setprecision(16);

   for(int j=0; j<degp1; j++)
   {
      cout << "coefficients of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "c[" << i << "][" << j << "] : ";     
         cout << solrehihihi[i][j] << "  " << solrelohihi[i][j] << endl << "  "
              << solrehilohi[i][j] << "  " << solrelolohi[i][j] << endl << "  "
              << solrehihilo[i][j] << "  " << solrelohilo[i][j] << endl << "  "
              << solrehilolo[i][j] << "  " << solrelololo[i][j] << endl << "  "
              << solimhihihi[i][j] << "  " << solimlohihi[i][j] << endl << "  "
              << solimhilohi[i][j] << "  " << solimlolohi[i][j] << endl << "  "
              << solimhihilo[i][j] << "  " << solimlohilo[i][j] << endl << "  "
              << solimhilolo[i][j] << "  " << solimlololo[i][j] << endl;
      }
   }

// compute the right hand sides via evaluation

   double **rhsrehihihi = new double*[dim];
   double **rhsrelohihi = new double*[dim];
   double **rhsrehilohi = new double*[dim];
   double **rhsrelolohi = new double*[dim];
   double **rhsrehihilo = new double*[dim];
   double **rhsrelohilo = new double*[dim];
   double **rhsrehilolo = new double*[dim];
   double **rhsrelololo = new double*[dim];
   double **rhsimhihihi = new double*[dim];
   double **rhsimlohihi = new double*[dim];
   double **rhsimhilohi = new double*[dim];
   double **rhsimlolohi = new double*[dim];
   double **rhsimhihilo = new double*[dim];
   double **rhsimlohilo = new double*[dim];
   double **rhsimhilolo = new double*[dim];
   double **rhsimlololo = new double*[dim];

   for(int i=0; i<dim; i++)
   {
      rhsrehihihi[i] = new double[degp1];
      rhsrelohihi[i] = new double[degp1];
      rhsrehilohi[i] = new double[degp1];
      rhsrelolohi[i] = new double[degp1];
      rhsrehihilo[i] = new double[degp1];
      rhsrelohilo[i] = new double[degp1];
      rhsrehilolo[i] = new double[degp1];
      rhsrelololo[i] = new double[degp1];
      rhsimhihihi[i] = new double[degp1];
      rhsimlohihi[i] = new double[degp1];
      rhsimhilohi[i] = new double[degp1];
      rhsimlolohi[i] = new double[degp1];
      rhsimhihilo[i] = new double[degp1];
      rhsimlohilo[i] = new double[degp1];
      rhsimhilolo[i] = new double[degp1];
      rhsimlololo[i] = new double[degp1];

      rhsrehihihi[i][0] = 1.0;     // initialize product to one
      rhsrelohihi[i][0] = 0.0;
      rhsrehilohi[i][0] = 0.0; 
      rhsrelolohi[i][0] = 0.0;
      rhsrehihilo[i][0] = 0.0; 
      rhsrelohilo[i][0] = 0.0;
      rhsrehilolo[i][0] = 0.0; 
      rhsrelololo[i][0] = 0.0;
      rhsimhihihi[i][0] = 0.0;
      rhsimlohihi[i][0] = 0.0;
      rhsimhilohi[i][0] = 0.0;
      rhsimlolohi[i][0] = 0.0;
      rhsimhihilo[i][0] = 0.0;
      rhsimlohilo[i][0] = 0.0;
      rhsimhilolo[i][0] = 0.0;
      rhsimlololo[i][0] = 0.0;

      for(int k=1; k<degp1; k++)
      {
         rhsrehihihi[i][k] = 0.0; rhsrelohihi[i][k] = 0.0;
         rhsrehilohi[i][k] = 0.0; rhsrelolohi[i][k] = 0.0;
         rhsrehihilo[i][k] = 0.0; rhsrelohilo[i][k] = 0.0;
         rhsrehilolo[i][k] = 0.0; rhsrelololo[i][k] = 0.0;
         rhsimhihihi[i][k] = 0.0; rhsimlohihi[i][k] = 0.0;
         rhsimhilohi[i][k] = 0.0; rhsimlolohi[i][k] = 0.0;
         rhsimhihilo[i][k] = 0.0; rhsimlohilo[i][k] = 0.0;
         rhsimhilolo[i][k] = 0.0; rhsimlololo[i][k] = 0.0;
      }
   }
   evaluate_complex8_monomials
      (dim,deg,rowsA,
       solrehihihi,solrelohihi,solrehilohi,solrelolohi,
       solrehihilo,solrelohilo,solrehilolo,solrelololo,
       solimhihihi,solimlohihi,solimhilohi,solimlolohi,
       solimhihilo,solimlohilo,solimhilolo,solimlololo,
       rhsrehihihi,rhsrelohihi,rhsrehilohi,rhsrelolohi,
       rhsrehihilo,rhsrelohilo,rhsrehilolo,rhsrelololo,
       rhsimhihihi,rhsimlohihi,rhsimhilohi,rhsimlolohi,
       rhsimhihilo,rhsimlohilo,rhsimhilolo,rhsimlololo);

   cout << "the evaluated right hand sides :" << endl;

   for(int j=0; j<degp1; j++)
   {
      cout << "coefficients of degree " << j << " :" << endl;

      for(int i=0; i<dim; i++)
      {
         cout << "r[" << i << "][" << j << "] : ";     
         cout << rhsrehihihi[i][j] << "  " << rhsrelohihi[i][j] << endl << "  "
              << rhsrehilohi[i][j] << "  " << rhsrelolohi[i][j] << endl << "  "
              << rhsrehihilo[i][j] << "  " << rhsrelohilo[i][j] << endl << "  "
              << rhsrehilolo[i][j] << "  " << rhsrelololo[i][j] << endl << "  "
              << rhsimhihihi[i][j] << "  " << rhsimlohihi[i][j] << endl << "  "
              << rhsimhilohi[i][j] << "  " << rhsimlolohi[i][j] << endl << "  "
              << rhsimhihilo[i][j] << "  " << rhsimlohilo[i][j] << endl << "  "
              << rhsimhilolo[i][j] << "  " << rhsimlololo[i][j] << endl;
      }
   }
   return 0;
}
