// Collects the addition jobs to evaluate and differentiate
// one polynomial in several variables.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "random_polynomials.h"
#include "addition_job.h"
#include "addition_jobs.h"

using namespace std;

void write_addition_jobs ( int dim, int nbr, int *nvr );
/*
 * DESCRIPTION :
 *   Writes all jobs to add forward, backward, and cross products
 *   to define a reduction tree.
 *
 * ON ENTRY :
 *   dim      total number of variables;
 *   nbr      number of monomials, excluding the constant;
 *   nvr      nbr integers with the number of variables in each monomial,
 *            nvr[k] is the number of variables in monomial k. */

int main ( void )
{
   cout << "Give the seed (0 for time) : ";
   int seed; cin >> seed;

   cout << "Give the dimension : ";
   int dim;  cin >> dim;


   cout << "Give the number of terms : ";
   int nbr; cin >> nbr;

   int seedused;

   if(seed != 0)
   {
      srand(seed);
      seedused = seed;
   }
   else
   {
      const int timevalue = time(NULL); // for a random seed
      srand(timevalue);
      seedused = timevalue;
   }
   const int deg = 0;
   const int pwr = 1;

   double *cst = new double[deg+1]; // constant coefficient series
   double **cff = new double*[nbr]; // coefficient series of terms
   for(int i=0; i<nbr; i++) cff[i] = new double[deg+1];
   int *nvr = new int[nbr]; // number of variables in each monomial

   make_supports(dim,nbr,nvr); // define supports of polynomial

   int **idx = new int*[nbr];  // indices of variables in monomials
   for(int i=0; i<nbr; i++) idx[i] = new int[nvr[i]];
   int **exp = new int*[nbr];  // exponents of the variables
   for(int i=0; i<nbr; i++) exp[i] = new int[nvr[i]];

   bool fail = make_real_polynomial(dim,nbr,pwr,deg,nvr,idx,exp,cst,cff);

   for(int i=0; i<nbr; i++)
   {
      cout << "Indices of monomial " << i << " :";
      for(int j=0; j<nvr[i]; j++) cout << " " << idx[i][j]; cout << endl;
   }
   write_addition_jobs(dim,nbr,nvr);

   AdditionJobs jobs(dim);

   jobs.make(nbr,nvr,true);

   cout << "number of addition jobs : " << jobs.get_count() << endl;
   cout << "number of layers : " << jobs.get_depth() << endl;
   cout << "frequency of layer counts :" << endl;
   int checksum = 0;
   for(int i=0; i<jobs.get_depth(); i++)
   {
      cout << i << " : " << jobs.get_layer_count(i) << endl;
      checksum = checksum + jobs.get_layer_count(i); 
   }
   cout << "layer count sum : " << checksum << endl;

   for(int k=0; k<jobs.get_depth(); k++)
   {
      cout << "jobs at layer " << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
         cout << jobs.get_job(k,i) << endl;
   }
   cout << "seed used : " << seedused << endl;

   return 0;
}

void write_addition_jobs ( int dim, int nbr, int *nvr )
{
   cout << "layer 0 : " << endl;
   {
      AdditionJob job(1,1,-1,0,-1);
      cout << job << endl;
   }
   if(nbr > 1)
   {
      for(int i=1; i<nbr-1; i=i+2) 
      {
         AdditionJob job(1,i+1,i,nvr[i+1]-1,nvr[i]-1);
         cout << job << endl;
      }
      int stride = 2;
      int laycnt = 1;
      int istart = 0;

      while(stride < nbr)
      {
         cout << "layer " << laycnt++ << " :" << endl;
    
         for(int i=istart; i<nbr-stride; i=i+2*stride) 
         {
            AdditionJob job(1,i+stride,i,nvr[i+stride]-1,nvr[i]-1);
            cout << job << endl;
         }

         istart = istart + stride;
         stride = 2*stride;
      }
   }
}
