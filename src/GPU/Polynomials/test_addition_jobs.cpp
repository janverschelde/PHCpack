// Collects the addition jobs to evaluate and differentiate
// one polynomial in several variables.

#include <iostream>
#include <cstdlib>
#include <ctime>
#include "random_polynomials.h"
#include "addition_job.h"
#include "addition_jobs.h"

using namespace std;

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
   cout << "  seed used : " << seedused << endl;

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
   bool dup = duplicate_supports(dim,nbr,nvr,idx,false);
   if(dup)
      cout << "Duplicate supports found." << endl;
   else
      cout << "No duplicate supports found." << endl;

   AdditionJobs jobs(dim,nbr); // initialize with the number of monomials

   jobs.make(nbr,nvr,idx,true);

   cout << "Index count :";
   for(int i=0; i<dim; i++) cout << " " << jobs.get_differential_count(i);
   cout << endl;
   cout << "The differential indices :" << endl;
   for(int i=0; i<dim; i++)
   {
      cout << "variable " << i << " :";
      for(int j=0; j<=jobs.get_differential_count(i); j++)
         cout << " " << jobs.get_differential_index(i,j);
      cout << endl;
   }
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
