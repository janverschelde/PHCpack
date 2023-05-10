// Collects the addition jobs to evaluate and differentiate
// one polynomial in several variables.

#include <iostream>
#include <cstdlib>
#include "addition_jobs.h"
#include "prompt_test_supports.h"

using namespace std;

int main ( void )
{
   int seedused,dim,nva,nbr;

   prompt_testpoly_dimensions(&seedused,&dim,&nva,&nbr);

   int *nvr = new int[nbr]; // number of variables in each monomial
   int **idx = new int*[nbr];  // indices of variables in monomials

   make_test_supports(dim,nva,nbr,nvr,idx);

   AdditionJobs jobs(dim,nbr); // initialize with the number of monomials

   jobs.make(nbr,nvr,idx,true);

   for(int k=0; k<jobs.get_depth(); k++)
   {
      cout << "jobs at layer " << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
         cout << jobs.get_job(k,i) << endl;
   }
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
   cout << "dimension : " << dim << endl;
   cout << "number of monomials : " << nbr << endl;
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

   cout << "seed used : " << seedused << endl;

   return 0;
}
