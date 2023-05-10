// Collects the convolution jobs to evaluate and differentiate
// one polynomial in several variables.

#include <iostream>
#include <cstdlib>
#include "convolution_jobs.h"
#include "prompt_test_supports.h"

using namespace std;

int main ( void )
{
   int seedused,dim,nva,nbr;

   prompt_testpoly_dimensions(&seedused,&dim,&nva,&nbr);

   int *nvr = new int[nbr]; // number of variables in each monomial
   int **idx = new int*[nbr];  // indices of variables in monomials

   make_test_supports(dim,nva,nbr,nvr,idx);

   ConvolutionJobs jobs(dim);

   jobs.make(nbr,nvr,idx,true);

   for(int k=0; k<jobs.get_depth(); k++)
   {
      cout << "jobs at layer " << k << " :" << endl;
      for(int i=0; i<jobs.get_layer_count(k); i++)
         cout << jobs.get_job(k,i) << endl;
   }
   cout << "dimension : " << dim << endl;
   cout << "number of monomials : " << nbr << endl;
   cout << "number of convolution jobs : " << jobs.get_count() << endl;
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
