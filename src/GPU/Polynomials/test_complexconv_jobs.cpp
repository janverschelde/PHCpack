// Tests the definition of convolution and increment jobs
// to evaluate and differentiate one polynomial 
// in several variables at complex numbers.

#include <iostream>
#include <cstdlib>
#include "prompt_test_supports.h"
#include "complexconv_jobs.h"
#include "complexinc_job.h"
#include "complexinc_jobs.h"

using namespace std;

int main ( void )
{
   int seedused,dim,nva,nbr;

   prompt_testpoly_dimensions(&seedused,&dim,&nva,&nbr);

   int *nvr = new int[nbr]; // number of variables in each monomial
   int **idx = new int*[nbr];  // indices of variables in monomials

   make_test_supports(dim,nva,nbr,nvr,idx);

   ComplexConvolutionJobs cnvjobs(dim);

   cnvjobs.make(nbr,nvr,idx,true);

   ComplexIncrementJobs incjobs(cnvjobs,true);

   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      cout << "jobs at layer " << k << " :" << endl;
      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
         cout << cnvjobs.get_job(k,i) << endl;
      for(int i=0; i<incjobs.get_layer_count(k); i++)
         cout << incjobs.get_job(k,i) << endl;
   }
   cout << "dimension : " << dim << endl;
   cout << "number of monomials : " << nbr << endl;
   cout << "number of convolution jobs : " << cnvjobs.get_count() << endl;
   cout << "number of increment jobs : " << incjobs.get_count() << endl;
   cout << "number of layers : " << cnvjobs.get_depth() << endl;
   cout << "frequency of layer counts :" << endl;
   int checksum = 0;
   for(int i=0; i<cnvjobs.get_depth(); i++)
   {
      cout << i << " : " << cnvjobs.get_layer_count(i)
                << " + " << incjobs.get_layer_count(i) << endl;
      checksum = checksum + cnvjobs.get_layer_count(i)
                          + incjobs.get_layer_count(i); 
   }
   cout << "total layer count sum : " << checksum << endl;
   cout << "number of convolutions jobs : " << cnvjobs.get_count() << endl;
   cout << "number of increment jobs : " << incjobs.get_count() << endl;

   cout << "seed used : " << seedused << endl;

   return 0;
}
