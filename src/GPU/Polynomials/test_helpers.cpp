/* The file test_helpers.cpp contains the definitions of functions
 * with prototypes in test_helpers.h. */

#include <iostream>
#include <iomanip>
#include "test_helpers.h"
#include "write_job_counts.h"

void make_all_jobs
 ( int dim, int nbr, int *nvr, int **idx,
   ConvolutionJobs *cnvjobs, AdditionJobs *addjobs, bool verbose )
{
   cnvjobs->make(nbr,nvr,idx,verbose);

   if(verbose)
   {
      write_convolution_counts(*cnvjobs);

      for(int k=0; k<cnvjobs->get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<cnvjobs->get_layer_count(k); i++)
            cout << cnvjobs->get_job(k,i) << endl;
      }
      cout << endl;
   }
   addjobs->make(nbr,nvr,idx,verbose);

   if(verbose)
   {
      cout << "The differential indices :" << endl;
      for(int i=0; i<dim; i++)
      {
         cout << "variable " << i << " :";
         for(int j=0; j<=addjobs->get_differential_count(i); j++)
            cout << " " << addjobs->get_differential_index(i,j);
         cout << endl;
      }
      write_addition_counts(*addjobs);
   
      for(int k=0; k<addjobs->get_depth(); k++)
      {
         cout << "jobs at layer " << k << " :" << endl;
         for(int i=0; i<addjobs->get_layer_count(k); i++)
            cout << addjobs->get_job(k,i) << endl;
      }
   }
}

void write_jobs_report
 ( int dim, int nva, int nbr, int deg,
   ConvolutionJobs cnvjobs, AdditionJobs addjobs )
{
   cout << "dimension : " << dim << endl;
   if(nva > 0) cout << "number of variables per monomial : " << nva << endl;
   cout << "number of monomials : " << nbr << endl;
   write_convolution_counts(cnvjobs);
   write_addition_counts(addjobs);
   write_operation_counts(deg,cnvjobs,addjobs);
}

void write_CPU_timings ( double lapsec1, double lapsec2 )
{
   cout << fixed << setprecision(3);
   cout << "Elapsed CPU time (Linux), Wall time (Windows) : " << endl;
   cout << "  (1) without jobs : " << lapsec1 << " seconds," << endl;
   cout << "  (2) cnv/add jobs : " << lapsec2 << " seconds." << endl;
   cout << scientific << setprecision(16);
}
