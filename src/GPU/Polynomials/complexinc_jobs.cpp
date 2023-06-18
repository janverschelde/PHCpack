// The file complexinc_jobs.cpp defines the methods of the class
// ComplexIncrementJobs, specified in "complexinc_jobs.h".

#include <cstdlib>
#include <iostream>
#include "complexconv_jobs.h"
#include "complexinc_jobs.h"

using namespace std;

ComplexIncrementJobs::ComplexIncrementJobs
 ( ComplexConvolutionJobs cnvjobs, bool verbose )
{
   dimension = cnvjobs.get_dimension();

   for(int i=0; i<dimension; i++) freqlaycnt.push_back(0);

   for(int i=0; i<dimension; i++)
   {
      vector<ComplexIncrementJob> jobvec;
      jobs.push_back(jobvec);
   }
   jobcount = 0;
   laydepth = dimension;

   for(int k=0; k<cnvjobs.get_depth(); k++)
   {
      if(verbose)
         cout << "adding jobs at layer " << k << " :" << endl;

      for(int i=0; i<cnvjobs.get_layer_count(k); i++)
      {
         ComplexConvolutionJob cnvjob = cnvjobs.get_job(k,i);
         const int cnvotp = cnvjob.get_output_type();
         const int cnvmon = cnvjob.get_monomial_index();
         const int cnvoix = cnvjob.get_output_index();

         if(cnvotp < 4)
         {
            ComplexIncrementJob job(cnvmon,cnvotp,cnvoix,true);

            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job
                                         << " : layer " << k << endl;
            jobs[k].push_back(job);
            freqlaycnt[k] = freqlaycnt[k] + 1;
         }
         if((cnvotp > 6) && (cnvotp < 10))
         {
            ComplexIncrementJob job(cnvmon,cnvotp-6,cnvoix,false);

            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job
                                         << " : layer " << k << endl;
            jobs[k].push_back(job);
            freqlaycnt[k] = freqlaycnt[k] + 1;
         }
      }
   }
}

int ComplexIncrementJobs::get_dimension ( void ) const
{
   return dimension;
}

int ComplexIncrementJobs::get_count ( void ) const
{
   return jobcount;
}

int ComplexIncrementJobs::get_layer_count ( int k ) const
{
   if(k >= dimension)
      return 0;
   else
      return freqlaycnt[k];
}

int ComplexIncrementJobs::get_depth ( void ) const
{
   return laydepth;
}

ComplexIncrementJob ComplexIncrementJobs::get_job ( int k, int i ) const
{
   return jobs[k][i];
}

ComplexIncrementJobs::~ComplexIncrementJobs ( void )
{ 
   dimension = 0;
   jobcount = 0;
   laydepth = 0;
}

ComplexIncrementJob ComplexIncrementJobs::make_job
 ( ComplexConvolutionJob cnvjob ) const
{
   const int monidx = cnvjob.get_monomial_index();
   const int outptp = cnvjob.get_output_type();
   const int outpix = cnvjob.get_output_index();

   if(outptp < 4)
   {
      ComplexIncrementJob job(monidx,outptp,outpix,true);
      return job;
   }
   else
   {
      ComplexIncrementJob job(monidx,outptp-6,outpix,false);
      return job;
   }
}
