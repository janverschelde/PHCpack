// The file addition_jobs.cpp defines the methods of the class AdditionJobs,
// specified in "addition_jobs.h".

#include <cstdlib>
#include <iostream>
#include "addition_jobs.h"

AdditionJobs::AdditionJobs ( int dim )
{
   dimension = dim;
   jobcount = 0;
   laydepth = 0;
}

void AdditionJobs::make ( int nbr, int *nvr, bool verbose )
{
   freqlaycnt = new int[dimension];
   for(int i=0; i<dimension; i++) freqlaycnt[i] = 0;

   for(int i=0; i<dimension; i++)
   {
      vector<AdditionJob> jobvec;
      jobs.push_back(jobvec);
   }
   if(verbose) cout << "layer 0 : " << endl;
   {
      jobcount = jobcount + 1; freqlaycnt[0] = freqlaycnt[0] + 1;
      AdditionJob job(1,0,-1,nvr[0]-1,-1);
      if(verbose) cout << jobcount << " : " << job
                                   << " : layer 0" << endl;
      jobs[0].push_back(job);
   }
   laydepth = 1; // we have one layer
   if(nbr > 1)
   {
      for(int i=1; i<nbr-1; i=i+2) 
      {
         jobcount = jobcount + 1; freqlaycnt[0] = freqlaycnt[0] + 1;
         AdditionJob job(1,i+1,i,nvr[i+1]-1,nvr[i]-1);
         if(verbose) cout << jobcount << " : " << job
                                      << " : layer 0" << endl;
         jobs[0].push_back(job);
      }
      int stride = 2;
      int laycnt = 1;
      int istart = 0;

      while(stride < nbr)
      {
         if(verbose) cout << "layer " << laycnt << " :" << endl;
    
         for(int i=istart; i<nbr-stride; i=i+2*stride) 
         {
            jobcount = jobcount + 1;
            freqlaycnt[laycnt] = freqlaycnt[laycnt] + 1;
            AdditionJob job(1,i+stride,i,nvr[i+stride]-1,nvr[i]-1);
            if(verbose) cout << jobcount << " : " << job
                             << " : layer " << laycnt << endl;
            jobs[laycnt].push_back(job);
         }
         laycnt = laycnt + 1;
         istart = istart + stride;
         stride = 2*stride;
      }
      laydepth = laycnt;
   }
}

int AdditionJobs::get_dimension ( void ) const
{
   return dimension;
}

int AdditionJobs::get_count ( void ) const
{
   return jobcount;
}

int AdditionJobs::get_layer_count ( int k ) const
{
   if(k >= dimension)
      return 0;
   else
      return freqlaycnt[k];
}

int AdditionJobs::get_depth ( void ) const
{
   return laydepth;
}

AdditionJob AdditionJobs::get_job ( int k, int i ) const
{
   return jobs[k][i];
}

AdditionJobs::~AdditionJobs ( void )
{
   dimension = 0;
}
