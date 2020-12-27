// The file addition_jobs.cpp defines the methods of the class AdditionJobs,
// specified in "addition_jobs.h".

#include <cstdlib>
#include <iostream>
#include "addition_jobs.h"

AdditionJobs::AdditionJobs ( int nbr )
{
   dimension = nbr;
   jobcount = 0;
   laydepth = 0;
}

void AdditionJobs::recursive_start ( int nbr, int *level, int *stride )
{
   int lvl = 0;
   int pwr = 1; // pwr = 2**lvl

   while(pwr <= nbr)
   {
      lvl = lvl+1;
      pwr = pwr*2;
   }
   *level = lvl-1;
   *stride = pwr/2;
}

void AdditionJobs::recursive_make
 ( int level, int stride, int nbr, int *nvr, bool verbose )
{
   const int ix1 = nbr - 1;
   const int ix2 = ix1 - stride;

   if(ix2 >= 0)
   {
      AdditionJob job(1,ix1,ix2,nvr[ix1]-1,nvr[ix2]-1);
      if(verbose) cout << "adding " << job << " to layer " << level << endl;
      jobs[level].push_back(job);
      freqlaycnt[level] = freqlaycnt[level] + 1;
   }
   else if(ix2 == -1)
   {
      AdditionJob job(1,ix1,-1,nvr[ix1]-1,-1);
      if(verbose) cout << "adding " << job << " to layer " << level << endl;
      jobs[level].push_back(job);
      freqlaycnt[level] = freqlaycnt[level] + 1;
   }
   if(level > 0)
   {
      recursive_make(level-1,stride/2,nbr,nvr,verbose);
      if(nbr > stride)
         recursive_make(level-1,stride/2,nbr-stride,nvr,verbose);
   }
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
   int level,stride;

   recursive_start(nbr,&level,&stride);
   recursive_make(level,stride,nbr,nvr,verbose);

   laydepth = level+1;
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
