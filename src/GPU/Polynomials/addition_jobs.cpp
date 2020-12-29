// The file addition_jobs.cpp defines the methods of the class AdditionJobs,
// specified in "addition_jobs.h".

#include <cstdlib>
#include <iostream>
#include "addition_jobs.h"

AdditionJobs::AdditionJobs ( int dim, int nbr )
{
   nbrvar = dim;
   nbrmon = nbr;
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
      jobcount = jobcount + 1;
      freqlaycnt[level] = freqlaycnt[level] + 1;
   }
   else if(ix2 == -1)
   {
      AdditionJob job(1,ix1,-1,nvr[ix1]-1,-1);
      if(verbose) cout << "adding " << job << " to layer " << level << endl;
      jobs[level].push_back(job);
      jobcount = jobcount + 1;
      freqlaycnt[level] = freqlaycnt[level] + 1;
   }
   if(level > 0)
   {
      recursive_make(level-1,stride/2,nbr,nvr,verbose);
      if(nbr > stride)
         recursive_make(level-1,stride/2,nbr-stride,nvr,verbose);
   }
}

void AdditionJobs::recursive_first_make
 ( int level, int stride, int nbr, int *nvr, bool verbose )
{
   const int ix0 = nbr - stride;
   const int ix1 = difidx[0][nbr];
   const int ix2 = difidx[0][ix0];

   if(ix0 > 0)
   {
      AdditionJob job(2,ix1,ix2,nvr[ix1]-2,nvr[ix2]-2);
      if(verbose) cout << "adding " << job << " to layer " << level << endl;
      jobs[level].push_back(job);
      jobcount = jobcount + 1;
      freqlaycnt[level] = freqlaycnt[level] + 1;
   }
   else if((ix0 == 0) and (ix2 != -1))
   {
      AdditionJob job(2,ix1,-1,nvr[ix1]-2,ix2);
      if(verbose) cout << "adding " << job << " to layer " << level << endl;
      jobs[level].push_back(job);
      jobcount = jobcount + 1;
      freqlaycnt[level] = freqlaycnt[level] + 1;
   }
   if(level > 0)
   {
      recursive_first_make(level-1,stride/2,nbr,nvr,verbose);
      if(nbr > stride)
         recursive_first_make(level-1,stride/2,nbr-stride,nvr,verbose);
   }
}

void AdditionJobs::differential_index_count
 ( int dim, int nbr, int *nvr, int **idx, int *cnt, bool verbose )
{
   for(int i=0; i<dim; i++)
   {
      if(verbose) cout << "Variable " <<  i << " occurs in monomials"; 

      cnt[i] = 0;
      for(int j=0; j<nbr; j++)
      {
         for(int k=0; k<nvr[j]; k++)
            if(idx[j][k] == i)
            {
               if(verbose)
               {
                  cout << " " << j;
                  if(nvr[j] == 1) cout << "(cff!)";
               }
               cnt[i] = cnt[i] + 1; break;
            }
      }
      if(verbose) cout << endl;
   }
}

void AdditionJobs::make_differential_indices
 ( int dim, int nbr, int *nvr, int **idx, int *cnt, int **difidx,
   bool verbose )
{
   int pos;

   for(int i=0; i<dim; i++)
   {
      if(verbose) cout << "Variable " <<  i << " occurs in monomials"; 
      
      difidx[i][0] = -1;
      pos = 1;
      for(int j=0; j<nbr; j++)
      {
         for(int k=0; k<nvr[j]; k++)
            if(idx[j][k] == i)
            {
               if(verbose)
               {
                  cout << " " << j;
                  if(nvr[j] == 1) cout << "(cff!)";
               }
               if(nvr[j] == 1)
               {
                  difidx[i][0] = j;
                  cnt[i] = cnt[i] - 1;
               }
               else
               {
                  difidx[i][pos++] = j;
               }
               break;
            }
      }
      if(verbose) cout << endl;
   }
}

void AdditionJobs::make ( int nbr, int *nvr, int **idx, bool verbose )
{
   freqlaycnt = new int[nbrmon];
   for(int i=0; i<nbrmon; i++) freqlaycnt[i] = 0;

   difcnt = new int[nbrmon];
   differential_index_count(nbrvar,nbr,nvr,idx,difcnt,verbose);

   difidx = new int*[nbrvar];
   for(int i=0; i<nbrvar; i++) difidx[i] = new int[difcnt[i]+1];
   make_differential_indices(nbrvar,nbr,nvr,idx,difcnt,difidx,verbose);

   for(int i=0; i<nbrmon; i++)
   {
      vector<AdditionJob> jobvec;
      jobs.push_back(jobvec);
   }
   int level,stride;

   if(verbose) cout << "-> adding jobs for the value ..." << endl;
   recursive_start(nbr,&level,&stride);
   recursive_make(level,stride,nbr,nvr,verbose);

   laydepth = level+1;

   if(verbose) cout << "-> adding jobs for the first derivative ..." << endl;
   recursive_start(difcnt[0],&level,&stride);
   recursive_first_make(level,stride,difcnt[0],nvr,verbose);
}

int AdditionJobs::get_number_of_variables ( void ) const
{
   return nbrvar;
}

int AdditionJobs::get_number_of_monomials ( void ) const
{
   return nbrmon;
}

int AdditionJobs::get_count ( void ) const
{
   return jobcount;
}

int AdditionJobs::get_layer_count ( int k ) const
{
   if(k >= nbrmon)
      return 0;
   else
      return freqlaycnt[k];
}

int AdditionJobs::get_depth ( void ) const
{
   return laydepth;
}

int AdditionJobs::get_differential_count ( int k ) const
{
   if(k >= nbrvar)
      return 0;
   else
      return difcnt[k];
}

int AdditionJobs::get_differential_index ( int k, int i ) const
{
   if(k >= nbrvar)
      return -1;
   else if(i > difcnt[k])
      return -1;
   else
      return difidx[k][i];
}

AdditionJob AdditionJobs::get_job ( int k, int i ) const
{
   return jobs[k][i];
}

AdditionJobs::~AdditionJobs ( void )
{
   nbrvar = 0;
   nbrmon = 0;
}
