// The file job_coordinates.cpp defines the functions specified
// in job_coordinates.h.

#include <iostream>
#include "job_coordinates.h"

using namespace std;

int coefficient_count ( int dim, int nbr, int deg, int *nvr )
{
   int count = 1 + nbr + dim;

   for(int i=0; i<nbr; i++)
   {
      count = count + nvr[i];
      if(nvr[i] == 2)
         count = count + 1;
      else if(nvr[i] > 2)
         count = count + 2*(nvr[i]-2);
   }
   count = count*(deg+1);

   return count;
}

void coefficient_indices
 ( int dim, int nbr, int deg, int *nvr,
   int *fsums, int *bsums, int *csums,
   int *fstart, int *bstart, int *cstart )
{
   fsums[0] = nvr[0]; bsums[0] = 0; csums[0] = 0;

   if(nvr[0] == 2)
   {
      bsums[0] = 1;
   }
   else if(nvr[0] > 2)
   {
      bsums[0] = nvr[0] - 2;
      csums[0] = nvr[0] - 2;
   }
   for(int i=1; i<nbr; i++)
   {
      fsums[i] = fsums[i-1] + nvr[i];
      if(nvr[i] < 2)
      {
         bsums[i] = bsums[i-1];
         csums[i] = csums[i-1];
      }
      else if(nvr[i] == 2)
      {
         bsums[i] = bsums[i-1] + 1;
         csums[i] = csums[i-1];
      }
      else // nvr[i] > 2
      {
         bsums[i] = bsums[i-1] + nvr[i] - 2;
         csums[i] = csums[i-1] + nvr[i] - 2;
      }
   }
   fstart[0] = (1+nbr+dim)*(deg+1);
   for(int i=1; i<nbr; i++) fstart[i] = fstart[0] + fsums[i-1]*(deg+1);

   bstart[0] = fstart[0] + fsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++) bstart[i] = bstart[0] + bsums[i-1]*(deg+1);

   cstart[0] = bstart[0] + bsums[nbr-1]*(deg+1);
   for(int i=1; i<nbr; i++) cstart[i] = cstart[0] + csums[i-1]*(deg+1);
}

void convjob_indices
 ( ConvolutionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int monidx = job.get_monomial_index();
   const int jobinp1tp = job.get_first_type();
   const int jobinp1ix = job.get_first_input();
   const int jobinp2tp = job.get_second_type();
   const int jobinp2ix = job.get_second_input();
   const int joboutptp = job.get_output_type();
   const int joboutidx = job.get_output_index();
   const int deg1 = deg+1;

   if(verbose)
   {
      cout << "  in1 type : " << jobinp1tp << ", idx : " << jobinp1ix;
      cout << "  in2 type : " << jobinp2tp << ", idx : " << jobinp2ix;
      cout << "  out type : " << joboutptp << ", idx : " << joboutidx;
      cout << endl;
   }
   if(jobinp1tp < 0)           // first input is coefficient
      *inp1ix = (1 + monidx)*deg1;
   else if(jobinp1tp == 0)     // first input is input series
      *inp1ix = (1 + nbr + jobinp1ix)*deg1;
   else if(jobinp1tp == 1)     // first input is forward product
      *inp1ix = fstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 2)     // first input is backward product
      *inp1ix = bstart[monidx] + jobinp1ix*deg1;
   else if(jobinp1tp == 3)     // first input is cross product
      *inp1ix = cstart[monidx] + jobinp1ix*deg1;

   if(jobinp2tp < 0)           // second input is coefficient
      *inp2ix = (1 + monidx)*deg1;
   else if(jobinp2tp == 0)     // second input is input series
      *inp2ix = (1 + nbr + jobinp2ix)*deg1;
   else if(jobinp2tp == 1)     // second input is forward product
      *inp2ix = fstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 2)     // second input is backward product
      *inp2ix = bstart[monidx] + jobinp2ix*deg1;
   else if(jobinp2tp == 3)     // second input is cross product
      *outidx = cstart[monidx] + jobinp2ix*deg1;

   if(joboutptp == 1)          // output is forward product
      *outidx = fstart[monidx] + joboutidx*deg1;
   else if(joboutptp == 2)    // output is backward product
   {
      if(joboutidx < nvr[monidx]-2) // last backward product is special ...
         *outidx = bstart[monidx] + joboutidx*deg1;
      else
      {
         if(nvr[monidx] == 2)
            *outidx = bstart[monidx];
         else
            *outidx = bstart[monidx] + (joboutidx-1)*deg1;

      }
   }
   else if(joboutptp == 3)    // output is cross product
      *outidx = cstart[monidx] + joboutidx*deg1;

   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void convjobs_coordinates
 ( ConvolutionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{ 
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      convjob_indices(jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
                      dim,nbr,deg,nvr,fstart,bstart,cstart,verbose);
}

void addjob_indices
 ( AdditionJob job, int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   const int adtype = job.get_addition_type();
   const int intype = job.get_increment_type();
   const int updmon = job.get_update_monomial();
   const int updidx = job.get_update_index();
   const int incmon = job.get_increment_monomial();
   const int incidx = job.get_increment_index();
   const int deg1 = deg+1;

   if(verbose)
   {
      cout << "  add type : " << adtype
           << ", mon : " << updmon << ", idx : " << updidx;
      cout << "  inc type : " << intype
           << ", mon : " << incmon << ", idx : " << incidx;
      cout << endl;
   }
   if(adtype == 1)
   {
      *inp1ix = fstart[updmon] + updidx*deg1;
   }
   else if(adtype == 2)  // on GPU, one backward item less
   {
      if(updidx == 0)
         *inp1ix = bstart[updmon];
      else
         *inp1ix = bstart[updmon] + (updidx-1)*deg1;
   }
   else if(adtype == 3)
   {
      *inp1ix = cstart[updmon] + updidx*deg1;
   }
   *outidx = *inp1ix;
   if(incmon < 0)
   {
      if(incidx < 0)
         *inp2ix = 0; // start with constant coefficient
      else
         *inp2ix = (1 + incidx)*deg1;
   }
   else
   {
      if(intype == 1)
      {
         *inp2ix = fstart[incmon] + incidx*deg1;
      }
      else if(intype == 2)
      {                                // on GPU, on backward item less
         if(incidx == 0)
            *inp2ix = bstart[incmon];
         else
            *inp2ix = bstart[incmon] + (incidx-1)*deg1;
      }
      else if(intype == 3)
      {
         *inp2ix = cstart[incmon] + incidx*deg1;
      }
   }
   if(verbose)
   {
      cout << "-> inp1ix : " << *inp1ix
           << ", inp2ix : " << *inp2ix
           << ", outidx : " << *outidx << endl;
   }
}

void addjobs_coordinates
 ( AdditionJobs jobs, int layer,
   int *inp1ix, int *inp2ix, int *outidx,
   int dim, int nbr, int deg, int *nvr,
   int *fstart, int *bstart, int *cstart, bool verbose )
{
   for(int i=0; i<jobs.get_layer_count(layer); i++)
      addjob_indices(jobs.get_job(layer,i),&inp1ix[i],&inp2ix[i],&outidx[i],
                     dim,nbr,deg,nvr,fstart,bstart,cstart,verbose);
}
