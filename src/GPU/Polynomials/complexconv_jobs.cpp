// The file complexconv_jobs.cpp defines the methods of the class
// ComplexConvolutionJobs, specified in "complexconv_jobs.h".

#include <cstdlib>
#include <iostream>
#include "complexconv_jobs.h"

using namespace std;

ComplexConvolutionJobs::ComplexConvolutionJobs ( int dim )
{
   dimension = dim;
   jobcount = 0;
   laydepth = 0;
}

void ComplexConvolutionJobs::make_monomial
 ( int nvr, int *idx, int monidx, bool verbose )
{
   int ix1 = idx[0];
   int ix2;
   int layer = 0; // determines the order of execution

   freqlaycnt[layer] = freqlaycnt[layer] + 4; // add 4 jobs to layer
   {
      // coefficient times the first variable to forward f[0]
      // first operand of real part
      ComplexConvolutionJob job1(monidx,-1,-1,0,ix1,1,0);
      jobcount = jobcount + 1;
      if(verbose) cout << jobcount << " : " << job1
                                   << " : layer " << layer << endl;
      jobs[layer].push_back(job1);
      // second operand of real part
      ComplexConvolutionJob job2(monidx,-2,-1,4,ix1,4,0);
      jobcount = jobcount + 1;
      if(verbose) cout << jobcount << " : " << job2
                                   << " : layer " << layer << endl;
      jobs[layer].push_back(job2);
      // first operand of imaginary part
      ComplexConvolutionJob job3(monidx,-1,-1,4,ix1,7,0);
      jobcount = jobcount + 1;
      if(verbose) cout << jobcount << " : " << job3
                                   << " : layer " << layer << endl;
      jobs[layer].push_back(job3);
      // second operand of imaginary part
      ComplexConvolutionJob job4(monidx,-2,-1,0,ix1,10,0);
      jobcount = jobcount + 1;
      if(verbose) cout << jobcount << " : " << job4
                                   << " : layer " << layer << endl;
      jobs[layer].push_back(job4);
   }
   layer = layer + 1;

   for(int i=1; i<nvr; i++)
   {                                                  // f[i] = f[i-1]*x[i]
      ix2 = idx[i];

      freqlaycnt[layer] = freqlaycnt[layer] + 4;
      {
         // first operand of real part
         ComplexConvolutionJob job1(monidx,1,i-1,0,ix2,1,i);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job1
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job1);
         // second operand of real part
         ComplexConvolutionJob job2(monidx,5,i-1,4,ix2,4,i);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job2
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job2);
         // first operand of imaginary part
         ComplexConvolutionJob job3(monidx,1,i-1,4,ix2,7,i);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job3
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job3);
         // second operand of imaginary part
         ComplexConvolutionJob job4(monidx,5,i-1,0,ix2,10,i);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job4
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job4);
      }
      layer = layer + 1;
   }
   // The layer is always incremented after writing,
   // so its value equals the numbers of layers.
   if(layer > laydepth) laydepth = layer;

   if(nvr > 2)
   {
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];           // b[0] = x[n-1]*x[n-2]

      layer = 0;                       // reset layer for backward products

      freqlaycnt[layer] = freqlaycnt[layer] + 4;
      {
         // first operand of real part
         ComplexConvolutionJob job1(monidx,0,ix1,0,ix2,2,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job1
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job1);
         // second operand of real part
         ComplexConvolutionJob job2(monidx,4,ix1,4,ix2,5,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job2
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job2);
         // first operand of imaginary part
         ComplexConvolutionJob job3(monidx,0,ix1,4,ix2,8,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job3
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job3);
         // second operand of imaginary part
         ComplexConvolutionJob job4(monidx,4,ix1,0,ix2,11,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job4
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job4);
      }
      layer = layer + 1;

      for(int i=1; i<nvr-2; i++)
      {                                           // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];

         freqlaycnt[layer] = freqlaycnt[layer] + 4;
         {
            // first operand of real part
            ComplexConvolutionJob job1(monidx,2,i-1,0,ix2,2,i);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job1
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job1);
            // second operand of real part
            ComplexConvolutionJob job2(monidx,6,i-1,4,ix2,5,i);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job2
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job2);
            // first operand of imaginary part
            ComplexConvolutionJob job3(monidx,2,i-1,4,ix2,8,i);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job3
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job3);
            // second operand of imaginary part
            ComplexConvolutionJob job4(monidx,6,i-1,0,ix2,11,i);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job4
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job4);
         }
         layer = layer + 1;
      }
      freqlaycnt[layer] = freqlaycnt[layer] + 4;
      {
         // first operand of real part
         ComplexConvolutionJob job1(monidx,2,nvr-3,-1,-1,2,nvr-2);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job1
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job1);
         // second operand of real part
         ComplexConvolutionJob job2(monidx,6,nvr-3,-2,-1,5,nvr-2);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job2
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job2);
         // first operand of imaginary part
         ComplexConvolutionJob job3(monidx,2,nvr-3,-2,-1,8,nvr-2);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job3
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job3);
         // second operand of imaginary part
         ComplexConvolutionJob job4(monidx,6,nvr-3,-1,-1,11,nvr-2);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job4
                                      << " : layer " << layer << endl;
         jobs[layer].push_back(job4);
      }
      layer = layer + 1;
      // host code uses cross[0] as work space,
      // kernels write to backward[nvr-2]

      if(nvr == 3)
      {                                                 // c[0] = f[0]*x[2]
         ix2 = idx[2];

         layer = 1;
         freqlaycnt[layer] = freqlaycnt[layer] + 4;
         {
            // first operand of real part
            ComplexConvolutionJob job1(monidx,1,0,0,ix2,3,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job1
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job1);
            // second operand of real part
            ComplexConvolutionJob job2(monidx,5,0,4,ix2,6,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job2
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job2);
            // first operand of imaginary part
            ComplexConvolutionJob job3(monidx,1,0,4,ix2,9,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job3
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job3);
            // second operand of imaginary part
            ComplexConvolutionJob job4(monidx,5,0,0,ix2,12,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job4
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job4);
         }
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                          // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;

            // layer = max(i,n-4-i) + 1;
            layer = (i > nvr-4-i) ? i+1 : nvr-4-i+1;
            freqlaycnt[layer] = freqlaycnt[layer] + 4;
            {
               // first operand of real part
               ComplexConvolutionJob job1(monidx,1,i,2,ix2,3,i);
               jobcount = jobcount + 1;
               if(verbose) cout << jobcount << " : " << job1
                                            << " : layer " << layer << endl;
               jobs[layer].push_back(job1);
               // second operand of real part
               ComplexConvolutionJob job2(monidx,5,i,6,ix2,6,i);
               jobcount = jobcount + 1;
               if(verbose) cout << jobcount << " : " << job2
                                            << " : layer " << layer << endl;
               jobs[layer].push_back(job2);
               // first operand of imaginary part
               ComplexConvolutionJob job3(monidx,1,i,6,ix2,9,i);
               jobcount = jobcount + 1;
               if(verbose) cout << jobcount << " : " << job3
                                            << " : layer " << layer << endl;
               jobs[layer].push_back(job3);
               // second operand of imaginary part
               ComplexConvolutionJob job4(monidx,5,i,2,ix2,12,i);
               jobcount = jobcount + 1;
               if(verbose) cout << jobcount << " : " << job4
                                            << " : layer " << layer << endl;
               jobs[layer].push_back(job4);
            }
         }
         ix2 = idx[nvr-1];                        // c[n-3] = f[n-3]*x[n-1]

         layer = nvr-2;
         freqlaycnt[layer] = freqlaycnt[layer] + 4;
         {
            // first operand of real part
            ComplexConvolutionJob job1(monidx,1,nvr-3,0,ix2,3,nvr-3);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job1
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job1);
            // second operand of real part
            ComplexConvolutionJob job2(monidx,5,nvr-3,4,ix2,6,nvr-3);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job2
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job2);
            // first operand of imaginary part
            ComplexConvolutionJob job3(monidx,1,nvr-3,4,ix2,9,nvr-3);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job3
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job3);
            // second operand of imaginary part
            ComplexConvolutionJob job4(monidx,5,nvr-3,0,ix2,12,nvr-3);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job4
                                         << " : layer " << layer << endl;
            jobs[layer].push_back(job4);
         }
      }
   }
}

void ComplexConvolutionJobs::make
 ( int nbr, int *nvr, int **idx, bool verbose )
{
   for(int i=0; i<dimension; i++) freqlaycnt.push_back(0);

   for(int i=0; i<dimension; i++)
   {
      vector<ComplexConvolutionJob> jobvec;
      jobs.push_back(jobvec);
   }
   jobcount = 0;

   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         ix1 = idx[i][0];

         freqlaycnt[0] = freqlaycnt[0] + 4;
         // first operand of real part
         ComplexConvolutionJob job1(i,-1,-1,0,ix1,1,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job1
                                      << " : layer 0" << endl;
         jobs[0].push_back(job1);
         // second operand of real part
         ComplexConvolutionJob job2(i,-2,-1,4,ix1,4,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job2
                                      << " : layer 0" << endl;
         jobs[0].push_back(job2);
         // first operand of imaginary part
         ComplexConvolutionJob job3(i,-1,-1,4,ix1,7,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job3
                                      << " : layer 0" << endl;
         jobs[0].push_back(job3);
         // second operand of imaginary part
         ComplexConvolutionJob job4(i,-2,-1,0,ix1,10,0);
         jobcount = jobcount + 1;
         if(verbose) cout << jobcount << " : " << job4
                                      << " : layer 0" << endl;
         jobs[0].push_back(job4);

         if(laydepth < 1) laydepth = 1; // we have one layer
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         freqlaycnt[0] = freqlaycnt[0] + 4;
         {
            // first operand of real part
            ComplexConvolutionJob job1(i,-1,-1,0,ix1,1,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job1
                                         << " : layer 0" << endl;
            jobs[0].push_back(job1);
            // second operand of real part
            ComplexConvolutionJob job2(i,-2,-1,4,ix1,4,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job2
                                         << " : layer 0" << endl;
            jobs[0].push_back(job2);
            // first operand of imaginary part
            ComplexConvolutionJob job3(i,-1,-1,4,ix1,7,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job3
                                         << " : layer 0" << endl;
            jobs[0].push_back(job3);
            // second operand of imaginary part
            ComplexConvolutionJob job4(i,-2,-1,0,ix1,10,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job4
                                         << " : layer 0" << endl;
            jobs[0].push_back(job4);
         }
         freqlaycnt[0] = freqlaycnt[0] + 4;
         {
            // first operand of real part
            ComplexConvolutionJob job1(i,-1,-1,0,ix2,2,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job1
                                         << " : layer 0" << endl;
            jobs[0].push_back(job1);
            // second operand of real part
            ComplexConvolutionJob job2(i,-2,-1,4,ix2,5,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job2
                                         << " : layer 0" << endl;
            jobs[0].push_back(job2);
            // first operand of imaginary part
            ComplexConvolutionJob job3(i,-1,-1,4,ix2,8,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job3
                                         << " : layer 0" << endl;
            jobs[0].push_back(job3);
            // second operand of real part
            ComplexConvolutionJob job4(i,-2,-1,0,ix2,11,0);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job4
                                         << " : layer 0" << endl;
            jobs[0].push_back(job4);
         }
         freqlaycnt[1] = freqlaycnt[1] + 4;
         {
            // first operand of real part
            ComplexConvolutionJob job1(i,1,0,0,ix2,1,1);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job1
                                         << " : layer 1" << endl;
            jobs[1].push_back(job1);
            // second operand of real part
            ComplexConvolutionJob job2(i,5,0,4,ix2,4,1);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job2
                                         << " : layer 1" << endl;
            jobs[1].push_back(job2);
            // first operand of imaginary part
            ComplexConvolutionJob job3(i,1,0,4,ix2,7,1);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job3
                                         << " : layer 1" << endl;
            jobs[1].push_back(job3);
            // second operand of imaginary part
            ComplexConvolutionJob job4(i,5,0,0,ix2,10,1);
            jobcount = jobcount + 1;
            if(verbose) cout << jobcount << " : " << job4
                                         << " : layer 1" << endl;
            jobs[1].push_back(job4);
         }
         if(laydepth < 2) laydepth = 2; // we have two layers
      }
      else if(nvr[i] > 2)
      {
         make_monomial(nvr[i],idx[i],i,verbose);
      }
   }
}

int ComplexConvolutionJobs::get_dimension ( void ) const
{
   return dimension;
}

int ComplexConvolutionJobs::get_count ( void ) const
{
   return jobcount;
}

int ComplexConvolutionJobs::get_layer_count ( int k ) const
{
   if(k >= dimension)
      return 0;
   else
      return freqlaycnt[k];
}

int ComplexConvolutionJobs::get_depth ( void ) const
{
   return laydepth;
}

ComplexConvolutionJob ComplexConvolutionJobs::get_job ( int k, int i ) const
{
   return jobs[k][i];
}

ComplexConvolutionJobs::~ComplexConvolutionJobs( void )
{
   dimension = 0;
   jobcount = 0;
   laydepth = 0;
}
