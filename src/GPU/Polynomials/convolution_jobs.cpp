// The file convolution_jobs.cpp defines the methods of the class
// ConvolutionJobs, specified in "convolution_jobs.h".

#include <cstdlib>
#include <iostream>
#include "convolution_jobs.h"

using namespace std;

ConvolutionJobs::ConvolutionJobs ( int dim )
{
   dimension = dim;
   jobcount = 0;
   laydepth = 0;
}

void ConvolutionJobs::make_monomial
 ( int nvr, int *idx, int monidx, bool verbose )
{
   int ix1 = idx[0];
   int ix2;
   int layer = 0; // determines the order of execution

   jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

   if(verbose)
   {
      cout << jobcount << " : ";
      cout << "monomial " << monidx << " : ";            // f[0] = cff*x[0]
      cout << "cff * input[" << ix1 << "] to f[0] : ";
      cout << "layer " << layer << endl;
   }
   layer = layer + 1;

   for(int i=1; i<nvr; i++)
   {                                                  // f[i] = f[i-1]*x[i]
      ix2 = idx[i];

      jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

      if(verbose)
      {
         cout << jobcount << " : ";
         cout << "monomial " << monidx << " : ";
         cout << "f[" << i-1 << "] * "
              << "input[" << ix2 << "] to f[" << i << "] : ";
         cout << "layer " << layer << endl;
      }
      layer = layer + 1;
   }
   if(nvr > 2)
   {
      ix1 = idx[nvr-1]; ix2 = idx[nvr-2];           // b[0] = x[n-1]*x[n-2]

      layer = 0;                       // reset layer for backward products

      jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

      if(verbose)
      {
         cout << jobcount << " : ";
         cout << "monomial " << monidx << " : ";
         cout << "input[" << ix1 << "] * "
              << "input[" << ix2 << "] to b[0] : ";
         cout << "layer " << layer << endl;
      }
      layer = layer + 1;

      for(int i=1; i<nvr-2; i++)
      {                                           // b[i] = b[i-1]*x[n-2-i]
         ix2 = idx[nvr-2-i];

         jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

         if(verbose)
         {
            cout << jobcount << " : ";
            cout << "monomial " << monidx << " : ";
            cout << "b[" << i-1 << "] * "
                 << "input[" << ix2 << "] to b[" << i << "] : ";
            cout << "layer " << layer << endl;
         }
      }
      layer = layer + 1;

      jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

      if(verbose)
      {
         cout << jobcount << " : ";
         cout << "monomial " << monidx << " : ";      // b[n-3] = cff*b[n-3]
         cout << "b[" << nvr-3 << "] * cff to "
              << "b[" << nvr-2 << "] : ";
         cout << "layer " << layer << endl;
      }
      layer = layer + 1;
      // host code uses cross[0] as work space,
      // kernels write to backward[nvr-2]

      if(nvr == 3)
      {                                                 // c[0] = f[0]*x[2]
         ix2 = idx[2];

         jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

         if(verbose)
         {
            cout << jobcount << " : ";
            cout << "monomial " << monidx << " : ";
            cout << "f[0] * input[" << ix2 << "] to c[0] : ";
            cout << "layer " << layer << endl;
         }
         layer = layer + 1;
      }
      else
      {
         for(int i=0; i<nvr-3; i++)
         {                                          // c[i] = f[i]*b[n-4-i]
            ix2 = nvr-4-i;

            jobcount = jobcount + 1;
            freqlaycnt[layer] = freqlaycnt[layer] + 1;

            if(verbose)
            {
               cout << jobcount << " : ";
               cout << "monomial " << monidx << " : ";
               cout << "f[" << i << "] * b[" << ix2
                    << "] to c[" << i << "] : ";
               cout << "layer " << layer << endl;
            }
            layer = layer + 1;
         }
         ix2 = idx[nvr-1];                        // c[n-3] = f[n-3]*x[n-1]

         jobcount = jobcount + 1; freqlaycnt[layer] = freqlaycnt[layer] + 1;

         if(verbose)
         {
            cout << jobcount << " : ";
            cout << "monomial " << monidx << " : ";
            cout << "f[" << nvr-3 << "] * input[" << ix2
                 << "] to c[" << nvr-3 << "] : ";
            cout << "layer " << layer << endl;
         }
         layer = layer + 1;
      }
   }
   // The layer is always incremented after writing,
   // so its value equals the numbers of layers.
   if(layer > laydepth) laydepth = layer;
}

void ConvolutionJobs::make ( int nbr, int *nvr, int **idx, bool verbose )
{
   int maxdepth = 2*dimension-2;

   freqlaycnt = new int[maxdepth];
   for(int i=0; i<maxdepth; i++) freqlaycnt[i] = 0;

   int ix1,ix2;

   for(int i=0; i<nbr; i++)
   {
      if(nvr[i] == 1)
      {
         jobcount = jobcount + 1; freqlaycnt[0] = freqlaycnt[0] + 1;

         if(verbose)
         {
            ix1 = idx[i][0];
            cout << jobcount << " : ";
            cout << "monomial " << i << " : ";
            cout << "input[" << ix1 << "] * cff to f[0] : ";
            cout << "layer 0" << endl;
         }
         if(laydepth < 1) laydepth = 1; // we have one layer
      }
      else if(nvr[i] == 2)
      {
         ix1 = idx[i][0]; ix2 = idx[i][1];

         jobcount = jobcount + 1; freqlaycnt[0] = freqlaycnt[0] + 1;

         if(verbose)
         {
            cout << jobcount << " : ";
            cout << "monomial " << i << " : ";
            cout << "cff * input[" << ix1 << "] to c[0] : ";
            cout << "layer 0" << endl;
         }
         jobcount = jobcount + 1; freqlaycnt[0] = freqlaycnt[0] + 1;

         if(verbose)
         {
            cout << jobcount << " : ";
            cout << "monomial " << i << " : ";
            cout << "cff * input[" << ix2 << "] to b[0] : ";
            cout << "layer 0" << endl;
         }
         jobcount = jobcount + 1; freqlaycnt[1] = freqlaycnt[1] + 1;

         if(verbose)
         {
            cout << jobcount << " : ";
            cout << "monomial " << i << " : ";
            cout << "cff * " << "input[" << ix2 << "] to f[0] : ";
            cout << "layer 1" << endl;
         }
         if(laydepth < 2) laydepth = 2; // we have two layers
      }
      else
      {
         make_monomial(nvr[i],idx[i],i,verbose);
      }
   }
}

int ConvolutionJobs::get_dimension ( void ) const
{
   return dimension;
}

int ConvolutionJobs::get_count ( void ) const
{
   return jobcount;
}

int ConvolutionJobs::get_layer_count ( int k ) const
{
   if(k >= 2*dimension - 2)
      return 0;
   else
      return freqlaycnt[k];
}

int ConvolutionJobs::get_depth ( void ) const
{
   return laydepth;
}

ConvolutionJobs::~ConvolutionJobs( void )
{
   free(freqlaycnt);
}
