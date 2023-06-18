// The file write_job_counts.cpp defines the functions
// specified in write_job_counts.h.

#include <iostream>
#include "write_job_counts.h"

using namespace std;

void write_convolution_counts ( ConvolutionJobs jobs )
{
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
}

void write_complexconv_counts ( ComplexConvolutionJobs jobs )
{
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
}

void write_addition_counts ( AdditionJobs jobs )
{
   cout << "number of addition jobs : " << jobs.get_count() << endl;
   cout << "number of layers : " << jobs.get_depth() << endl;
   cout << "frequency of layer counts :" << endl;

   int checksum = 0;

   for(int i=0; i<jobs.get_depth(); i++)
   {
      cout << i << " : " << jobs.get_layer_count(i) << endl;
      checksum = checksum + jobs.get_layer_count(i); 
   }
   cout << "layer count sum : " << checksum << endl;
}

void write_complexadd_counts ( ComplexAdditionJobs jobs )
{
   cout << "number of addition jobs : " << jobs.get_count() << endl;
   cout << "number of layers : " << jobs.get_depth() << endl;
   cout << "frequency of layer counts :" << endl;

   int checksum = 0;

   for(int i=0; i<jobs.get_depth(); i++)
   {
      cout << i << " : " << jobs.get_layer_count(i) << endl;
      checksum = checksum + jobs.get_layer_count(i); 
   }
   cout << "layer count sum : " << checksum << endl;
}

void convolution_operation_counts
 ( int deg, ConvolutionJobs cnvjobs,
   long long int *addcnt, long long int *mulcnt, int vrblvl )
{
   const long long int nbrcnv = cnvjobs.get_count();
   const long long int degp1 = deg+1;
   const long long int cnvaddcnt = degp1*deg*nbrcnv;
   const long long int cnvmulcnt = degp1*degp1*nbrcnv;

   if(vrblvl > 0)
   {
      const long long int totalcnt = cnvaddcnt + cnvmulcnt;

      cout << "truncation degree : " << deg << endl;
      cout << "number of operations in " << nbrcnv
           << " convolution jobs : " << endl;
      cout << "  #additions : " << cnvaddcnt;
      cout << ", #multiplications : " << cnvmulcnt << endl;
      cout << "total #operations : " << totalcnt << endl;
   }
   *addcnt = cnvaddcnt;
   *mulcnt = cnvmulcnt;
}

void complexconv_operation_counts
 ( int deg, ComplexConvolutionJobs cnvjobs,
   long long int *addcnt, long long int *mulcnt, int vrblvl )
{
   const long long int nbrcnv = cnvjobs.get_count();
   const long long int degp1 = deg+1;
   const long long int cnvaddcnt = degp1*deg*nbrcnv;
   const long long int cnvmulcnt = degp1*degp1*nbrcnv;

   if(vrblvl > 0)
   {
      const long long int totalcnt = cnvaddcnt + cnvmulcnt;

      cout << "truncation degree : " << deg << endl;
      cout << "number of operations in " << nbrcnv
           << " convolution jobs : " << endl;
      cout << "  #additions : " << cnvaddcnt;
      cout << ", #multiplications : " << cnvmulcnt << endl;
      cout << "total #operations : " << totalcnt << endl;
   }
   *addcnt = cnvaddcnt;
   *mulcnt = cnvmulcnt;
}

void complexinc_operation_counts
 ( int deg, ComplexIncrementJobs incjobs,
   long long int *addcnt, int vrblvl )
{
   const long long int nbrinc = incjobs.get_count();
   const long long int degp1 = deg+1;
   const long long int incaddcnt = degp1*nbrinc;

   if(vrblvl > 0)
   {
      cout << "truncation degree : " << deg << endl;
      cout << "number of additions in " << nbrinc
           << " increment jobs : " << incaddcnt << endl;
   }
   *addcnt = incaddcnt;
}

void write_operation_counts
 ( int deg, ConvolutionJobs cnvjobs, AdditionJobs addjobs )
{
   const int nbrcnv = cnvjobs.get_count();
   const int nbradd = addjobs.get_count();

   const long int cnvaddcnt = (deg+1)*deg*nbrcnv;
   const long int cnvmulcnt = (deg+1)*(deg+1)*nbrcnv;
   const long int updaddcnt = (deg+1)*nbradd;
   const long int totalcnt = cnvaddcnt + cnvmulcnt + updaddcnt;

   cout << "truncation degree : " << deg << endl;
   cout << "number of operations in " << nbrcnv
        << " convolution jobs : " << endl;
   cout << "  #additions : " << cnvaddcnt;
   cout << ", #multiplications : " << cnvmulcnt << endl;

   cout << "number of operations in " << nbradd
        << " addition jobs : " << endl;
   cout << "  #additions : " << updaddcnt << endl;

   cout << "total #additions : " << cnvaddcnt+updaddcnt << endl;
   cout << "total #multiplications : " << cnvmulcnt << endl;
   cout << "total #operations : " << totalcnt << endl;
}

void write_complexop_counts
 ( int deg, ComplexConvolutionJobs cnvjobs, ComplexAdditionJobs addjobs )
{
   const int nbrcnv = cnvjobs.get_count();
   const int nbradd = addjobs.get_count();

   const long int cnvaddcnt = (deg+1)*deg*nbrcnv;
   const long int cnvmulcnt = (deg+1)*(deg+1)*nbrcnv;
   const long int updaddcnt = (deg+1)*nbradd;
   const long int totalcnt = cnvaddcnt + cnvmulcnt + updaddcnt;

   cout << "truncation degree : " << deg << endl;
   cout << "number of operations in " << nbrcnv
        << " convolution jobs : " << endl;
   cout << "  #additions : " << cnvaddcnt;
   cout << ", #multiplications : " << cnvmulcnt << endl;

   cout << "number of operations in " << nbradd
        << " addition jobs : " << endl;
   cout << "  #additions : " << updaddcnt << endl;

   cout << "total #additions : " << cnvaddcnt+updaddcnt << endl;
   cout << "total #multiplications : " << cnvmulcnt << endl;
   cout << "total #operations : " << totalcnt << endl;
}
