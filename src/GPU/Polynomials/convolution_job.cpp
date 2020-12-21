// The file convolution_job.cpp defines the methods of the class
// ConvolutionJob, specified in "convolution_job.h".

#include "convolution_job.h"

ConvolutionJob::ConvolutionJob
 ( int monomialix,
   int input1tp, int input1ix,
   int input2tp, int input2ix,
   int outputtp, int outputix )
{
   monidx = monomialix;  // monomial index
   inp1tp = input1tp;    // type of the first input
   inp1ix = input1ix;    // index of the first input
   inp2tp = input2tp;    // type of the second input
   inp2ix = input2ix;    // index of the second input
   outptp = outputtp;    // type of output
   outidx = outputix;    // output index
}

int ConvolutionJob::get_monomial_index ( void ) const
{
   return monidx;
}

int ConvolutionJob::get_first_type ( void ) const
{
   return inp1tp;
}

int ConvolutionJob::get_first_input ( void ) const
{
   return inp1ix;
}

int ConvolutionJob::get_second_type ( void ) const
{
   return inp2tp;
}

int ConvolutionJob::get_second_input ( void ) const
{
   return inp2ix;
}

int ConvolutionJob::get_output_type ( void ) const
{
   return outptp;
}

int ConvolutionJob::get_output_index ( void ) const
{
   return outidx;
}

std::ostream& operator<< ( std::ostream& os, const ConvolutionJob& job )
{
   os << "monomial " << job.monidx << " : ";
   if(job.inp1tp == -1) os << "cff * ";
   if(job.inp1tp == 0) os << "input[" << job.inp1ix << "] * ";
   if(job.inp1tp == 1) os << "f[" << job.inp1ix << "] * ";
   if(job.inp1tp == 2) os << "b[" << job.inp1ix << "] * ";
   if(job.inp1tp == 3) os << "c[" << job.inp1ix << "] * ";
   if(job.inp2tp == -1) os << "cff to ";
   if(job.inp2tp == 0) os << "input[" << job.inp2ix << "] to ";
   if(job.inp2tp == 1) os << "f[" << job.inp2ix << "] to ";
   if(job.inp2tp == 2) os << "b[" << job.inp2ix << "] to ";
   if(job.inp2tp == 3) os << "c[" << job.inp2ix << "] to ";
   if(job.outptp == 1) os << "f[" << job.outidx << "]";
   if(job.outptp == 2) os << "b[" << job.outidx << "]";
   if(job.outptp == 3) os << "c[" << job.outidx << "]";

   return os;
}
