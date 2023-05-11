// The file complexconv_job.cpp defines the methods of
// the class ComplexConvolutionJob, specified in "complexconv_job.h".

#include "complexconv_job.h"

ComplexConvolutionJob::ComplexConvolutionJob
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

int ComplexConvolutionJob::get_monomial_index ( void ) const
{
   return monidx;
}

int ComplexConvolutionJob::get_first_type ( void ) const
{
   return inp1tp;
}

int ComplexConvolutionJob::get_first_input ( void ) const
{
   return inp1ix;
}

int ComplexConvolutionJob::get_second_type ( void ) const
{
   return inp2tp;
}

int ComplexConvolutionJob::get_second_input ( void ) const
{
   return inp2ix;
}

int ComplexConvolutionJob::get_output_type ( void ) const
{
   return outptp;
}

int ComplexConvolutionJob::get_output_index ( void ) const
{
   return outidx;
}

std::ostream& operator<<
 ( std::ostream& os, const ComplexConvolutionJob& job )
{
   os << "monomial " << job.monidx << " : ";
   if(job.inp1tp == -2) os << "cffim * ";
   if(job.inp1tp == -1) os << "cffre * ";
   if(job.inp1tp ==  0) os << "input[" << job.inp1ix << "]re * ";
   if(job.inp1tp ==  1) os << "f[" << job.inp1ix << "]re * ";
   if(job.inp1tp ==  2) os << "b[" << job.inp1ix << "]re * ";
   if(job.inp1tp ==  3) os << "c[" << job.inp1ix << "]re * ";
   if(job.inp1tp ==  4) os << "input[" << job.inp1ix << "]im * ";
   if(job.inp1tp ==  5) os << "f[" << job.inp1ix << "]im * ";
   if(job.inp1tp ==  6) os << "b[" << job.inp1ix << "]im * ";
   if(job.inp1tp ==  7) os << "c[" << job.inp1ix << "]im * ";
   if(job.inp2tp == -2) os << "cffim to ";
   if(job.inp2tp == -1) os << "cffre to ";
   if(job.inp2tp ==  0) os << "input[" << job.inp2ix << "]re to ";
   if(job.inp2tp ==  1) os << "f[" << job.inp2ix << "]re to ";
   if(job.inp2tp ==  2) os << "b[" << job.inp2ix << "]re to ";
   if(job.inp2tp ==  3) os << "c[" << job.inp2ix << "]re to ";
   if(job.inp2tp ==  4) os << "input[" << job.inp2ix << "]im to ";
   if(job.inp2tp ==  5) os << "f[" << job.inp2ix << "]im to ";
   if(job.inp2tp ==  6) os << "b[" << job.inp2ix << "]im to ";
   if(job.inp2tp ==  7) os << "c[" << job.inp2ix << "]im to ";
   if(job.outptp ==  1) os << "f[" << job.outidx << "]re^a";
   if(job.outptp ==  2) os << "b[" << job.outidx << "]re^a";
   if(job.outptp ==  3) os << "c[" << job.outidx << "]re^a";
   if(job.outptp ==  4) os << "f[" << job.outidx << "]re^b";
   if(job.outptp ==  5) os << "b[" << job.outidx << "]re^b";
   if(job.outptp ==  6) os << "c[" << job.outidx << "]re^b";
   if(job.outptp ==  7) os << "f[" << job.outidx << "]im^a";
   if(job.outptp ==  8) os << "b[" << job.outidx << "]im^a";
   if(job.outptp ==  9) os << "c[" << job.outidx << "]im^a";
   if(job.outptp == 10) os << "f[" << job.outidx << "]im^b";
   if(job.outptp == 11) os << "b[" << job.outidx << "]im^b";
   if(job.outptp == 12) os << "c[" << job.outidx << "]im^b";

   return os;
}
