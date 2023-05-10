// The file convolution_job.h defines the class ConvolutionJob
// to store the data of one convolution job.

#ifndef __convolution_job_h__
#define __convolution_job_h__

#include <iostream>

class ConvolutionJob
{
   public:

      ConvolutionJob ( int monomialix,
                       int input1tp, int input1ix,
                       int input2tp, int input2ix,
                       int outputtp, int outputix );
      /*
       * DESCRIPTION :
       *   One convolution job is uniquely determined by five integers.
       *
       * ON ENTRY :
       *   monomialix   index of the monomial;
       *   input1tp     type of the first input:
       *                  -1 for monomial coefficient,
       *                   0 for an input series,
       *                   1 for forward, 2 for backward, 3 for cross;
       *   input1ix     index of the values of the first input series, 
       *                if input1tp equals 0,
       *                or the index in forward, backward or cross,
       *                if input1tp is 1, 2, or 3 respectively;
       *   input2tp     type of the second input:
       *                  -1 for monomial coefficient,
       *                   0 for an input series,
       *                   1 for forward, 2 for backward, 3 for cross;
       *   input2ix     index of the values of the second input series,
       *                if input2tp equals 0,
       *                or the index in forward, backward or cross,
       *                if input2tp is 1, 2, or 3 respectively;
       *   outputtp     type of output, either forward, backward,
       *                or cross product, respectively 1, 2, or 3;
       *   outputix     index in the forward, backward, or cross product. */

      int get_monomial_index ( void ) const;
      // Returns the monomial index of the job.

      int get_first_type ( void ) const;
      // Returns the type of the first input of the job.

      int get_first_input ( void ) const;
      // Returns the first input index of the job.

      int get_second_type ( void ) const;
      // Returns the type of the second input of the job.

      int get_second_input ( void ) const;
      // Returns the second input index of the job.

      int get_output_type ( void ) const;
      // Returns the output type of the job.

      int get_output_index ( void ) const;
      // Returns the output index of the job.

      friend std::ostream& operator<<
         ( std::ostream& os, const ConvolutionJob& job );
      // Defines the output of a job.

   private:

      int monidx; // monomial index
      int inp1tp; // first input type
      int inp1ix; // index of the first input
      int inp2tp; // second input type
      int inp2ix; // index of the second input
      int outptp; // type of output
      int outidx; // output index
};

#endif
