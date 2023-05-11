// The file complexconv_job.h defines the class ComplexConvolutionJob
// to store the data of one convolution job on complex data.

#ifndef __complexconv_job_h__
#define __complexconv_job_h__

#include <iostream>

class ComplexConvolutionJob
{
   public:

      ComplexConvolutionJob ( int monomialix,
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
       *                  -2 for the imaginary part of a monomial coefficient,
       *                  -1 for the real part of a monomial coefficient,
       *                   0 for the real parts of an input series,
       *                   1 for forward, 2 for backward, 3 for cross,
       *                   4 for the imaginary parts of an input series,
       *                   5 for forward, 6 for backward, 7 for cross;
       *   input1ix     index of the values of the first input series, 
       *                if input1tp equals 0 (or 4),
       *                or the index in forward, backward or cross,
       *                if input1tp is 1, 2, or 3 respectively,
       *                or 5, 6, 7 for the imaginary parts;
       *   input2tp     type of the second input:
       *                  -2 for the imaginary part of a monomial coefficient,
       *                  -1 for the real part of a monomial coefficient,
       *                   0 for the real parts of an input series,
       *                   1 for forward, 2 for backward, 3 for cross;
       *                   4 for the imaginary parts of an input series,
       *                   5 for forward, 6 for backward, 7 for cross;
       *   input2ix     index of the values of the second input series,
       *                if input2tp equals 0 (or 4),
       *                or the index in forward, backward or cross,
       *                if input2tp is 1, 2, or 3 respectively,
       *                or 5, 6, 7 for the imaginary parts;
       *   outputtp     type of output, either forward, backward,
       *                or cross product, respectively 1, 2, or 3,
       *                for first operand of the real parts of the output,
       *                4, 5, or 6 for the second operand of real parts,
       *                7, 8, or 9 for the first operand of imaginary parts,
       *                10, 11, or 12 for the second operand
       *                of imaginary parts;
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
         ( std::ostream& os, const ComplexConvolutionJob& job );
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
