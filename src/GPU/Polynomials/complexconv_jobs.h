// The file complexconv_jobs.h defines the class ComplexConvolutionJobs
// to setup the layers of convolution jobs on complex data.

#ifndef __complexconv_jobs_h__
#define __complexconv_jobs_h__

#include <vector>
#include "complexconv_job.h"

using namespace std;

class ComplexConvolutionJobs
{
   public:

      ComplexConvolutionJobs ( int dim );
      /*
       * DESCRIPTION :
       *   Sets the dimension to dim.
       *   The dimension is the number of variables in a polynomial. */

      void make ( int nbr, int *nvr, int **idx, bool verbose );
      /*
       * DESCRIPTION :
       *   Makes the frequency table of convolution jobs,
       *   given the supports of a polynomial in several variables.
       *
       * REQUIRED : nbr <= get_dimension().
       *
       * ON ENTRY :
       *   nbr     number of monomials in the polynomial;
       *   nvr     array of nbr counts of the variables in monomials;
       *   idx     array of nbr indices to the participating variables;
       *   verbose if true, then one line is written for each job,
       *           if false, then the constructor remains silent. */

      int get_dimension ( void ) const;
      // Returns the dimension.

      int get_count ( void ) const;
      // Returns the number of convolution jobs.

      int get_layer_count ( int k ) const;
      // Returns the number of jobs in layer k.

      int get_depth ( void ) const;
      // Returns the number of layers of convolution jobs.

      ComplexConvolutionJob get_job ( int k, int i ) const;
      /*
       * DESCRIPTION :
       *   Returns the i-th convolution job at layer k.
       *
       * REQUIRED : k < get_depth() and i < get_layer_count(k). */

      ~ComplexConvolutionJobs ( void );
      /*
       * DESCRIPTION :
       *   Frees the memory occupied by the frequency table. */
 
   private:

      int dimension; // the dimension is the total number of variables
      int jobcount;  // total number of convolution jobs
      int laydepth;  // number of layers is bounded by the dimension

      vector<int> freqlaycnt; // frequency table of jobs at each layer

      vector< vector<ComplexConvolutionJob> > jobs;

      void make_monomial ( int nvr, int *idx, int monidx, bool verbose );
      /*
       * DESCRIPTION :
       *   Updates the frequency table for one monomial.
       *
       * ON ENTRY :
       *   nvr     number of variables in the monomial;
       *   idx     array of nvr indices to the participating variables;
       *   monidx  index of the monomial;
       *   verbose if true, then one line is written for each job,
       *           if false, then the constructor remains silent. */
};

#endif
