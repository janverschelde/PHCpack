// The file complexinc_jobs.h defines the class ComplexIncrementJobs
// to setup the layers of increment jobs on complex data.

#ifndef __complexinc_jobs_h__
#define __complexinc_jobs_h__

#include <vector>
#include "complexconv_job.h"
#include "complexconv_jobs.h"
#include "complexinc_job.h"

using namespace std;

class ComplexIncrementJobs
{
   public:

      ComplexIncrementJobs ( ComplexConvolutionJobs cnvjobs, bool verbose );
      /*
       * DESCRIPTION :
       *   Makes the table of all increment jobs corresponding
       *   to the convolution jobs in cnvjobs.
       *   If verbose, then one line is written for each job. */

      int get_dimension ( void ) const;
      // Returns the dimension.

      int get_count ( void ) const;
      // Returns the number of jobs.

      int get_layer_count ( int k ) const;
      // Returns the number of jobs in layer k.

      int get_depth ( void ) const;
      // Returns the number of layers of jobs.

      ComplexIncrementJob get_job ( int k, int i ) const;
      /*
       * DESCRIPTION :
       *   Returns the i-th job at layer k.
       *
       * REQUIRED : k < get_depth() and i < get_layer_count(k). */

      ~ComplexIncrementJobs ( void );
      /*
       * DESCRIPTION :
       *   Frees the memory occupied by the frequency table. */
 
   private:

      int dimension; // the dimension is the total number of variables
      int jobcount;  // total number of jobs
      int laydepth;  // number of layers is bounded by the dimension

      vector<int> freqlaycnt; // frequency table of jobs at each layer

      vector< vector<ComplexIncrementJob> > jobs;

      ComplexIncrementJob make_job ( ComplexConvolutionJob cnvjob ) const;
      /*
       * DESCRIPTION :
       *   Returns the increment job corresponding to the convolution job.
       *
       * REQUIRED :
       *   cnvjob.outputtp is in { 1, 2, 3, 7, 8, 9 }.
       *
       * ON ENTRY :
       *   cnvjob   convolution job. */
};

#endif
