// The file addition_jobs.h defines the class AdditionJobs
// to setup the layers of addition jobs.

#ifndef __addition_jobs_h__
#define __addition_jobs_h__

#include <vector>
#include "addition_job.h"

using namespace std;

class AdditionJobs
{
   public:

      AdditionJobs ( int nbr );
      /*
       * DESCRIPTION :
       *   Sets the dimension to nbr.
       *   The dimension is the number of monomials. */

      void make ( int nbr, int *nvr, bool verbose );
      /*
       * DESCRIPTION :
       *   Makes the reduction tree of addition jobs,
       *   given the supports of a polynomial in several variables.
       *
       * REQUIRED : nbr <= get_dimension().
       *
       * ON ENTRY :
       *   nbr     number of monomials in the polynomial;
       *   nvr     array of nbr counts of the variables in monomials;
       *   verbose if true, then one line is written for each job,
       *           if false, then the constructor remains silent. */

      int get_dimension ( void ) const;
      // Returns the dimension, the number of terms.

      int get_count ( void ) const;
      // Returns the number of addition jobs.

      int get_layer_count ( int k ) const;
      // Returns the number of jobs in layer k.

      int get_depth ( void ) const;
      // Returns the number of layers of addition jobs.

      AdditionJob get_job ( int k, int i ) const;
      /*
       * DESCRIPTION :
       *   Returns the i-th addition job at layer k.
       *
       * REQUIRED : k < get_depth() and i < get_layer_count(k). */

      ~AdditionJobs ( void );
      /*
       * DESCRIPTION :
       *   Frees the memory occupied by the frequency table. */
 
   private:

      int dimension;
      int jobcount;
      int laydepth;
      int *freqlaycnt;

      vector< vector<AdditionJob> > jobs;

      void recursive_start ( int nbr, int *level, int *stride );
      /*
       * DESCRIPTION :
       *   Determines the value for the level and stride,
       *   given the number of monomials in nbr,
       *   to start the recursive_make. */

      void recursive_make
       ( int level, int stride, int nbr, int *nvr, bool verbose );
      /*
       * DESCRIPTION :
       *   Adds the jobs recursively, starting at the top.
       *
       * ON ENTRY :
       *   level    current layer of jobs;
       *   stride   current stride;
       *   nbr      number of monomials, excluding the constant;
       *   nvr      nbr integers count the variables in each monomial,
       *            nvr[k] is the number of variables in monomial k;
       *   verbose  if true, writes one line per job added. */
};

#endif
