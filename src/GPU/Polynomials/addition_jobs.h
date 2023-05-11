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

      AdditionJobs ( int dim, int nbr );
      /*
       * DESCRIPTION :
       *   Sets the dimensions to dim and nbr.
       *
       * ON ENTRY :
       *   dim      the total number of variables;
       *   nbr      the number of monomials. */

      void make ( int nbr, int *nvr, int **idx, bool verbose );
      /*
       * DESCRIPTION :
       *   Makes the reduction tree of addition jobs,
       *   given the supports of a polynomial in several variables.
       *
       * REQUIRED : nbr <= get_dimension().
       *
       * ON ENTRY :
       *   nbr      number of monomials in the polynomial;
       *   nvr      array of nbr counts of the variables in monomials;
       *   idx      array of nbr support vectors,
       *            idx[k] has nvr[k] integers, idx[k][i] is the index
       *            of the i-th variable in monomial k;
       *   verbose  if true, then one line is written for each job,
       *            if false, then the constructor remains silent. */

      int get_number_of_variables ( void ) const;
      // Returns the number of variables.

      int get_number_of_monomials ( void ) const;
      // Returns the number of monomials.

      int get_count ( void ) const;
      // Returns the number of addition jobs.

      int get_layer_count ( int k ) const;
      // Returns the number of jobs in layer k.

      int get_depth ( void ) const;
      // Returns the number of layers of addition jobs.

      int get_differential_count ( int k ) const;
      // Returns the number of monomials that contain variable k.

      int get_differential_index ( int k, int i ) const;
      /*
       * DESCRIPTION :
       *   Returns the index of the i-th monomial that contains
       *   the variable k.
       *
       * REQUIRED : i <= get_differential_count(k). */

      int position ( int n, int *idx, int k );
      /*
       * DESCRIPTION :
       *   Returns the position of k in the array idx, or
       *   return -1 if k does not occur in the values of idx.
       *   Needed to find which derivative a cross product stores.
       *
       * ON ENTRY :
       *   n        number of entries in idx;
       *   idx      sequence of n indices;
       *   k        one index. */

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

      int nbrvar;
      int nbrmon;
      int jobcount;
      int laydepth;

      vector<int> freqlaycnt; // frequency table of jobs at each layer

      vector<int> difcnt; // counts of the differential indices

      vector< vector<int> > difidx;

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

      void recursive_first_make
       ( int level, int stride, int nbr, int *nvr, bool verbose );
      /*
       * DESCRIPTION :
       *   Adds the jobs recursively for the first derivative,
       *   starting at the top.
       *
       * ON ENTRY :
       *   level    current layer of jobs;
       *   stride   current stride;
       *   nbr      value of difcnt[0];
       *   nvr      nbr integers count the variables in each monomial,
       *            nvr[k] is the number of variables in monomial k;
       *   verbose  if true, writes one line per job added. */

      void recursive_other_make
       ( int level, int stride, int nbr, int *nvr, int **idx, int varidx,
         bool verbose );
      /*
       * DESCRIPTION :
       *   Adds the jobs recursively for the derivative, other than
       *   the first one, starting at the top.
       *
       * ON ENTRY :
       *   level    current layer of jobs;
       *   stride   current stride;
       *   nbr      value of difcnt[varidx];
       *   nvr      nbr integers count the variables in each monomial,
       *            nvr[k] is the number of variables in monomial k;
       *   idx      array of nbr support vectors,
       *            idx[k] has nvr[k] integers, idx[k][i] is the index
       *            of the i-th variable in monomial k;
       *   vardidx  index of the variable for the derivative;
       *   verbose  if true, writes one line per job added. */

      void differential_index_count
       ( int dim, int nbr, int *nvr, int **idx,
         vector<int> &cnt, bool verbose );
      /*
       * DESCRIPTION :
       *   Counts the number of monomials where the variables appear.
       *
       * ON ENTRY :
       *   dim      the total number of variables;
       *   nbr      number of monomials;
       *   nvr      nbr integers, nvr[k] counts the number of variables
       *            in the k-th monomial;
       *   idx      array of nbr support vectors,
       *            idx[k] has nvr[k] integers, idx[k][i] is the index
       *            of the i-th variable in monomial k;
       *   cnt      space for dim integers;
       *   verbose  indicates if output will be written during the count.
       *
       * ON RETURN :
       *   cnt      cnt[k] counts the number of monomials that contain k. */ 

      void make_differential_indices
       ( int dim, int nbr, int *nvr, int **idx,
         vector<int> &cnt, vector< vector<int> > &difidx, bool verbose );
      /*
       * DESCRIPTION :
       *   Defines the rearrangment of the monomial indices according
       *   to which monomials the variables appear.
       *
       * ON ENTRY :
       *   dim      the total number of variables;
       *   nbr      number of monomials;
       *   nvr      nbr integers, nvr[k] counts the number of variables
       *            in the k-th monomial;
       *   idx      array of nbr support vectors,
       *            idx[k] has nvr[k] integers, idx[k][i] is the index
       *            of the i-th variable in monomial k;
       *   cnt      dim counts, cnt[k] counts the number of monomials
       *            that contain k;
       *   difidx   array of dim integer vectors,
       *            difidx[k] has space for cnt[k]+1 numbers;
       *   verbose  indicates if output needs to be written;
       *
       * ON RETURN :
       *   cnt      cnt[k] is reduced by one if there is a monomial
       *            where only the variable k occurs;
       *   difidx   differential indices,
       *            difidx[k] lists all indices to those monomials
       *            that contain variable k, with difidx[k][0]
       *            == -1, if monomial coefficients play no role,
       *            ==  i, if cff[i] contributes to derivative k,
       *            difidx[k][i] is the index of the i-th monomial
       *            that contains the variable k. */
};

#endif
