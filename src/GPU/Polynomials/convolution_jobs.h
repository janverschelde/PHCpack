// The file convolution_jobs.h defines the class ConvolutionJobs
// to setup the layers of convolution jobs.

#ifndef __convolution_jobs_h__
#define __convolution_jobs_h__

class ConvolutionJobs
{
   public:

      ConvolutionJobs ( int dim );
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
      /*
       * DESCRIPTION :
       *   Returns the dimension. */

      int get_count ( void ) const;
      /*
       * DESCRIPTION :
       *   Returns the number of convolution jobs. */

      int get_layer_count ( int k ) const;
      /*
       * DESCRIPTION :
       *   Returns the number of jobs in layer k. */

      int get_depth ( void ) const;
      /*
       * DESCRIPTION :
       *   Returns the number of layers of convolution jobs. */

      ~ConvolutionJobs ( void );
      /*
       * DESCRIPTION :
       *   Frees the memory occupied by the frequency table. */
 
   private:

      int dimension;
      int jobcount;
      int laydepth;
      int *freqlaycnt;

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
