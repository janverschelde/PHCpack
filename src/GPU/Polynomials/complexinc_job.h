// The file complexinc_job.h defines the class ComplexIncrementJob
// to store the data of one job to increment the first operand,
// of either the real or imaginary part of complex data.

#ifndef __complexinc_job_h__
#define __complexinc_job_h__

#include <iostream>

class ComplexIncrementJob
{
   public:

      ComplexIncrementJob ( int monomix, int kind, int ix, bool isreal );
      /*
       * DESCRIPTION :
       *   One increment job is uniquely determined by three integers
       *   and one boolean.
       *
       * ON ENTRY :
       *   monomix  index of the monomial;
       *   kind     defines the kind of the increment
       *             1 for forward products, 
       *             2 for backward products,
       *             3 for cross products;
       *   ix       index in forward, backward, or cross products;
       *   isreal   if true, then the real first operand is incremented,
       *            otherwise, the increment happens on imaginary parts. */

      int get_monomial_index ( void ) const;
      // Returns the monomial index of the job.

      int get_increment_kind ( void ) const;
      // Returns the kind of the increment.

      int get_increment_index ( void ) const;
      // Returns the index of the increment.

      bool get_isreal ( void ) const;
      // Returns the isreal value of the job.

      friend std::ostream& operator<<
         ( std::ostream& os, const ComplexIncrementJob& job );
      // Defines the output of a job.

   private:

      int monidx; // monomial index
      int inkind; // kind of increment 
      int incidx; // index of increment
      bool rornot; // real or not?
};

#endif
