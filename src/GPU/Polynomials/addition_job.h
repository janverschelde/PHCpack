// The file addition_job.h defines the class AdditionJob
// to store the data of one addition job.

#ifndef __addition_job_h__
#define __addition_job_h__

#include <iostream>

class AdditionJob
{
   public:

      AdditionJob ( int atp, int itp,
                    int monidx1, int monidx2, int ix1, int ix2 );
      /*
       * DESCRIPTION :
       *   An addition job is defined by a type and four indices.
       *
       * ON ENTRY :
       *   atp      defines the type of the update, 
       *            1 for forward products,
       *            2 for backward products,
       *            3 for cross products;
       *   itp      defines the type of the increment in the update, 
       *            1 for forward products,
       *            2 for backward products,
       *            3 for cross products;
       *   monix1   index of the update monomial;
       *   monix2   index of the increment monomial,
       *            -1 if the increment is the constant;
       *   ix1      index of the update;
       *   ix2      index of the increment,
       *            -1 if the increment is a coefficient. */

      int get_addition_type ( void ) const;
      // Returns 1, 2, or 3, for f, b, or c respectively.

      int get_increment_type ( void ) const;
      // Returns 1, 2, or 3, for f, b, or c respectively.

      int get_update_monomial ( void ) const;
      // Returns the index of the monomial in the update.

      int get_update_index ( void ) const;
      // Returns the index of the update.

      int get_increment_monomial ( void ) const;
      // Returns the index of the monomial in the increment.

      int get_increment_index ( void ) const;
      // Returns the index of the increment.

      friend std::ostream& operator<<
         ( std::ostream& os, const AdditionJob& job );
      // Defines the output of a job.

   private:

      int adtype; // update type is 1, 2, or 3 for f, b, or c
      int intype; // increment type is 1, 2, or 3 for f, b, or c
      int updmon; // monomial index of the update
      int incmon; // monomial index of the increment
      int updidx; // index of the update
      int incidx; // index of the increment
};

#endif
