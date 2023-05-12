// The file complexadd_job.h defines the class ComplexAdditionJob
// to store the data of one addition job on complex data.

#ifndef __complexadd_job_h__
#define __complexadd_job_h__

#include <iostream>

class ComplexAdditionJob
{
   public:

      ComplexAdditionJob ( int atp, int itp,
                           int monidx1, int monidx2, int ix1, int ix2 );
      /*
       * DESCRIPTION :
       *   An addition job is defined by a type and four indices.
       *
       * ON ENTRY :
       *   atp      defines the type of the update, 
       *            1 for forward products, real part,
       *            2 for backward products, real part,
       *            3 for cross products, real part,
       *            4 for forward products, imaginary part,
       *            5 for backward products, imaginary part,
       *            6 for cross products, imaginary part;
       *   itp      defines the type of the increment in the update, 
       *             1 for forward products, first real operand,
       *             2 for backward products, first real operand,
       *             3 for cross products, first real operand;
       *             4 for forward products, second real operand,
       *             5 for backward products, second real operand,
       *             6 for cross products, second real operand,
       *             7 for forward products, first imaginary operand,
       *             8 for backward products, first imaginary operand,
       *             9 for cross products, first imaginary operand;
       *            10 for forward products, second imaginary operand,
       *            11 for backward products, second imaginary operand,
       *            12 for cross products, second imaginary operand;
       *   monix1   index of the update monomial;
       *   monix2   index of the increment monomial,
       *            -1 if the increment is the constant, real part,
       *            -2 if the increment is the constant, imaginary part;
       *   ix1      index of the update;
       *   ix2      index of the increment,
       *            -1 if the increment is a coefficient, real part,
       *            -2 if the increment is a coefficient, imaginary part. */

      int get_addition_type ( void ) const;
      // Returns the type of the update.

      int get_increment_type ( void ) const;
      // Returns the type of the increment in the update.

      int get_update_monomial ( void ) const;
      // Returns the index of the monomial in the update.

      int get_update_index ( void ) const;
      // Returns the index of the update.

      int get_increment_monomial ( void ) const;
      // Returns the index of the monomial in the increment.

      int get_increment_index ( void ) const;
      // Returns the index of the increment.

      friend std::ostream& operator<<
         ( std::ostream& os, const ComplexAdditionJob& job );
      // Defines the output of a job.

   private:

      int adtype; // update type is in 1, 2, .., 6 
      int intype; // increment type is in 1, 2, .., 12
      int updmon; // monomial index of the update
      int incmon; // monomial index of the increment
      int updidx; // index of the update
      int incidx; // index of the increment
};

#endif
