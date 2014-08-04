with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_VecVecs;
with Generic_Lists_of_Vectors;

package QuadDobl_Complex_VecLists is 
  new Generic_Lists_of_Vectors(QuadDobl_Complex_Ring,
                               QuadDobl_Complex_Vectors,
                               QuadDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines lists of links to vectors of quad double complex numbers.
