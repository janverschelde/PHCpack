with Generic_VecVecs;
with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;

package QuadDobl_Complex_VecVecs is 
  new Generic_VecVecs(QuadDobl_Complex_Ring,QuadDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of quad double complex numbers.
