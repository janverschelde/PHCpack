with Generic_VecVecs;
with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;

package DoblDobl_Complex_VecVecs is 
  new Generic_VecVecs(DoblDobl_Complex_Ring,DoblDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of double double complex numbers.
