with Generic_VecVecs;
with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;

package PentDobl_Complex_VecVecs is 
  new Generic_VecVecs(PentDobl_Complex_Ring,PentDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of penta double complex numbers.
