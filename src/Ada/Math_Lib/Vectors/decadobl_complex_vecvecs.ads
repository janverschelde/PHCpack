with Generic_VecVecs;
with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;

package DecaDobl_Complex_VecVecs is 
  new Generic_VecVecs(DecaDobl_Complex_Ring,DecaDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of deca double complex numbers.
