with Generic_VecVecs;
with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;

package OctoDobl_Complex_VecVecs is 
  new Generic_VecVecs(OctoDobl_Complex_Ring,OctoDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of octo double complex numbers.
