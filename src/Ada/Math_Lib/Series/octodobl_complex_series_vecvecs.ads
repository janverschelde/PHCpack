with Generic_VecVecs;
with OctoDobl_Complex_Series_Ring;
with OctoDobl_Complex_Series_Vectors;

package OctoDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(OctoDobl_Complex_Series_Ring,
                      OctoDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with octo double complex numbers as coefficients.
