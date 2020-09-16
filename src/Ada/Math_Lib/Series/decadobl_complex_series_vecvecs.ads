with Generic_VecVecs;
with DecaDobl_Complex_Series_Ring;
with DecaDobl_Complex_Series_Vectors;

package DecaDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(DecaDobl_Complex_Series_Ring,
                      DecaDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with deca double complex numbers as coefficients.
