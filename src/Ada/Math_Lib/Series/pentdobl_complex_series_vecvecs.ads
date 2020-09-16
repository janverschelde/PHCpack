with Generic_VecVecs;
with PentDobl_Complex_Series_Ring;
with PentDobl_Complex_Series_Vectors;

package PentDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(PentDobl_Complex_Series_Ring,
                      PentDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with penta double complex numbers as coefficients.
