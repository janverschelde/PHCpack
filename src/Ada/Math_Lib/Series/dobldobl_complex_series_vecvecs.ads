with Generic_VecVecs;
with DoblDobl_Complex_Series_Ring;
with DoblDobl_Complex_Series_Vectors;

package DoblDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(DoblDobl_Complex_Series_Ring,
                      DoblDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with double double complex numbers as coefficients.
