with Generic_VecVecs;
with QuadDobl_Complex_Series_Ring;
with QuadDobl_Complex_Series_Vectors;

package QuadDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(QuadDobl_Complex_Series_Ring,
                      QuadDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with quad double complex numbers as coefficients.
