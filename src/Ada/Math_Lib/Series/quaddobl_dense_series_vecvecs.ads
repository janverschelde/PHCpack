with Generic_VecVecs;
with QuadDobl_Dense_Series_Ring;
with QuadDobl_Dense_Series_Vectors;

package QuadDobl_Dense_Series_VecVecs is 
  new Generic_VecVecs(QuadDobl_Dense_Series_Ring,
                      QuadDobl_Dense_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with quad double complex numbers as coefficients.
