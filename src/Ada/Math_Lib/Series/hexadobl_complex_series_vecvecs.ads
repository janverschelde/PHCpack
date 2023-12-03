with Generic_VecVecs;
with HexaDobl_Complex_Series_Ring;
with HexaDobl_Complex_Series_Vectors;

package HexaDobl_Complex_Series_VecVecs is 
  new Generic_VecVecs(HexaDobl_Complex_Series_Ring,
                      HexaDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines vectors of vectors over the ring of truncated power series
--   with hexa double complex numbers as coefficients.
