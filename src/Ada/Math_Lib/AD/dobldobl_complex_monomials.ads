with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_VecVecs;
with Generic_Monomials;

package DoblDobl_Complex_Monomials is
  new Generic_Monomials(DoblDobl_Complex_Ring,
                        DoblDobl_Complex_Vectors,
                        DoblDobl_Complex_VecVecs);

-- DESCRIPTION :
--   Defines monomials over the ring of double double complex numbers.
