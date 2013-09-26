with DoblDobl_Complex_Ring;              use DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with Generic_Matrices;

package DoblDobl_Complex_Matrices is
  new Generic_Matrices(DoblDobl_Complex_Ring,
                       DoblDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of double double complex numbers.
