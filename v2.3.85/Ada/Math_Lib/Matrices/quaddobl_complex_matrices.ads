with QuadDobl_Complex_Ring;              use QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with Generic_Matrices;

package QuadDobl_Complex_Matrices is
  new Generic_Matrices(QuadDobl_Complex_Ring,
                       QuadDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of quad double complex numbers.
