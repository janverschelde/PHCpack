with DecaDobl_Complex_Ring;              use DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with Generic_Matrices;

package DecaDobl_Complex_Matrices is
  new Generic_Matrices(DecaDobl_Complex_Ring,DecaDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of deca double complex numbers.
