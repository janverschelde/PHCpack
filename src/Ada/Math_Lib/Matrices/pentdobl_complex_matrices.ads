with PentDobl_Complex_Ring;              use PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with Generic_Matrices;

package PentDobl_Complex_Matrices is
  new Generic_Matrices(PentDobl_Complex_Ring,PentDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of penta double complex numbers.
