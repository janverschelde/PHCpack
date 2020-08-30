with OctoDobl_Complex_Ring;              use OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with Generic_Matrices;

package OctoDobl_Complex_Matrices is
  new Generic_Matrices(OctoDobl_Complex_Ring,OctoDobl_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of octo double complex numbers.
