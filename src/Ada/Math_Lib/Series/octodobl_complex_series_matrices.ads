with OctoDobl_Complex_Series_Ring;
with OctoDobl_Complex_Series_Vectors;
with Generic_Matrices;

package OctoDobl_Complex_Series_Matrices is 
  new Generic_Matrices(OctoDobl_Complex_Series_Ring,
                       OctoDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of octo double complex series.
