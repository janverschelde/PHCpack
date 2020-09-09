with OctoDobl_Complex_Numbers;            use OctoDobl_Complex_Numbers;
with OctoDobl_Complex_Ring;
with OctoDobl_Complex_Vectors;
with Generic_Dense_Series;

package OctoDobl_Complex_Series is
  new Generic_Dense_Series(OctoDobl_Complex_Ring,
                           OctoDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in octo double precision.
