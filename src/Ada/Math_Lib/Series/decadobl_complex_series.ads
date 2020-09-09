with DecaDobl_Complex_Numbers;            use DecaDobl_Complex_Numbers;
with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Vectors;
with Generic_Dense_Series;

package DecaDobl_Complex_Series is
  new Generic_Dense_Series(DecaDobl_Complex_Ring,
                           DecaDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in deca double precision.
