with QuadDobl_Complex_Numbers;            use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with Generic_Dense_Series;

package QuadDobl_Complex_Series is
  new Generic_Dense_Series(QuadDobl_Complex_Ring,
                           QuadDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in quad double precision.
