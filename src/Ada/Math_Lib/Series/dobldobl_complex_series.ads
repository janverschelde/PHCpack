with DoblDobl_Complex_Numbers;            use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with Generic_Dense_Series;

package DoblDobl_Complex_Series is
  new Generic_Dense_Series(DoblDobl_Complex_Ring,
                           DoblDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in double double precision.
