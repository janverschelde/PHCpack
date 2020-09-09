with PentDobl_Complex_Numbers;            use PentDobl_Complex_Numbers;
with PentDobl_Complex_Ring;
with PentDobl_Complex_Vectors;
with Generic_Dense_Series;

package PentDobl_Complex_Series is
  new Generic_Dense_Series(PentDobl_Complex_Ring,
                           PentDobl_Complex_Vectors,"/",Div);

-- DESCRIPTION :
--   Defines series using vectors over the ring of complex numbers,
--   with the division, in penta double precision.
