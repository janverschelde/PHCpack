with PentDobl_Complex_Series_Ring;
with PentDobl_Complex_Series_Vectors;
with Generic_Matrices;

package PentDobl_Complex_Series_Matrices is 
  new Generic_Matrices(PentDobl_Complex_Series_Ring,
                       PentDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of penta double complex series.
