with DoblDobl_Complex_Series_Ring;
with DoblDobl_Complex_Series_Vectors;
with Generic_Matrices;

package DoblDobl_Complex_Series_Matrices is 
  new Generic_Matrices(DoblDobl_Complex_Series_Ring,
                       DoblDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of double double complex series.
