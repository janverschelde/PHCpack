with QuadDobl_Complex_Series_Ring;
with QuadDobl_Complex_Series_Vectors;
with Generic_Matrices;

package QuadDobl_Complex_Series_Matrices is 
  new Generic_Matrices(QuadDobl_Complex_Series_Ring,
                       QuadDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of quad double complex series.
