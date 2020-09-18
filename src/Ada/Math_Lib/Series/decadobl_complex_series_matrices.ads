with DecaDobl_Complex_Series_Ring;
with DecaDobl_Complex_Series_Vectors;
with Generic_Matrices;

package DecaDobl_Complex_Series_Matrices is 
  new Generic_Matrices(DecaDobl_Complex_Series_Ring,
                       DecaDobl_Complex_Series_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of deca double complex series.
