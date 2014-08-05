with QuadDobl_Complex_Poly_Ring_io;
with QuadDobl_Complex_Poly_Vectors;
with Generic_Vectors_io;

package QuadDobl_Complex_Poly_Vectors_io is 
  new Generic_Vectors_io(QuadDobl_Complex_Poly_Ring_io,
                         QuadDobl_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of quad double complex polynomials.
