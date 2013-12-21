with Multprec_Complex_Poly_Ring_io;
with Multprec_Complex_Poly_Vectors;
with Generic_Vectors_io;

package Multprec_Complex_Poly_Vectors_io is 
  new Generic_Vectors_io(Multprec_Complex_Poly_Ring_io,
                         Multprec_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines input/output of vectors of multiprecision complex polynomials.
