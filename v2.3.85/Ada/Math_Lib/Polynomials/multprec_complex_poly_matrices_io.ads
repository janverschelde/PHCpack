with Multprec_Complex_Poly_Ring_io;
with Multprec_Complex_Poly_Vectors;
with Multprec_Complex_Poly_Matrices;
with Generic_Matrices_io;

package Multprec_Complex_Poly_Matrices_io is 
  new Generic_Matrices_io(Multprec_Complex_Poly_Ring_io,
                          Multprec_Complex_Poly_Vectors,
                          Multprec_Complex_Poly_Matrices);

-- DESCRIPTION :
--   Defines input/output of matrices of multiprecision complex polynomials.
