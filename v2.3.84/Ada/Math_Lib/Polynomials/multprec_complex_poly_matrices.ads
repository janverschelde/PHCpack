with Multprec_Complex_Poly_Ring;         use Multprec_Complex_Poly_Ring;
with Multprec_Complex_Poly_Vectors;
with Generic_Matrices;

package Multprec_Complex_Poly_Matrices is
  new Generic_Matrices(Multprec_Complex_Poly_Ring,
                       Multprec_Complex_Poly_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multiprecision complex polynomials.
