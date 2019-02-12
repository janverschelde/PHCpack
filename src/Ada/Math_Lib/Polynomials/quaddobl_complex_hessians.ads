with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Systems;
with Generic_Hessian_Matrices;

package QuadDobl_Complex_Hessians is
  new Generic_Hessian_Matrices(QuadDobl_Complex_Ring,
                               QuadDobl_Complex_Vectors,
                               QuadDobl_Complex_Matrices,
                               QuadDobl_Complex_Polynomials,
                               QuadDobl_Complex_Poly_Functions,
                               QuadDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Hessian matrices for
--   polynomials over the quad double complex numbers.
