with QuadDobl_Complex_Ring;
with QuadDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package QuadDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(QuadDobl_Complex_Ring,
                                 QuadDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the quad double complex numbers.
