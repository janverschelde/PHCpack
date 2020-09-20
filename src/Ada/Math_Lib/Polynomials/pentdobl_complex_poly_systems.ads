with PentDobl_Complex_Ring;
with PentDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package PentDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(PentDobl_Complex_Ring,
                                 PentDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the penta double complex numbers.
