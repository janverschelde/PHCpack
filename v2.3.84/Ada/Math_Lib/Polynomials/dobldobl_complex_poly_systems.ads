with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package DoblDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(DoblDobl_Complex_Ring,
                                 DoblDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the double double complex numbers.
