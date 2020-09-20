with DecaDobl_Complex_Ring;
with DecaDobl_Complex_Polynomials;
with Generic_Polynomial_Systems;

package DecaDobl_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(DecaDobl_Complex_Ring,
                                 DecaDobl_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the deca double complex numbers.
