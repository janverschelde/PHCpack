with Standard_Complex_Ring;
with Standard_Complex_Polynomials;
with Generic_Polynomial_Systems;

package Standard_Complex_Poly_Systems is
  new Generic_Polynomial_Systems(Standard_Complex_Ring,
                                 Standard_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems over the standard complex numbers.
