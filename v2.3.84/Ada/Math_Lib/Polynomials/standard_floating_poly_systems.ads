with Standard_Floating_Ring;
with Standard_Floating_Polynomials;
with Generic_Polynomial_Systems;

package Standard_Floating_Poly_Systems is
  new Generic_Polynomial_Systems(Standard_Floating_Ring,
                                 Standard_Floating_Polynomials);

-- DESCRIPTION :
--   Defines polynomial systems with double double coefficients.
