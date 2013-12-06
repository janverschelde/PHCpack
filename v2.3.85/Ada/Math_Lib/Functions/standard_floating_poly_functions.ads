with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Standard_Floating_Polynomials;
with Generic_Polynomial_Functions;

package Standard_Floating_Poly_Functions is
  new Generic_Polynomial_Functions(Standard_Floating_Ring,
                                   Standard_Floating_Vectors,
                                   Standard_Floating_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for standard floating-point coefficients.
