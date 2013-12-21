with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_Polynomials;
with Generic_Polynomial_Functions;

package Standard_Complex_Poly_Functions is
  new Generic_Polynomial_Functions(Standard_Complex_Ring,
                                   Standard_Complex_Vectors,
                                   Standard_Complex_Polynomials);

-- DESCRIPTION :
--   Defines polynomial functions for standard complex numbers.
