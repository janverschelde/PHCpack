with Standard_Dense_Series_Ring;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_VecVecs;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
with Standard_Series_Poly_Systems;
with Generic_Poly_System_Functions;

package Standard_Series_Poly_SysFun is
  new Generic_Poly_System_Functions(Standard_Dense_Series_Ring,
                                    Standard_Dense_Series_Vectors,
                                    Standard_Dense_Series_VecVecs,
                                    Standard_Series_Polynomials,
                                    Standard_Series_Poly_Functions,
                                    Standard_Series_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   standard complex numbers as coefficients.
