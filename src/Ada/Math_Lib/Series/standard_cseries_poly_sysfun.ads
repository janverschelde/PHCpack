with Standard_Complex_Series_Ring;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Functions;
with Standard_CSeries_Poly_Systems;
with Generic_Poly_System_Functions;

package Standard_CSeries_Poly_SysFun is
  new Generic_Poly_System_Functions(Standard_Complex_Series_Ring,
                                    Standard_Complex_Series_Vectors,
                                    Standard_Complex_Series_VecVecs,
                                    Standard_CSeries_Polynomials,
                                    Standard_CSeries_Poly_Functions,
                                    Standard_CSeries_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of multivariate polynomials.
--   The polynomials have as coefficients truncated power series with
--   double precision complex numbers as coefficients.
