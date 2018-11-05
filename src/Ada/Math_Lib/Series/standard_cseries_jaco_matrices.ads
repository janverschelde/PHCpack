with Standard_Complex_Series_Ring;
with Standard_Complex_Series_Vectors;
with Standard_Complex_Series_VecVecs;
with Standard_Complex_Series_Matrices;
with Standard_CSeries_Polynomials;
with Standard_CSeries_Poly_Functions;
with Standard_CSeries_Poly_Systems;
with Standard_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Standard_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Standard_Complex_Series_Ring,
                                Standard_Complex_Series_Vectors,
                                Standard_Complex_Series_VecVecs,
                                Standard_Complex_Series_Matrices,
                                Standard_CSeries_Polynomials,
                                Standard_CSeries_Poly_Functions,
                                Standard_CSeries_Poly_Systems,
                                Standard_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with double precision complex numbers.
