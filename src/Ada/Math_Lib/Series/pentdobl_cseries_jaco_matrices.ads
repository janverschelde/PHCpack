with PentDobl_Complex_Series_Ring;
with PentDobl_Complex_Series_Vectors;
with PentDobl_Complex_Series_VecVecs;
with PentDobl_Complex_Series_Matrices;
with PentDobl_CSeries_Polynomials;
with PentDobl_CSeries_Poly_Functions;
with PentDobl_CSeries_Poly_Systems;
with PentDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package PentDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(PentDobl_Complex_Series_Ring,
                                PentDobl_Complex_Series_Vectors,
                                PentDobl_Complex_Series_VecVecs,
                                PentDobl_Complex_Series_Matrices,
                                PentDobl_CSeries_Polynomials,
                                PentDobl_CSeries_Poly_Functions,
                                PentDobl_CSeries_Poly_Systems,
                                PentDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with penta double precision complex numbers.
