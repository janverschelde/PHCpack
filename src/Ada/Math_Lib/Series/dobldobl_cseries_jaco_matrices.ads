with DoblDobl_Complex_Series_Ring;
with DoblDobl_Complex_Series_Vectors;
with DoblDobl_Complex_Series_VecVecs;
with DoblDobl_Complex_Series_Matrices;
with DoblDobl_CSeries_Polynomials;
with DoblDobl_CSeries_Poly_Functions;
with DoblDobl_CSeries_Poly_Systems;
with DoblDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package DoblDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(DoblDobl_Complex_Series_Ring,
                                DoblDobl_Complex_Series_Vectors,
                                DoblDobl_Complex_Series_VecVecs,
                                DoblDobl_Complex_Series_Matrices,
                                DoblDobl_CSeries_Polynomials,
                                DoblDobl_CSeries_Poly_Functions,
                                DoblDobl_CSeries_Poly_Systems,
                                DoblDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with double double precision complex numbers.
