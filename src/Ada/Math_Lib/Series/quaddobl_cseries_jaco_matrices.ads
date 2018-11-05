with QuadDobl_Complex_Series_Ring;
with QuadDobl_Complex_Series_Vectors;
with QuadDobl_Complex_Series_VecVecs;
with QuadDobl_Complex_Series_Matrices;
with QuadDobl_CSeries_Polynomials;
with QuadDobl_CSeries_Poly_Functions;
with QuadDobl_CSeries_Poly_Systems;
with QuadDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package QuadDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(QuadDobl_Complex_Series_Ring,
                                QuadDobl_Complex_Series_Vectors,
                                QuadDobl_Complex_Series_VecVecs,
                                QuadDobl_Complex_Series_Matrices,
                                QuadDobl_CSeries_Polynomials,
                                QuadDobl_CSeries_Poly_Functions,
                                QuadDobl_CSeries_Poly_Systems,
                                QuadDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with quad double precision complex numbers.
