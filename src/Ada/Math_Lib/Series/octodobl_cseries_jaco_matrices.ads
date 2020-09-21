with OctoDobl_Complex_Series_Ring;
with OctoDobl_Complex_Series_Vectors;
with OctoDobl_Complex_Series_VecVecs;
with OctoDobl_Complex_Series_Matrices;
with OctoDobl_CSeries_Polynomials;
with OctoDobl_CSeries_Poly_Functions;
with OctoDobl_CSeries_Poly_Systems;
with OctoDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package OctoDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(OctoDobl_Complex_Series_Ring,
                                OctoDobl_Complex_Series_Vectors,
                                OctoDobl_Complex_Series_VecVecs,
                                OctoDobl_Complex_Series_Matrices,
                                OctoDobl_CSeries_Polynomials,
                                OctoDobl_CSeries_Poly_Functions,
                                OctoDobl_CSeries_Poly_Systems,
                                OctoDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with octo double precision complex numbers.
