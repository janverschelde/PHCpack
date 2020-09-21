with DecaDobl_Complex_Series_Ring;
with DecaDobl_Complex_Series_Vectors;
with DecaDobl_Complex_Series_VecVecs;
with DecaDobl_Complex_Series_Matrices;
with DecaDobl_CSeries_Polynomials;
with DecaDobl_CSeries_Poly_Functions;
with DecaDobl_CSeries_Poly_Systems;
with DecaDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package DecaDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(DecaDobl_Complex_Series_Ring,
                                DecaDobl_Complex_Series_Vectors,
                                DecaDobl_Complex_Series_VecVecs,
                                DecaDobl_Complex_Series_Matrices,
                                DecaDobl_CSeries_Polynomials,
                                DecaDobl_CSeries_Poly_Functions,
                                DecaDobl_CSeries_Poly_Systems,
                                DecaDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with deca double precision complex numbers.
