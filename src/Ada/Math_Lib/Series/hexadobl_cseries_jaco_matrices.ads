with HexaDobl_Complex_Series_Ring;
with HexaDobl_Complex_Series_Vectors;
with HexaDobl_Complex_Series_VecVecs;
with HexaDobl_Complex_Series_Matrices;
with HexaDobl_CSeries_Polynomials;
with HexaDobl_CSeries_Poly_Functions;
with HexaDobl_CSeries_Poly_Systems;
with HexaDobl_CSeries_Poly_SysFun;
with Generic_Jacobian_Matrices;

package HexaDobl_CSeries_Jaco_Matrices is
  new Generic_Jacobian_Matrices(HexaDobl_Complex_Series_Ring,
                                HexaDobl_Complex_Series_Vectors,
                                HexaDobl_Complex_Series_VecVecs,
                                HexaDobl_Complex_Series_Matrices,
                                HexaDobl_CSeries_Polynomials,
                                HexaDobl_CSeries_Poly_Functions,
                                HexaDobl_CSeries_Poly_Systems,
                                HexaDobl_CSeries_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials in several variables over the field of
--   truncated power series with hexa double precision complex numbers.
