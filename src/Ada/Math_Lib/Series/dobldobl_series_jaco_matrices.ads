with DoblDobl_Dense_Series_Ring;
with DoblDobl_Dense_Series_Vectors;
with DoblDobl_Dense_Series_VecVecs;
with DoblDobl_Dense_Series_Matrices;
with DoblDobl_Series_Polynomials;
with DoblDobl_Series_Poly_Functions;
with DoblDobl_Series_Poly_Systems;
with DoblDobl_Series_Poly_SysFun;
with Generic_Jacobian_Matrices;

package DoblDobl_Series_Jaco_Matrices is
  new Generic_Jacobian_Matrices(DoblDobl_Dense_Series_Ring,
                                DoblDobl_Dense_Series_Vectors,
                                DoblDobl_Dense_Series_VecVecs,
                                DoblDobl_Dense_Series_Matrices,
                                DoblDobl_Series_Polynomials,
                                DoblDobl_Series_Poly_Functions,
                                DoblDobl_Series_Poly_Systems,
                                DoblDobl_Series_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices
--   for systems of series polynomials over the complex numbers,
--   in double double precision.
