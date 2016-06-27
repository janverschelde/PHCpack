with QuadDobl_Dense_Series_Ring;
with QuadDobl_Dense_Series_Vectors;
with QuadDobl_Dense_Series_VecVecs;
with QuadDobl_Dense_Series_Matrices;
with QuadDobl_Series_Polynomials;
with QuadDobl_Series_Poly_Functions;
with QuadDobl_Series_Poly_Systems;
with QuadDobl_Series_Poly_SysFun;
with Generic_Jacobian_Matrices;

package QuadDobl_Series_Jaco_Matrices is
  new Generic_Jacobian_Matrices(QuadDobl_Dense_Series_Ring,
                                QuadDobl_Dense_Series_Vectors,
                                QuadDobl_Dense_Series_VecVecs,
                                QuadDobl_Dense_Series_Matrices,
                                QuadDobl_Series_Polynomials,
                                QuadDobl_Series_Poly_Functions,
                                QuadDobl_Series_Poly_Systems,
                                QuadDobl_Series_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices
--   for systems of series polynomials over the complex numbers,
--   in quad double precision.
