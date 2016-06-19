with Standard_Dense_Series_Ring;
with Standard_Dense_Series_Vectors;
with Standard_Dense_Series_VecVecs;
with Standard_Dense_Series_Matrices;
with Standard_Series_Polynomials;
with Standard_Series_Poly_Functions;
with Standard_Series_Poly_Systems;
with Standard_Series_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Standard_Series_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Standard_Dense_Series_Ring,
                                Standard_Dense_Series_Vectors,
                                Standard_Dense_Series_VecVecs,
                                Standard_Dense_Series_Matrices,
                                Standard_Series_Polynomials,
                                Standard_Series_Poly_Functions,
                                Standard_Series_Poly_Systems,
                                Standard_Series_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of series polynomials over the standard complex numbers.
