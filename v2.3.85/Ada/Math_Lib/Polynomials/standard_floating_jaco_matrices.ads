with Standard_Floating_Ring;
with Standard_Floating_Vectors;
with Standard_Floating_VecVecs;
with Standard_Floating_Matrices;
with Standard_Floating_Polynomials;
with Standard_Floating_Poly_Functions;
with Standard_Floating_Poly_Systems;
with Standard_Floating_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Standard_Floating_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Standard_Floating_Ring,
                                Standard_Floating_Vectors,
                                Standard_Floating_VecVecs,
                                Standard_Floating_Matrices,
                                Standard_Floating_Polynomials,
                                Standard_Floating_Poly_Functions,
                                Standard_Floating_Poly_Systems,
                                Standard_Floating_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials with standard floating-point coefficients.
