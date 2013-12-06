with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_VecVecs;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Standard_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Standard_Complex_Ring,
                                Standard_Complex_Vectors,
                                Standard_Complex_VecVecs,
                                Standard_Complex_Matrices,
                                Standard_Complex_Polynomials,
                                Standard_Complex_Poly_Functions,
                                Standard_Complex_Poly_Systems,
                                Standard_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the standard complex numbers.
