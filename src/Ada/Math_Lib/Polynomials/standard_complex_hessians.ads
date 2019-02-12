with Standard_Complex_Ring;
with Standard_Complex_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with Generic_Hessian_Matrices;

package Standard_Complex_Hessians is
  new Generic_Hessian_Matrices(Standard_Complex_Ring,
                               Standard_Complex_Vectors,
                               Standard_Complex_Matrices,
                               Standard_Complex_Polynomials,
                               Standard_Complex_Poly_Functions,
                               Standard_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Hessian matrices for
--   polynomials over the standard complex numbers.
