with DoblDobl_Complex_Ring;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Systems;
with Generic_Hessian_Matrices;

package DoblDobl_Complex_Hessians is
  new Generic_Hessian_Matrices(DoblDobl_Complex_Ring,
                               DoblDobl_Complex_Vectors,
                               DoblDobl_Complex_Matrices,
                               DoblDobl_Complex_Polynomials,
                               DoblDobl_Complex_Poly_Functions,
                               DoblDobl_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Hessian matrices for
--   polynomials over the double double complex numbers.
