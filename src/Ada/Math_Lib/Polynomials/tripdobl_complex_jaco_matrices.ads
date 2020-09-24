with TripDobl_Complex_Ring;
with TripDobl_Complex_Vectors;
with TripDobl_Complex_VecVecs;
with TripDobl_Complex_Matrices;
with TripDobl_Complex_Polynomials;
with TripDobl_Complex_Poly_Functions;
with TripDobl_Complex_Poly_Systems;
with TripDobl_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package TripDobl_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(TripDobl_Complex_Ring,
                                TripDobl_Complex_Vectors,
                                TripDobl_Complex_VecVecs,
                                TripDobl_Complex_Matrices,
                                TripDobl_Complex_Polynomials,
                                TripDobl_Complex_Poly_Functions,
                                TripDobl_Complex_Poly_Systems,
                                TripDobl_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the triple double complex numbers.
