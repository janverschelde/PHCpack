with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Multprec_Complex_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Multprec_Complex_Ring,
                                Multprec_Complex_Vectors,
                                Multprec_Complex_VecVecs,
                                Multprec_Complex_Matrices,
                                Multprec_Complex_Polynomials,
                                Multprec_Complex_Poly_Functions,
                                Multprec_Complex_Poly_Systems,
                                Multprec_Complex_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials over the multi-precision complex numbers.
