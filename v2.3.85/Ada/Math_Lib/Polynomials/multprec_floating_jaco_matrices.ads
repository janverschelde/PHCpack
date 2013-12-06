with Multprec_Floating_Ring;
with Multprec_Floating_Vectors;
with Multprec_Floating_VecVecs;
with Multprec_Floating_Matrices;
with Multprec_Floating_Polynomials;
with Multprec_Floating_Poly_Functions;
with Multprec_Floating_Poly_Systems;
with Multprec_Floating_Poly_SysFun;
with Generic_Jacobian_Matrices;

package Multprec_Floating_Jaco_Matrices is
  new Generic_Jacobian_Matrices(Multprec_Floating_Ring,
                                Multprec_Floating_Vectors,
                                Multprec_Floating_VecVecs,
                                Multprec_Floating_Matrices,
                                Multprec_Floating_Polynomials,
                                Multprec_Floating_Poly_Functions,
                                Multprec_Floating_Poly_Systems,
                                Multprec_Floating_Poly_SysFun);

-- DESCRIPTION :
--   Defines functions for creating and evaluating Jacobian matrices for
--   systems of polynomials with multiprecision floating-point coefficients.
