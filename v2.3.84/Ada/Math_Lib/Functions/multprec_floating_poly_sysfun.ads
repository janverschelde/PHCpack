with Multprec_Floating_Ring;
with Multprec_Floating_Vectors;
with Multprec_Floating_VecVecs;
with Multprec_Floating_Polynomials;
with Multprec_Floating_Poly_Functions;
with Multprec_Floating_Poly_Systems;
with Generic_Poly_System_Functions;

package Multprec_Floating_Poly_SysFun is
  new Generic_Poly_System_Functions(Multprec_Floating_Ring,
                                    Multprec_Floating_Vectors,
                                    Multprec_Floating_VecVecs,
                                    Multprec_Floating_Polynomials,
                                    Multprec_Floating_Poly_Functions,
                                    Multprec_Floating_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials with
--   multi-precision real floating-point coefficients.
