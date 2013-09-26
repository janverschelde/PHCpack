with Multprec_Complex_Ring;
with Multprec_Complex_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Functions;
with Multprec_Complex_Poly_Systems;
with Generic_Poly_System_Functions;

package Multprec_Complex_Poly_SysFun is
  new Generic_Poly_System_Functions(Multprec_Complex_Ring,
                                    Multprec_Complex_Vectors,
                                    Multprec_Complex_VecVecs,
                                    Multprec_Complex_Polynomials,
                                    Multprec_Complex_Poly_Functions,
                                    Multprec_Complex_Poly_Systems);

-- DESCRIPTION :
--   Defines functions for evaluating systems of polynomials over the
--   multi-precision complex numbers.
