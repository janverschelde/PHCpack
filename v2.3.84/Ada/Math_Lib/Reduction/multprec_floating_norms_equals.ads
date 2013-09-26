with Multprec_Floating_Ring;
with Multprec_Floating_Ring.FField;
with Multprec_Floating_Vectors;
with Multprec_Floating_Matrices;
with Generic_Norms_Equals;

package Multprec_Floating_Norms_Equals is
  new Generic_Norms_Equals(Multprec_Floating_Ring,
                           Multprec_Floating_Ring.FField,
                           Multprec_Floating_Vectors,
                           Multprec_Floating_Matrices);

-- DESCRIPTION :
--   Defines norms and equalities for multi-precision floating-point numbers.
