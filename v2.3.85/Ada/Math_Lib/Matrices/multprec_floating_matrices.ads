with Multprec_Floating_Ring;              use Multprec_Floating_Ring;
with Multprec_Floating_Vectors;
with Generic_Matrices;

package Multprec_Floating_Matrices is
  new Generic_Matrices(Multprec_Floating_Ring,
                       Multprec_Floating_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multi-precision floating-point numbers.
