with Multprec_Floating_Ring;
with Generic_Vectors;

package Multprec_Floating_Vectors is
  new Generic_Vectors(Multprec_Floating_Ring);

-- DESCRIPTION :
--   Defines vectors over the ring of multi-precision floating-point numbers.
