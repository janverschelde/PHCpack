with Multprec_Complex_Ring;              use Multprec_Complex_Ring;
with Multprec_Complex_Vectors;
with Generic_Matrices;

package Multprec_Complex_Matrices is
  new Generic_Matrices(Multprec_Complex_Ring,
                       Multprec_Complex_Vectors);

-- DESCRIPTION :
--   Defines matrices over the ring of multi-precision complex numbers.
