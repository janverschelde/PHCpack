with Multprec_Integer64_Ring;
with Multprec_Integer64_Ring.DDomain;
with Greatest_Common_Divisors;

package Multprec64_Common_Divisors is
  new Greatest_Common_Divisors(Multprec_Integer64_Ring,
                               Multprec_Integer64_Ring.DDomain);

-- DESCRIPTION :
--   Defines greatest common divisors over the euclidean domain of
--   multi-precision integer numbers, using 64-bit arithmetic.
