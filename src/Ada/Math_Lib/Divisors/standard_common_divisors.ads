with Standard_Integer_Ring;
with Standard_Integer_Ring.DDomain;
with Greatest_Common_Divisors;

package Standard_Common_Divisors is
  new Greatest_Common_Divisors(Standard_Integer_Ring,
                               Standard_Integer_Ring.DDomain);

-- DESCRIPTION :
--   Defines greatest common divisors over the Euclidean domain of
--   standard integer numbers.
