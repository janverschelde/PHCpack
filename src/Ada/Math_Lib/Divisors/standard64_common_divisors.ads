with Standard_Integer64_Ring;
with Standard_Integer64_Ring.DDomain;
with Greatest_Common_Divisors;

package Standard64_Common_Divisors is
  new Greatest_Common_Divisors(Standard_Integer64_Ring,
                               Standard_Integer64_Ring.DDomain);

-- DESCRIPTION :
--   Defines greatest common divisors over the Euclidean domain of
--   standard integer numbers of type long long integer.
