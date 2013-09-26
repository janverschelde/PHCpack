with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;

package Multprec_Common_Divisors is

-- DESCRIPTION :
--   This package contains routines for the computation
--   of the greatest common divisor and least common multiple,
--   for multiprecision integer numbers.

  function gcd ( a,b : Integer_Number ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the greatest common divisor of a and b.

  function lcm ( a,b : Integer_Number ) return Integer_Number;

  -- DESCRIPTION :
  --   Returns the least common multiple of a and b.

  procedure gcd ( a,b : in Integer_Number; k,l,d : out Integer_Number );

  -- DESCRIPTION :
  --   Computes the greatest common divisor d of a and b;
  --   After gcd(a,b,k,l,d), there holds: k*a + l*b = d.

end Multprec_Common_Divisors;
