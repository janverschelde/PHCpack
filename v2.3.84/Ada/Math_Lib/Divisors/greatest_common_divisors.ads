with Abstract_Ring;
with Abstract_Ring.Domain;

generic

  with package Ring is new Abstract_Ring(<>);
  with package Euclidean_Domain is new Ring.Domain(<>);

package Greatest_Common_Divisors is

-- DESCRIPTION :
--   This package contains routines for the computation
--   of the greatest common divisor and least common multiple.

  use Ring;
  use Euclidean_Domain;

  function gcd ( a,b : number ) return number;

  -- DESCRIPTION :
  --   Returns the greatest common divisor of a and b.

  function lcm ( a,b : number ) return number;

  -- DESCRIPTION :
  --   Returns the least common multiple of a and b.

  procedure gcd ( a,b : in number; k,l,d : out number );

  -- DESCRIPTION :
  --   Computes the greatest common divisor d of a and b;
  --   After gcd(a,b,k,l,d), there holds: k*a + l*b = d.

end Greatest_Common_Divisors;
