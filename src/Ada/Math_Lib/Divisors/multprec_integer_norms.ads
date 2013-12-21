with Multprec_Integer_Numbers;           use Multprec_Integer_Numbers;
with Multprec_Integer_Vectors;

package Multprec_Integer_Norms is

-- DESCRIPTION :
--   This package defines the norm of an integer vector as the gcd of
--   all its entries.  An integer vector is normalized if its elements
--   are relative prime with respect to each other.

  function gcd ( v : Multprec_Integer_Vectors.Vector )
               return Integer_Number;

  -- DESCRIPTION :
  --   Returns the positive gcd(v(v'first),..,v(v'last)).

  procedure Normalize ( v : in out Multprec_Integer_Vectors.Vector );

  -- DESCRIPTION :
  --   All elements in the vector are divided by gcd(v).

end Multprec_Integer_Norms;
