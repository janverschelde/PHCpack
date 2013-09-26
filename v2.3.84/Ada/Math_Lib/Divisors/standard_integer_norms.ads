with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer64_Vectors;

package Standard_Integer_Norms is

-- DESCRIPTION :
--   This package defines the norm of an integer vector as the gcd of
--   all its entries.  An integer vector is normalized if its elements
--   are relative prime with respect to each other.

  function gcd ( v : Standard_Integer_Vectors.Vector ) return integer32;
  function gcd ( v : Standard_Integer64_Vectors.Vector ) return integer64;

  -- DESCRIPTION :
  --   Returns the positive gcd(v(v'first),..,v(v'last)).

  procedure Normalize ( v : in out Standard_Integer_Vectors.Vector );
  procedure Normalize ( v : in out Standard_Integer64_Vectors.Vector );

  -- DESCRIPTION :
  --   All elements in the vector are divided by gcd(v).

end Standard_Integer_Norms;
