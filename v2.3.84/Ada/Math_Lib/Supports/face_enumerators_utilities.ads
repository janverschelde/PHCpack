with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Vectors;           use Standard_Integer_Vectors;
with Standard_Integer_VecVecs;           use Standard_Integer_VecVecs;

package Face_Enumerators_Utilities is

-- DESCRIPTION :
--   This package contains utilities for the face enumerators.

  function Is_Zero ( v : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the given vector equals the zero vector.

  function gcd ( v : Vector ) return integer32;

  -- DESCRIPTION :
  --   Returns the greatest common divisor gcd(v(v'first),..,v(v'last)).

  procedure Scale ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides each element in v by gcd(v), when gcd(v) /= 0 of course.

  function Is_In ( x : integer32; v : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if there exists an entry in v, say v(k),
  --   such that v(k) = x.

  function In_Edge ( x,a,b : Vector ) return boolean;

  -- DESCRIPTION :
  --   Returns true if the vector x lies between a and b.

  function In_Face ( k : in integer32; face,x : Vector; pts : VecVec )
                   return boolean;

  -- DESCRIPTION :
  --   Returns true if x lies in the given k-face, which contains entries
  --   to its elements in pts.

end Face_Enumerators_Utilities;
