with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;

package Exponent_Indices is

-- DESCRIPTION :
--   An exponent index stores the indices of the nonzero entries
--   in an exponent vector.

  function Exponent_Index
             ( xp : Standard_Integer_Vectors.Vector )
             return Standard_Integer_Vectors.Vector;

  -- DESCRIPTION :
  --   Returns the vector of positions k in xp for which xp(k) = 1.

  function Exponent_Index
             ( xp : Standard_Integer_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec;

  -- DESCRIPTION :
  --   Returns the vector of exponent indices for xp.

end Exponent_Indices;
