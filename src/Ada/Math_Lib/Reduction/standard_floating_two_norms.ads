with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;
with Standard_Floating_VecVecs;          use Standard_Floating_VecVecs;

package Standard_Floating_Two_Norms is

-- DESCRIPTION :
--   Provides implementation of 2-norm of standard-floating point vectors.

  function Norm2 ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector.

  procedure Normalize ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides every component of v by the 2-norm of v.

  procedure Normalize ( v : in out VecVec );

  -- DESCRIPTION :
  --   Normalizes every vector in v.

end Standard_Floating_Two_Norms;
