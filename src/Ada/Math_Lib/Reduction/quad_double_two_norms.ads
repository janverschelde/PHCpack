with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Vectors;                use Quad_Double_Vectors;
with Quad_Double_VecVecs;                use Quad_Double_VecVecs;

package Quad_Double_Two_Norms is

-- DESCRIPTION :
--   Provides implementation of 2-norm of double double vectors.

  function Norm2 ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector.

  procedure Normalize ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides every component of v by the 2-norm of v.

  procedure Normalize ( v : in out VecVec );

  -- DESCRIPTION :
  --   Normalizes every vector in v.

end Quad_Double_Two_Norms;
