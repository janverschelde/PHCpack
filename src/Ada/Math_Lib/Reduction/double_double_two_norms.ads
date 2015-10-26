with Double_Double_Numbers;              use Double_Double_Numbers;
with Double_Double_Vectors;              use Double_Double_Vectors;
with Double_Double_VecVecs;              use Double_Double_VecVecs;

package Double_Double_Two_Norms is

-- DESCRIPTION :
--   Provides implementation of 2-norm of double double vectors.

  function Norm2 ( v : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector.

  procedure Normalize ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides every component of v by the 2-norm of v.

  procedure Normalize ( v : in out VecVec );

  -- DESCRIPTION :
  --   Normalizes every vector in v.

end Double_Double_Two_Norms;
