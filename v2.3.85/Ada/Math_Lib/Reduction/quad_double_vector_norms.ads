with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Quad_Double_Vectors;                use Quad_Double_Vectors;

package Quad_Double_Vector_Norms is

-- DESCRIPTION :
--   Provides function to compute various norms of vectors
--   of quad doubles.

  function Norm2 ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the complex vector v.

  function Sum_Norm ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Quad_Double_Vector_Norms;
