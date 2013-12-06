with Quad_Double_Numbers;                use Quad_Double_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;

package QuadDobl_Complex_Vector_Norms is

-- DESCRIPTION :
--   Provides function to compute various norms of vectors
--   of complex quad doubles.

  function Norm2 ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the complex vector v.

  function Sum_Norm ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return quad_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end QuadDobl_Complex_Vector_Norms;
