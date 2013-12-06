with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Floating_Vectors;          use Standard_Floating_Vectors;

package Standard_Floating_Vector_Norms is

-- DESCRIPTION :
--   Provides function to compute various norms of vectors
--   of double doubles.

  function Norm2 ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the 2-norm of the complex vector v.

  function Sum_Norm ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Standard_Floating_Vector_Norms;
