with Hexa_Double_Numbers;                use Hexa_Double_Numbers;
with Hexa_Double_Vectors;                use Hexa_Double_Vectors;

package Hexa_Double_Vector_Norms is

-- DESCRIPTION :
--   Functions to compute various norms of vectors of hexa doubles.

  function Norm2 ( v : Vector ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector v.

  function Sum_Norm ( v : Vector ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Hexa_Double_Vector_Norms;
