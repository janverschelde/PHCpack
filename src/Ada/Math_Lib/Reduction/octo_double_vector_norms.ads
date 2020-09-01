with Octo_Double_Numbers;              use Octo_Double_Numbers;
with Octo_Double_Vectors;              use Octo_Double_Vectors;

package Octo_Double_Vector_Norms is

-- DESCRIPTION :
--   Functions to compute various norms of vectors of octo doubles.

  function Norm2 ( v : Vector ) return octo_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector v.

  function Sum_Norm ( v : Vector ) return octo_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return octo_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Octo_Double_Vector_Norms;
