with Penta_Double_Numbers;              use Penta_Double_Numbers;
with Penta_Double_Vectors;              use Penta_Double_Vectors;

package Penta_Double_Vector_Norms is

-- DESCRIPTION :
--   Functions to compute various norms of vectors of penta doubles.

  function Norm2 ( v : Vector ) return penta_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector v.

  function Sum_Norm ( v : Vector ) return penta_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return penta_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Penta_Double_Vector_Norms;
