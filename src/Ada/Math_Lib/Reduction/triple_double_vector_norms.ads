with Triple_Double_Numbers;              use Triple_Double_Numbers;
with Triple_Double_Vectors;              use Triple_Double_Vectors;

package Triple_Double_Vector_Norms is

-- DESCRIPTION :
--   Functions to compute various norms of vectors of triple doubles.

  function Norm2 ( v : Vector ) return triple_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector v.

  function Sum_Norm ( v : Vector ) return triple_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return triple_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Triple_Double_Vector_Norms;
