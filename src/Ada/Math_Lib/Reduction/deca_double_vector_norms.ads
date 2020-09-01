with Deca_Double_Numbers;                use Deca_Double_Numbers;
with Deca_Double_Vectors;                use Deca_Double_Vectors;

package Deca_Double_Vector_Norms is

-- DESCRIPTION :
--   Functions to compute various norms of vectors of deca doubles.

  function Norm2 ( v : Vector ) return deca_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the vector v.

  function Sum_Norm ( v : Vector ) return deca_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return deca_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end Deca_Double_Vector_Norms;
