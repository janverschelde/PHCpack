with Double_Double_Numbers;              use Double_Double_Numbers;
with DoblDobl_Complex_Numbers;           use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;           use DoblDobl_Complex_Vectors;

package DoblDobl_Complex_Vector_Norms is

-- DESCRIPTION :
--   Provides function to compute various norms of vectors
--   of complex double doubles.

  function Conjugated_Inner_Product ( v,w : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the complex inner product, taking complex conjugates
  --   of the components of v before multiplying with components of w.

  function Norm2 ( v : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the 2-norm of the complex vector v.

  procedure Normalize ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides the vector v by its 2-norm, if its 2-norm is nonzero.

  function Sum_Norm ( v : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the components of v.

  function Max_Norm ( v : Vector ) return double_double;

  -- DESCRIPTION :
  --   Returns the absolute value of the largest element in v.

end DoblDobl_Complex_Vector_Norms;
