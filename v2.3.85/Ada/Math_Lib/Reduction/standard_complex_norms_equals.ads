with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;
with Standard_Complex_Matrices;          use Standard_Complex_Matrices;

package Standard_Complex_Norms_Equals is

-- DESCRIPTION :
--   Provides norms of vectors and decision routines for equalities,
--   for standard complex numbers.

  function Conjugated_Inner_Product ( v,w : Vector ) return Complex_Number;

  -- DESCRIPTION :
  --   Returns the complex inner product, taking complex conjugates
  --   of the components of v before multiplying with components of w.

  function Norm2 ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the 2-norm (the normal Euclidean norm) of 
  --   the complex vector v.

  procedure Normalize ( v : in out Vector );

  -- DESCRIPTION :
  --   Divides the vector v by its 2-norm, if its 2-norm is nonzero.

  function Max_Norm ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the absolute value of the greatest element in v.

  function Sum_Norm ( v : Vector ) return double_float;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the elements in v.

  function Max_Norm ( m : Matrix ) return double_float;

  -- DESCRIPTION :
  --   Returns the absolute value of the maximal element of m.

  function Equal ( x,y : Complex_Number; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if abs(x-y) < tol, otherwise false is returned.

  function Equal ( x,y : Vector; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Equal(x(i),y(i),tol), for i in x'range=y'range,
  --   otherwise false is returned.

  -- REQUIRED : x'range = y'range.

end Standard_Complex_Norms_Equals;
