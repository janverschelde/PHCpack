with Multprec_Floating_Numbers;          use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;           use Multprec_Complex_Numbers;
with Multprec_Complex_Vectors;           use Multprec_Complex_Vectors;

package Multprec_Complex_Norms_Equals is

-- DESCRIPTION :
--   Provides norms of vectors and decision routines for equalities,
--   for multi-precision complex numbers.

  function Norm2 ( v : Vector ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the 2-norm (the normal Euclidean norm) of 
  --   the complex vector v.

  function Max_Norm ( v : Vector ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the absolute value of the greatest element in v.

  function Sum_Norm ( v : Vector ) return Floating_Number;

  -- DESCRIPTION :
  --   Returns the sum of all absolute values of the elements in v.

  function Equal ( x,y : Complex_Number; tol : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if abs(x-y) < tol, otherwise false is returned.

  function Equal ( x,y : Vector; tol : Floating_Number ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Equal(x(i),y(i),tol), for i in x'range=y'range,
  --   otherwise false is returned.

  -- REQUIRED : x'range = y'range.

end Multprec_Complex_Norms_Equals;
