with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with QuadDobl_Complex_Numbers;           use QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;           use QuadDobl_Complex_Vectors;

package QuadDobl_Complex_Equality_Tests is

-- DESCRIPTION :
--   Provides routines to decide whether two quad double complex numbers
--   and vectors are equal up to some given tolerance.

  function Equal ( x,y : Complex_Number; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if abs(x-y) < tol, otherwise false is returned.

  function Equal ( x,y : Vector; tol : double_float ) return boolean;

  -- DESCRIPTION :
  --   Returns true if Equal(x(i),y(i),tol), for i in x'range=y'range,
  --   otherwise false is returned.

  -- REQUIRED : x'range = y'range.

end QuadDobl_Complex_Equality_Tests;
