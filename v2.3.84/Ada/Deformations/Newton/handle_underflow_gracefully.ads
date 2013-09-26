with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Vectors;           use Standard_Complex_Vectors;

package Handle_Underflow_Gracefully is

-- DESCRIPTION :
--   Applying Newton's method might give floating point numbers
--   that will cause underflow when used in future calculations.
--   The operations in this package are to be applied in case of
--   exceptions occurring during the application of Newton's method.

  function Underflow_to_Zero ( f : double_float ) return double_float;

  -- DESCRIPTION :
  --   Returns 0.0 if f + 1.0 = 1.0 holds, otherwise f is returned.

  procedure Underflow_to_Zero ( x : in out Complex_Number );

  -- DESCRIPTION :
  --   Applies the test f + 1.0 = 1.0 to both real and imaginary parts
  --   of the complex number x.

  procedure Underflow_to_Zero ( x : in out Vector );

  -- DESCRIPTION :
  --   Sets floating-point numbers in all numbers of x to zero.

end Handle_Underflow_Gracefully;
