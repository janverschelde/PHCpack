with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Vectors;

package Handle_Underflow_Gracefully is

-- DESCRIPTION :
--   Applying Newton's method might give floating point numbers
--   that will cause underflow when used in future calculations.
--   The operations in this package are to be applied in case of
--   exceptions occurring during the application of Newton's method.

  function Underflow_to_Zero ( f : double_float ) return double_float;
  function Underflow_to_Zero ( f : double_double ) return double_double;
  function Underflow_to_Zero ( f : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns 0.0 if f + 1.0 = 1.0 holds, otherwise f is returned.

  procedure Underflow_to_Zero 
              ( x : in out Standard_Complex_Numbers.Complex_Number );
  procedure Underflow_to_Zero 
              ( x : in out DoblDobl_Complex_Numbers.Complex_Number );
  procedure Underflow_to_Zero 
              ( x : in out QuadDobl_Complex_Numbers.Complex_Number );

  -- DESCRIPTION :
  --   Applies the test f + 1.0 = 1.0 to both real and imaginary parts
  --   of the complex number x.

  procedure Underflow_to_Zero
              ( x : in out Standard_Complex_Vectors.Vector );
  procedure Underflow_to_Zero
              ( x : in out DoblDobl_Complex_Vectors.Vector );
  procedure Underflow_to_Zero
              ( x : in out QuadDobl_Complex_Vectors.Vector );

  -- DESCRIPTION :
  --   Sets floating-point numbers in all numbers of x to zero.

end Handle_Underflow_Gracefully;
