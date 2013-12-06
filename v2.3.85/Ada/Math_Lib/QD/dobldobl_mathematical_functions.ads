with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;

package DoblDobl_Mathematical_Functions is

-- DESCRIPTION :
--   This package provides some special mathematical functions for
--   double double numbers.  Note that exp and log are already provided
--   in the package double_double_numbers as needed for output.

  function SQRT ( x : double_float ) return double_double;
  function SQRT ( x : double_double ) return double_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function SIN ( x : double_double ) return double_double;
  function COS ( x : double_double ) return double_double;

  -- DESCRIPTION :
  --   Returns sine and cosine of x.
  --   Instead of NaN, returns -2.

  procedure SINCOS ( x : in double_double; 
                     sin_x,cos_x : out double_double );

  -- DESCRIPTION :
  --   Returns sin(x) and cos(x) in sin_x and cos_x.

  function TAN ( x : double_double ) return double_double;

  -- DESCRIPTION :
  --   Returns the tangent of x as sin(x)/cos(x).

  function ARCTAN ( x : double_double ) return double_double;

  -- DESCRIPTION :
  --   Returns the inverse of the tangent of x.

  function ARCTAN ( y,x : double_double ) return double_double;

  -- DESCRIPTION :
  --   If x equals 1, then returns ARCTAN(y)
  --   else returns the angle.

  function ARCSIN ( x : double_double ) return double_double;
  function ARCCOS ( x : double_double ) return double_double;

  -- DESCRIPTION :
  --   Returns the inverse of sine and cosine function.

  function Radius ( x,y : double_double ) return double_double;
  function Angle ( x,y : double_double ) return double_double;

  -- DESCRIPTION :
  --   Radius and angle return the polar representation of a complex number.

end DoblDobl_Mathematical_Functions;
