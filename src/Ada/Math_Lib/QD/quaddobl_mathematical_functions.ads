with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;

package QuadDobl_Mathematical_Functions is

-- DESCRIPTION :
--   This package provides some special mathematical functions for
--   quad double numbers.  Note that exp and log are already provided
--   in the package quad_double_numbers as needed for output.

  function SQRT ( x : double_float ) return quad_double;
  function SQRT ( x : double_double ) return quad_double;
  function SQRT ( x : quad_double ) return quad_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function SIN ( x : quad_double ) return quad_double;
  function COS ( x : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns sin(x) and cos(x).  Instead of NaN, returns -2.

  procedure SINCOS ( x : in quad_double; sin_x,cos_x : out quad_double );

  -- DESCRIPTION :
  --   Returns in sin_x and cos_x the sine and cosine of x respectively.

  function TAN ( x : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns the tangent of x as sin(x)/cos(x).

  function ARCTAN ( x : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns the inverse of the tangent of x.

  function ARCTAN ( y,x : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   if x equals 1, then returns ARCTAN(y) 
  --   else returns the angle.

  function ARCSIN ( x : quad_double ) return quad_double;
  function ARCCOS ( x : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Returns the inverse of sine and cosine function of x.

  function Radius ( x,y : quad_double ) return quad_double;
  function Angle ( x,y : quad_double ) return quad_double;

  -- DESCRIPTION :
  --   Radius and angle return the polar representation of a complex number.

end QuadDobl_Mathematical_Functions;
