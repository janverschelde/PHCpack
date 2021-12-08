with Hexa_Double_Numbers;                use Hexa_Double_Numbers;

package HexaDobl_Mathematical_Functions is

-- DESCRIPTION :
--   Mathematical functions for hexa double numbers.
--   The exp and log are already defined in Hexa_Double_Numbers.

  function SQRT ( x : hexa_double ) return hexa_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function SIN ( x : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Returns sin(x).

  function COS ( x : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Returns cos(x).

  procedure SINCOS ( x : in hexa_double; sin_x,cos_x : out hexa_double );

  -- DESCRIPTION :
  --   Returns in sin_x and cos_x the sine and cosine of x respectively.

  function TAN ( x : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the tangent of x as sin(x)/cos(x).

  function ARCTAN ( x : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the inverse of the tangent of x.

  function ARCTAN ( y,x : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   if x equals 1, then returns ARCTAN(y) 
  --   else returns the angle.

  function ARCSIN ( x : hexa_double ) return hexa_double;
  function ARCCOS ( x : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the inverse of sine and cosine function of x.

  function Radius ( x,y : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   If x and y are the real and imaginary part of a complex number z
  --   then Radius(x,y) is the radius of z, computed as SQRT(x*x + y*y).

  function Angle ( x,y : hexa_double ) return hexa_double;

  -- DESCRIPTION :
  --   Returns the angle in the polar representation of the complex
  --   number with real part in x and imaginary part in y.

end HexaDobl_Mathematical_Functions;
