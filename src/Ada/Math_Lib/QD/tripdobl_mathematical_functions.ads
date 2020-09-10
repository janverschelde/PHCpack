with Triple_Double_Numbers;              use Triple_Double_Numbers;

package TripDobl_Mathematical_Functions is

-- DESCRIPTION :
--   Mathematical functions for triple double numbers.
--   The exp and log are already defined in Triple_Double_Numbers.

  function SQRT ( x : triple_double ) return triple_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function SIN ( x : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   Returns sin(x).

  function COS ( x : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   Returns cos(x).

  procedure SINCOS ( x : in triple_double; sin_x,cos_x : out triple_double );

  -- DESCRIPTION :
  --   Returns in sin_x and cos_x the sine and cosine of x respectively.

  function TAN ( x : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   Returns the tangent of x as sin(x)/cos(x).

  function ARCTAN ( x : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   Returns the inverse of the tangent of x.

  function ARCTAN ( y,x : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   if x equals 1, then returns ARCTAN(y) 
  --   else returns the angle.

  function ARCSIN ( x : triple_double ) return triple_double;
  function ARCCOS ( x : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   Returns the inverse of sine and cosine function of x.

  function Radius ( x,y : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   If x and y are the real and imaginary part of a complex number z
  --   then Radius(x,y) is the radius of z, computed as SQRT(x*x + y*y).

  function Angle ( x,y : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   Returns the angle in the polar representation of the complex
  --   number with real part in x and imaginary part in y.

end TripDobl_Mathematical_Functions;
