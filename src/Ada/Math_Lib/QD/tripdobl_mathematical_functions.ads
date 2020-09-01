with Triple_Double_Numbers;              use Triple_Double_Numbers;

package TripDobl_Mathematical_Functions is

-- DESCRIPTION :
--   Mathematical functions for triple double numbers.
--   The exp and log are already defined in Triple_Double_Numbers.

  function SQRT ( x : triple_double ) return triple_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function Radius ( x,y : triple_double ) return triple_double;

  -- DESCRIPTION :
  --   If x and y are the real and imaginary part of a complex number z
  --   then Radius(x,y) is the radius of z, computed as SQRT(x*x + y*y).

end TripDobl_Mathematical_Functions;
