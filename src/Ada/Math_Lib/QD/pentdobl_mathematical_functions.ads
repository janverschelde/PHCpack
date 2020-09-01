with Penta_Double_Numbers;               use Penta_Double_Numbers;

package PentDobl_Mathematical_Functions is

-- DESCRIPTION :
--   Mathematical functions for penta double numbers.
--   The exp and log are already defined in penta_Double_Numbers.

  function SQRT ( x : penta_double ) return penta_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function Radius ( x,y : penta_double ) return penta_double;

  -- DESCRIPTION :
  --   If x and y are the real and imaginary part of a complex number z
  --   then Radius(x,y) is the radius of z, computed as SQRT(x*x + y*y).

end PentDobl_Mathematical_Functions;
