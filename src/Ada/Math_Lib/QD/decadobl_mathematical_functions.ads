with Deca_Double_Numbers;                use Deca_Double_Numbers;

package DecaDobl_Mathematical_Functions is

-- DESCRIPTION :
--   Mathematical functions for deca double numbers.
--   The exp and log are already defined in Deca_Double_Numbers.

  function SQRT ( x : deca_double ) return deca_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function Radius ( x,y : deca_double ) return deca_double;

  -- DESCRIPTION :
  --   If x and y are the real and imaginary part of a complex number z
  --   then Radius(x,y) is the radius of z, computed as SQRT(x*x + y*y).

end DecaDobl_Mathematical_Functions;
