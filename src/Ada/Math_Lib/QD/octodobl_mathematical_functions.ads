with Octo_Double_Numbers;                use Octo_Double_Numbers;

package OctoDobl_Mathematical_Functions is

-- DESCRIPTION :
--   Mathematical functions for octo double numbers.
--   The exp and log are already defined in Octo_Double_Numbers.

  function SQRT ( x : octo_double ) return octo_double;

  -- DSECRIPTION :
  --   Returns the square root of x, if x >= 0.
  --   If x < 0, then -1 is returned.

  function Radius ( x,y : octo_double ) return octo_double;

  -- DESCRIPTION :
  --   If x and y are the real and imaginary part of a complex number z
  --   then Radius(x,y) is the radius of z, computed as SQRT(x*x + y*y).

end OctoDobl_Mathematical_Functions;
