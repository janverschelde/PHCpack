with Standard_Floating_Numbers;          use Standard_Floating_Numbers;
with Standard_Mathematical_Functions;

package body PentDobl_Mathematical_Functions is

  function SQRT ( x : penta_double ) return penta_double is

  --  Perform the following Newton iteration:
  --    y' = y + (1 - x * y^2) * y / 2;
  --  which converges to 1/sqrt(x), starting with the
  --  double precision approximation to 1/sqrt(x).
  --  Since Newton's iteration more or less doubles the
  --  number of correct digits, we only need to perform it twice.

    res,h : penta_double;
    f : double_float;
    one : constant penta_double := create(1.0);

  begin
    if is_zero(x) then
      res := Create(0.0);
    elsif is_negative(x) then  -- original code returns NaN
      res := Create(-1.0);
    else
      f := Standard_Mathematical_Functions.SQRT(thumb_part(x));
      res := Create(f);
      res := one/res;
      h := mul_pwr2(x,0.5);
      res := res + ((0.5 - h*(res*res))*res);
      res := res + ((0.5 - h*(res*res))*res);
      res := res + ((0.5 - h*(res*res))*res);
      res := res*x;
    end if;
    return res;
  end SQRT;

  function Radius ( x,y : penta_double ) return penta_double is
  begin
    return SQRT(x*x + y*y);
  end Radius;

end PentDobl_Mathematical_Functions;
