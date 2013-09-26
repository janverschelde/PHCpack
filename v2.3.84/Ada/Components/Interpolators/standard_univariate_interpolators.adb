with Standard_Integer_Numbers;           use Standard_Integer_Numbers;

package body Standard_Univariate_Interpolators is

-- CONSTRUCTORS :

  function Create ( x,y : Vector ) return Vector is

    res : Vector(y'range) := y;

  begin
    for i in 1..y'last loop
      for j in 0..(i-1) loop
        res(i) := (res(j) - res(i))/(x(j) - x(i));
      end loop;
    end loop;
    return res;
  end Create;

  function Expand ( f,x : Vector ) return Vector is

    res : Vector(f'range);

  begin
    res(0) := f(f'last);
    for i in reverse 0..f'last-1 loop
      res(f'last-i) := res(f'last-i-1);     -- at degree f'last-i
      for j in reverse 1..f'last-i-1 loop   -- multiply with x-x(i)
        res(j) := res(j-1) - x(i)*res(j);
      end loop;
      res(0) := -x(i)*res(0) + f(i);
    end loop;
    return res;
  end Expand;

-- EVALUATORS :

  function Evalf ( f,x : Vector; a : Complex_Number ) return Complex_Number is

    res : Complex_Number := f(f'last);

  begin
    for i in reverse 0..f'last-1 loop
      res := res*(a - x(i));
      res := res + f(i);
    end loop;
    return res;
  end Evalf;

  function Evalc ( c : Vector; x : Complex_Number ) return Complex_Number is

    res : Complex_Number := c(c'last);

  begin
    for i in reverse 0..c'last-1 loop
      res := res*x;
      res := res + c(i);
    end loop;
    return res;
  end Evalc;

end Standard_Univariate_Interpolators;
