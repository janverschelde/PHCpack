with Standard_Integer_Numbers;          use Standard_Integer_Numbers;

package body Multprec_Univariate_Interpolators is

-- CONSTRUCTORS :

  function Create ( x,y : Vector ) return Vector is

    res : Vector(y'range);
    dif : Complex_Number;

  begin
    Copy(y,res);
    for i in 1..y'last loop
      for j in 0..(i-1) loop
        dif := x(j) - x(i);
        Min(res(i));
        Add(res(i),res(j));
        Div(res(i),dif);
        Clear(dif);
      end loop;
    end loop;
    return res;
  end Create;

  function Expand ( f,x : Vector ) return Vector is

    res : Vector(f'range);

  begin
    Copy(f(f'last),res(0));
    for i in reverse 0..f'last-1 loop
      Copy(res(f'last-i-1),res(f'last-i));  -- at degree f'last-i
      for j in reverse 1..f'last-i-1 loop   -- multiply with x-x(i)
        Mul(res(j),x(i));
        Min(res(j));
        Add(res(j),res(j-1));
      end loop;
      Mul(res(0),x(i));
      Min(res(0));
      Add(res(0),f(i));
    end loop;
    return res;
  end Expand;

-- EVALUATORS :

  function Evalf ( f,x : Vector; a : Complex_Number ) return Complex_Number is

    res,dif : Complex_Number;

  begin
    Copy(f(f'last),res);
    for i in reverse 0..f'last-1 loop
      dif := a - x(i);
      Mul(res,dif);
      Add(res,f(i));
      Clear(dif);
    end loop;
    return res;
  end Evalf;

  function Evalc ( c : Vector; x : Complex_Number ) return Complex_Number is

    res : Complex_Number;

  begin
    Copy(c(c'last),res);
    for i in reverse 0..c'last-1 loop
      Mul(res,x);
      Add(res,c(i));
    end loop;
    return res;
  end Evalc;

end Multprec_Univariate_Interpolators;
