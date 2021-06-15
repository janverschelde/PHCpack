with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Standard_Univariate_Interpolators;

package body Standard_Probe_Kernel is

  function Maximal_Degree ( p : Poly_Sys ) return integer32 is

    res : integer32 := Degree(p(p'first));
    dpi : integer32;

  begin
    for i in p'first+1..p'last loop
      dpi := Degree(p(i));
      if dpi > res
       then res := dpi;
      end if;
    end loop;
    return dpi;
  end Maximal_Degree;

  function Random_Vector_in_Kernel
              ( V : Matrix; corank : natural32 ) return Vector is

    res : Vector(V'range(1));
    r : Complex_Number;

  begin
    for i in res'range loop       -- initialize with last column in V
      res(i) := V(i,V'last(2));
    end loop;
    for j in 1..integer32(corank)-1 loop  -- add the (j+1)-to-last column of V
      r := Random1;                       -- multiplied with a random number
      for i in res'range loop
        res(i) := res(i) + r*V(i,V'last(2)-j);
      end loop;
    end loop;
    return res;
  end Random_Vector_in_Kernel;

  function Sample_Sum_on_Line
              ( f : Eval_Poly_Sys; z,w,t : Vector ) return Vector is

    res : Vector(t'range);
    x : Vector(z'range);
    y : Vector(f'range);

  begin
    for i in t'range loop
      x := z + t(i)*w;
      y := Eval(f,x);
      res(i) := Sum(y);
    end loop;
    return res;
  end Sample_Sum_on_Line;

  function Interpolation_Coefficients ( x,y : Vector ) return Vector is

    f : Vector(x'range) := Standard_Univariate_Interpolators.Create(x,y);

  begin
    return Standard_Univariate_Interpolators.Expand(f,x);
  end Interpolation_Coefficients;

  function Numerical_Order
              ( c : Vector; tol : double_float ) return natural32 is
  begin
    for i in c'range loop
      if AbsVal(c(i)) > tol
       then return natural32(i);
      end if;
    end loop;
    return natural32(c'last)+1;
  end Numerical_Order;

end Standard_Probe_Kernel;
