with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Complex_Numbers;          use Multprec_Complex_Numbers;
with Multprec_Random_Numbers;           use Multprec_Random_Numbers;
with Multprec_Complex_Polynomials;      use Multprec_Complex_Polynomials;
with Multprec_Univariate_Interpolators;

package body Multprec_Probe_Kernel is

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
              ( V : Matrix; corank,s : natural32 ) return Vector is

    res : Vector(V'range(1));
    r,t : Complex_Number;

  begin
    for i in res'range loop       -- initialize with last column in V
      Copy(V(i,V'last(2)),res(i));
    end loop;
    for j in 1..integer32(corank)-1 loop -- add the (j+1)-to-last column of V
      r := Random(s);                    -- multiplied with a random number
      for i in res'range loop
        t := r*V(i,V'last(2)-j);
        Add(res(i),t);
        Clear(t);
      end loop;
      Clear(r);
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
      x := t(i)*w;
      Add(x,z);
      y := Eval(f,x);
      res(i) := Sum(y);
      Clear(x); Clear(y);
    end loop;
    return res;
  end Sample_Sum_on_Line;

  function Interpolation_Coefficients ( x,y : Vector ) return Vector is

    res : Vector(y'range);
    f : Vector(x'range) := Multprec_Univariate_Interpolators.Create(x,y);

  begin
    res := Multprec_Univariate_Interpolators.Expand(f,x);
    Clear(f);
    return res;
  end Interpolation_Coefficients;

  function Numerical_Order
              ( c : Vector; tol : double_float ) return natural32 is

    res : natural32 := natural32(c'last)+1;
    v : Floating_Number;

  begin
    for i in c'range loop
      v := AbsVal(c(i));
      if v > tol
       then res := natural32(i);
      end if;
      Clear(v);
      exit when (res < natural32(c'last)+1);
    end loop;
    return res;
  end Numerical_Order;

end Multprec_Probe_Kernel;
