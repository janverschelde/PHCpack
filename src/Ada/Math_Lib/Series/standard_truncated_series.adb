with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Norms_Equals;

package body Standard_Truncated_Series is

  function "*" ( a,b : Vector ) return Vector is

    res : Vector(0..a'last);

  begin
    for i in res'range loop
      res(i) := a(0)*b(i);
      for j in 1..i loop
        res(i) := res(i) + a(j)*b(i-j);
      end loop;
    end loop;
    return res;
  end "*";

  function Inverse ( c : Vector ) return Vector is

    res : Vector(0..c'last);

  begin
    res(0) := 1.0/c(0);
    for i in 1..res'last loop
      res(i) := -c(1)*res(i-1);
      for j in 2..i loop
        res(i) := res(i) - c(j)*res(i-j);
      end loop;
      res(i) := res(i)/c(0);
    end loop;
    return res;
  end Inverse;

  function "/" ( a,b : Vector ) return Vector is

    invb : constant Vector(b'range) := Inverse(b); 

  begin
    return a*invb;
  end "/";

  function Eval ( c : Vector; t : double_float ) return Complex_Number is

    res : Complex_Number := c(0);
    pwt : double_float := t;

  begin
    for i in 1..(c'last-1) loop
      res := res + c(i)*pwt;
      pwt := pwt*t;
    end loop;
    res := res + c(c'last)*pwt;
    return res;
  end Eval;

  function Eval ( c : Vector; t : Complex_Number ) return Complex_Number is

    res : Complex_Number := c(0);
    pwt : Complex_Number := t;

  begin
    for i in 1..(c'last-1) loop
      res := res + c(i)*pwt;
      pwt := pwt*t;
    end loop;
    res := res + c(c'last)*pwt;
    return res;
  end Eval;

  function sqrt ( c : Vector ) return Vector is

  -- ALGORITHM : apply Newton's method to x^2 - c = 0,
  --   compute dx = (x^2 - c)/(2*x), followed by x = x - dx.
  --   Stop when the norm of x^2 - c is small enough.

    x : Vector(c'range) := c;
    wrk,dx : Vector(c'range);
    nrm : double_float;
    fac : constant Complex_Number := Create(0.5);
    nit : constant integer32 := 5 + c'last;

  begin
    for k in 1..nit loop
      wrk := Standard_Truncated_Series."*"(x,x) - c;
      nrm := Standard_Complex_Norms_Equals.Max_Norm(wrk);
      exit when (nrm < 1.0e-14);
      dx := fac*wrk/x;
      x := x - dx;
      nrm := Standard_Complex_Norms_Equals.Max_Norm(dx);
      exit when (nrm < 1.0e-14);
    end loop;
    return x;
  end sqrt;

end Standard_Truncated_Series;
