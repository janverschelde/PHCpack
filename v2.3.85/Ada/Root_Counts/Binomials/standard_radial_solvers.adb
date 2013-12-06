with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Complex_Numbers_Polar;    use Standard_Complex_Numbers_Polar;
with Multprec_Integer_Numbers;          use Multprec_Integer_Numbers;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Vectors;         use Multprec_Floating_Vectors;

-- for testing:
with text_io;                           use text_io;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;     

package body Standard_Radial_Solvers is

  function Radii ( c : Standard_Complex_Vectors.Vector )
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := Radius(c(i));
    end loop;
    return res;
  end Radii;

  function Scale ( c : Standard_Complex_Vectors.Vector;
                   r : Standard_Floating_Vectors.Vector )
                 return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := c(i)/r(i);
    end loop;
    return res;
  end Scale;

  function Log10 ( r : Standard_Floating_Vectors.Vector ) 
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(r'range);

  begin
    for i in r'range loop
      res(i) := LOG10(r(i));
    end loop;
    return res;
  end Log10;

  function Exp10 ( r : Standard_Floating_Vectors.Vector ) 
                 return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(r'range);

  begin
    for i in r'range loop
      res(i) := 10.0**r(i);
    end loop;
    return res;
  end Exp10;

  function Radial_Upper_Solve
              ( U : Standard_Integer_Matrices.Matrix;
                logr : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(logr'range)
        := (logr'range => 0.0);
    eva : double_float;

  begin
    for j in U'range loop
      eva := logr(j);
      for i in 1..j-1 loop
        eva := eva - double_float(U(i,j))*res(i);
      end loop;
      res(j) := eva/double_float(U(j,j));
    end loop;
    return res;
  end Radial_Upper_Solve;

  function to_double_float ( x : Integer_Number ) return double_float is

  -- DESCRIPTION :
  --   Converts the multiprecision integer number to a double float.

    f : Floating_Number := create(x);
    res : constant double_float := Round(f);

  begin
    Clear(f);
    return res;
  end to_double_float;

  function Decimal_Places
              ( U : Multprec_Integer_Matrices.Matrix )
              return natural32 is

  -- DESCRIPTION :
  --   Returns the largest number of decimal places of all numbers in U.

    res : natural32 := 0;
    dp : natural32;

  begin
    for i in U'range(1) loop
      for j in U'range(2) loop
        dp := Decimal_Places(U(i,j));
        if dp > res
         then res := dp;
        end if;
      end loop;
    end loop;
    return res;
  end Decimal_Places;

  function Radial_Upper_Solve
              ( U : Multprec_Integer_Matrices.Matrix;
                logr : Standard_Floating_Vectors.Vector )
              return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(logr'range)
        := (logr'range => 0.0);
    eva : double_float;
    dp : constant natural32 := Decimal_Places(U);

  begin
    put("decimal places of U is "); put(dp,1); new_line;
    for j in U'range loop
      eva := logr(j);
      for i in 1..j-1 loop
        eva := eva - to_double_float(U(i,j))*res(i);
      end loop;
      res(j) := eva/to_double_float(U(j,j));
    end loop;
    return res;
  end Radial_Upper_Solve;

  function Multiply ( A : in Standard_Integer_Matrices.Matrix;
                      x : in Standard_Floating_Vectors.Vector )
                    return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(x'range) := (x'range => 0.0);

  begin
    for j in A'range(1) loop
      for i in A'range(2) loop
        res(j) := res(j) + double_float(A(i,j))*x(i);
      end loop;
    end loop;
    return res;
  end Multiply;

  function Multiply ( A : in Multprec_Integer_Matrices.Matrix;
                      x : in Standard_Floating_Vectors.Vector )
                    return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(x'range) := (x'range => 0.0);
    mp_res : Multprec_Floating_Vectors.Vector(x'range);
    mp_x : Multprec_Floating_Vectors.Vector(x'range);

  begin
    for i in x'range loop
      mp_x(i) := Create(x(i));
      mp_res(i) := Create(0.0);
    end loop;
    for j in A'range(1) loop
      for i in A'range(2) loop
       -- res(j) := res(j) + to_double_float(A(i,j))*x(i);
        declare
          mpAij : Floating_Number := Create(A(i,j));
          acc : Floating_Number :=  mpAij*mp_x(i);
        begin
          Add(mp_res(j),acc);
          Clear(mpAij); Clear(acc);
        end;
      end loop;
    end loop;
    for i in res'range loop
      res(i) := Round(mp_res(i));
      Clear(mp_res(i)); Clear(mp_x(i));
    end loop;
    return res;
  end Multiply;

  function Eval ( A : in Standard_Integer_Matrices.Matrix;
                  x : in Standard_Floating_Vectors.Vector )
                return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(x'range) := (x'range => 1.0);

  begin
    for j in A'range(2) loop
      for i in A'range(1) loop
        res(j) := res(j)*(x(i)**integer(A(i,j)));
      end loop;
    end loop;
    return res;
  end Eval;

  function Binary_Exponentiation
             ( x : double_float; e : Integer_Number ) return double_float is

  -- DESCRIPTION :
  --   Returns x^e using binary exponentiation.

    res,acc : double_float;
    abse : Integer_Number;

  begin
    if Equal(e,0) then
      res := 1.0;
    else
      if e > 0
       then Copy(e,abse);
       else abse := -e;
      end if;
      res := x; acc := 1.0;
      if abse > 1 then          -- use binary exponentiation
        while abse > 0 loop
          if Rmd(abse,2) = 1
           then acc := acc*res;
          end if;
          Div(abse,2);
          if abse > 0
           then res := res*res;
          end if;
        end loop;
      else
        acc := res;
      end if;
      if e < 0
       then res := 1.0/acc;          -- compute reciprocal
       else res := acc;
      end if;
    end if;
    return res;
  end Binary_Exponentiation;

  function Eval ( A : in Multprec_Integer_Matrices.Matrix;
                  x : in Standard_Floating_Vectors.Vector )
                return Standard_Floating_Vectors.Vector is

    res : Standard_Floating_Vectors.Vector(x'range) := (x'range => 1.0);
    pwr : double_float;

  begin
    for j in A'range(2) loop
      for i in A'range(1) loop
        pwr := Binary_Exponentiation(x(i),A(i,j));
        res(j) := res(j)*pwr; -- res(j) := res(j)*(x(i)**A(i,j));
      end loop;
    end loop;
    return res;
  end Eval;

  procedure Multiply ( s : in out Standard_Complex_Vectors.Vector;
                       r : in Standard_Floating_Vectors.Vector ) is
  begin
    for i in s'range loop
      s(i) := s(i)*r(i);
    end loop;
  end Multiply;

  procedure Multiply ( s : in out Solution;
                       r : in Standard_Floating_Vectors.Vector )  is
  begin
    Multiply(s.v,r);
  end Multiply;

  procedure Multiply ( s : in out Solution_List;
                       r : in Standard_Floating_Vectors.Vector ) is

    tmp : Solution_List := s;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      Multiply(ls.v,r);
      Set_Head(tmp,ls);
      tmp := Tail_Of(tmp);
    end loop;
  end Multiply;

end Standard_Radial_Solvers;
