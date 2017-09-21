with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Random_Numbers;           use Standard_Random_Numbers;
with Standard_Complex_Numbers; 
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Numbers_Polar;
with DoblDobl_Complex_Numbers_Polar;
with QuadDobl_Complex_Numbers_Polar;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Lexicographic_Root_Enumeration;    use Lexicographic_Root_Enumeration;

package body Total_Degree_Start_Systems is

  procedure Total_Degree_Info is

  -- DESCRIPTION :
  --   Displays information about the total degree on screen.

    i : array(1..5) of string(1..65);

  begin
    i(1):="  The  total  degree  is  the  product  of  the  degrees  of  the";
    i(2):="polynomials in the system.  The i-th equation of the start system";
    i(3):="is a univariate polynomial in the i-th unknown of the same degree";
    i(4):="as  the i-th polynomial in the system that has to be solved.  The";
    i(5):="total degree equals the number of solutions of the start system. ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Total_Degree_Info;

  function Degrees ( p : Standard_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := natural32(Standard_Complex_Polynomials.Degree(p(i)));
    end loop;
    return res;
  end Degrees;

  function Degrees ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := natural32(DoblDobl_Complex_Polynomials.Degree(p(i)));
    end loop;
    return res;
  end Degrees;

  function Degrees ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
                   return Standard_Natural_Vectors.Vector is

    res : Standard_Natural_Vectors.Vector(p'range);

  begin
    for i in p'range loop
      res(i) := natural32(QuadDobl_Complex_Polynomials.Degree(p(i)));
    end loop;
    return res;
  end Degrees;

  function Product ( d : Standard_Natural_Vectors.Vector ) return natural32 is

    res : natural32 := d(d'first);

  begin
    for i in d'first+1..d'last loop
      res := res*d(i);
    end loop;
    return res;
  end Product;

  function Total_Degree
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return natural32 is
  begin
    return Product(Degrees(p));
  end Total_Degree;

  function Total_Degree
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return natural32 is
  begin
    return Product(Degrees(p));
  end Total_Degree;

  function Total_Degree
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return natural32 is
  begin
    return Product(Degrees(p));
  end Total_Degree;

-- CREATE THE SYSTEM :

  function Start_System
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    c : Standard_Complex_Vectors.Vector(p'range);

  begin
    for i in c'range loop
      c(i) := Random1;
    end loop;
    return Start_System(p,c);
  end Start_System;

  function Start_System 
             ( p : Standard_Complex_Poly_Systems.Poly_Sys;
               c : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    d : Standard_Natural_Vectors.Vector(p'range);

  begin
    for i in d'range loop
      d(i) := natural32(Standard_Complex_Polynomials.Degree(p(i)));
    end loop;
    return Start_System(d,c);
  end Start_System;

  function Start_System 
             ( d : Standard_Natural_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    c : Standard_Complex_Vectors.Vector(d'range);

  begin
    for i in c'range loop
      c(i) := Random1;
    end loop;
    return Start_System(d,c);
  end Start_System;

  function Start_System 
             ( d : Standard_Natural_Vectors.Vector;
               c : Standard_Complex_Vectors.Vector )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    res : Standard_Complex_Poly_Systems.Poly_Sys(d'range);
    n : constant integer32 := d'last;
    use Standard_Complex_Numbers,Standard_Complex_Polynomials;
    t : Term;

  begin
    t.dg := new Standard_Natural_Vectors.Vector'(1..n => 0);
    for i in d'range loop
      t.dg(i) := d(i);
      t.cf := Create(1.0);
      res(i) := Create(t);
      t.dg(i) := 0;
      t.cf := -c(i);
      Add(res(i),t);
    end loop;
    Clear(t);
    return res;
  end Start_System;

  function Coefficients
             ( q : in Standard_Complex_Polynomials.Poly )
             return Standard_Complex_Numbers.Complex_Number is

    use Standard_Complex_Numbers;
    res : Complex_Number := Create(0.0);
    iszero : boolean := true;
    use Standard_Complex_Polynomials;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if iszero
       then res := t.cf; continue := true; iszero := false;
       else res := -t.cf/res; continue := false;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(q);
    return res;
  end Coefficients;

  function Coefficients
             ( q : in DoblDobl_Complex_Polynomials.Poly )
             return DoblDobl_Complex_Numbers.Complex_Number is

    use DoblDobl_Complex_Numbers;
    res : Complex_Number := Create(integer(0));
    iszero : boolean := true;
    use DoblDobl_Complex_Polynomials;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if iszero
       then res := t.cf; continue := true; iszero := false;
       else res := -t.cf/res; continue := false;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(q);
    return res;
  end Coefficients;

  function Coefficients
             ( q : in QuadDobl_Complex_Polynomials.Poly )
             return QuadDobl_Complex_Numbers.Complex_Number is

    use QuadDobl_Complex_Numbers;
    res : Complex_Number := Create(integer(0));
    iszero : boolean := true;
    use QuadDobl_Complex_Polynomials;

    procedure Scan_Term ( t : in Term; continue : out boolean ) is
    begin
      if iszero
       then res := t.cf; continue := true; iszero := false;
       else res := -t.cf/res; continue := false;
      end if;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator(Scan_Term);

  begin
    Scan_Terms(q);
    return res;
  end Coefficients;

  function Coefficients
             ( q : in Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(q'range);

  begin
    for i in q'range loop
      res(i) := Coefficients(q(i));
    end loop;
    return res;
  end Coefficients;

  function Coefficients
             ( q : in DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(q'range);

  begin
    for i in q'range loop
      res(i) := Coefficients(q(i));
    end loop;
    return res;
  end Coefficients;

  function Coefficients
             ( q : in QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(q'range);

  begin
    for i in q'range loop
      res(i) := Coefficients(q(i));
    end loop;
    return res;
  end Coefficients;

-- PARTICULAR SOLVERS :

  function Eval ( d : in Standard_Natural_Vectors.Vector;
                  c,x : in Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector is

    use Standard_Complex_Numbers;
    res : Standard_Complex_Vectors.Vector(x'range);

  begin
    for i in res'range loop
      if d(i) = 0 then
        res(i) := Create(1.0);
      else
        res(i) := x(i);
        for j in 2..d(i) loop
          Mul(res(i),x(i));
        end loop;
      end if;
      Sub(res(i),c(i));
    end loop;
    return res;
  end Eval;

  function Root ( d : in Standard_Natural_Vectors.Vector;
                  s : in Standard_Natural_Vectors.Vector;
                  c : in Standard_Complex_Vectors.Vector )
                return Standard_Complex_Vectors.Vector is

    res : Standard_Complex_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := Standard_Complex_Numbers_Polar.Root(c(i),d(i),s(i));
    end loop;
    return res;
  end Root;

  function Root ( d : in Standard_Natural_Vectors.Vector;
                  s : in Standard_Natural_Vectors.Vector;
                  c : in DoblDobl_Complex_Vectors.Vector )
                return DoblDobl_Complex_Vectors.Vector is

    res : DoblDobl_Complex_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := DoblDobl_Complex_Numbers_Polar.Root(c(i),d(i),s(i));
    end loop;
    return res;
  end Root;

  function Root ( d : in Standard_Natural_Vectors.Vector;
                  s : in Standard_Natural_Vectors.Vector;
                  c : in QuadDobl_Complex_Vectors.Vector )
                return QuadDobl_Complex_Vectors.Vector is

    res : QuadDobl_Complex_Vectors.Vector(c'range);

  begin
    for i in c'range loop
      res(i) := QuadDobl_Complex_Numbers_Polar.Root(c(i),d(i),s(i));
    end loop;
    return res;
  end Root;

  function Create ( r : Standard_Complex_Vectors.Vector )
                  return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    s : Solution(r'last);

  begin
    s.t := Standard_Complex_Numbers.Create(0.0);
    s.m := 1; s.v := r;
    s.err := 0.0; s.rco := 1.0; s.res := 0.0;
    return s;
  end Create;

  function Create ( r : Standard_Complex_Vectors.Vector;
                    rcond : double_float )
                  return Standard_Complex_Solutions.Solution is

    use Standard_Complex_Solutions;
    s : Solution(r'last);

  begin
    s.t := Standard_Complex_Numbers.Create(0.0);
    s.m := 1; s.v := r;
    s.err := 0.0; s.rco := rcond; s.res := 0.0;
    return s;
  end Create;

  function Create ( r : DoblDobl_Complex_Vectors.Vector )
                  return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Solutions;
    s : Solution(r'last);

  begin
    s.t := DoblDobl_Complex_Numbers.Create(integer(0));
    s.m := 1; s.v := r;
    s.err := create(integer(0)); s.rco := create(integer(1)); s.res := create(integer(0));
    return s;
  end Create;

  function Create ( r : DoblDobl_Complex_Vectors.Vector;
                    rcond : double_double )
                  return DoblDobl_Complex_Solutions.Solution is

    use DoblDobl_Complex_Solutions;
    s : Solution(r'last);

  begin
    s.t := DoblDobl_Complex_Numbers.Create(integer(0));
    s.m := 1; s.v := r;
    s.err := create(integer(0)); s.rco := rcond; s.res := create(integer(0));
    return s;
  end Create;

  function Create ( r : QuadDobl_Complex_Vectors.Vector )
                  return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Solutions;
    s : Solution(r'last);

  begin
    s.t := QuadDobl_Complex_Numbers.Create(integer(0));
    s.m := 1; s.v := r;
    s.err := create(integer(0));
    s.rco := create(integer(1));
    s.res := create(integer(0));
    return s;
  end Create;

  function Create ( r : QuadDobl_Complex_Vectors.Vector;
                    rcond : quad_double )
                  return QuadDobl_Complex_Solutions.Solution is

    use QuadDobl_Complex_Solutions;
    s : Solution(r'last);

  begin
    s.t := QuadDobl_Complex_Numbers.Create(integer(0));
    s.m := 1; s.v := r;
    s.err := create(integer(0)); s.rco := rcond; s.res := create(integer(0));
    return s;
  end Create;

  function Solve ( d : in Standard_Natural_Vectors.Vector;
                   c : in Standard_Complex_Vectors.Vector )
                 return Standard_Complex_Solutions.Solution_List is

    use Standard_Complex_Solutions;
    res,res_last : Solution_List;
    n : constant natural32 := natural32(d'last);
    cont : boolean := true;
    acc : Standard_Natural_Vectors.Vector(d'range);

    procedure Compute_Root ( acc : in Standard_Natural_Vectors.Vector;
                             continue : out boolean ) is

      r : Standard_Complex_Vectors.Vector(c'range);

    begin
      r := Root(d,acc,c);
      Append(res,res_last,Create(r));
      continue := true;
    end Compute_Root;
    procedure Enum is new Lexicographic_Enumeration(Compute_Root);

  begin
    Enum(1,n,d,acc,cont);
    return res;
  end Solve;

-- COMPREHENSIVE SOLVERS :

  procedure Start_System 
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 c : in Standard_Complex_Vectors.Vector;
                 qsols : out Standard_Complex_Solutions.Solution_List ) is
  
    d : constant Standard_Natural_Vectors.Vector(p'range) := Degrees(p);

  begin
    q := Start_System(p,c);
    qsols := Solve(d,c);
  end Start_System;
 
  procedure Start_System
               ( p : in Standard_Complex_Poly_Systems.Poly_Sys;
                 q : out Standard_Complex_Poly_Systems.Poly_Sys;
                 qsols : out Standard_Complex_Solutions.Solution_List ) is

    c : Standard_Complex_Vectors.Vector(p'range);

  begin
    for i in c'range loop
      c(i) := Random1;
    end loop;
    Start_System(p,q,c,qsols);
  end Start_System;

-- RESIDUAL CALCULATION :

  function Sum_Residuals
             ( q : Standard_Complex_Poly_Systems.Poly_Sys;
               qsols : Standard_Complex_Solutions.Solution_List )
             return double_float is

    use Standard_Complex_Solutions;
    res : double_float := 0.0;
    y : Standard_Complex_Vectors.Vector(q'range);
    tmp : Solution_List := qsols;
    ls : Link_to_Solution;

  begin
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      y := Eval(q,ls.v);
      res := res + Max_Norm(y);
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Sum_Residuals;

  procedure Write_Residuals
              ( file : in file_type;
                q : in Standard_Complex_Poly_Systems.Poly_Sys;
                qsols : in Standard_Complex_Solutions.Solution_List;
                sum : out double_float ) is

    use Standard_Complex_Solutions;
    y : Standard_Complex_Vectors.Vector(q'range);
    tmp : Solution_List := qsols;
    cnt : natural32 := 0;
    nrm : double_float;
    ls : Link_to_Solution;

  begin
    sum := 0.0;
    while not Is_Null(tmp) loop
      ls := Head_Of(tmp);
      cnt := cnt + 1;
      y := Eval(q,ls.v);
      nrm := Max_Norm(y);
      put(file,cnt,4); put(file," : ");
      put(file,nrm); new_line(file);
      sum := sum + nrm;
      tmp := Tail_Of(tmp);
    end loop;
  end Write_Residuals;

end Total_Degree_Start_Systems;
