with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Complex_Numbers;
with DoblDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers;
with Standard_Complex_Norms_Equals;
with DoblDobl_Complex_Vector_Norms;
with QuadDobl_Complex_Vector_Norms;
with Standard_Complex_Singular_Values;
with DoblDobl_Complex_Singular_Values;
with QuadDobl_Complex_Singular_Values;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Complex_Poly_SysFun;
with DoblDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Poly_SysFun;
with Planes_and_Polynomials;

package body Start_Flag_Homotopies is

  function Inconsistent
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return boolean is

    use Standard_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if degree(p(i)) = 0
         then return true;
        end if;
      end if;
    end loop;
    return false;
  end Inconsistent;

  function Inconsistent
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return boolean is

    use DoblDobl_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if degree(p(i)) = 0
         then return true;
        end if;
      end if;
    end loop;
    return false;
  end Inconsistent;

  function Inconsistent
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return boolean is

    use QuadDobl_Complex_Polynomials;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if degree(p(i)) = 0
         then return true;
        end if;
      end if;
    end loop;
    return false;
  end Inconsistent;

  function Linear_Equations
             ( p : Standard_Complex_Poly_Systems.Poly_Sys )
             return Standard_Complex_Poly_Systems.Poly_Sys is

    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
    cnt : integer32 := res'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if Degree(p(i)) = 1 then
          cnt := cnt + 1;
          res(cnt) := p(i);
        end if;
      end if;
    end loop;
    return res(res'first..cnt);
  end Linear_Equations;

  function Linear_Equations
             ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys )
             return DoblDobl_Complex_Poly_Systems.Poly_Sys is

    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
    cnt : integer32 := res'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if Degree(p(i)) = 1 then
          cnt := cnt + 1;
          res(cnt) := p(i);
        end if;
      end if;
    end loop;
    return res(res'first..cnt);
  end Linear_Equations;

  function Linear_Equations
             ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys )
             return QuadDobl_Complex_Poly_Systems.Poly_Sys is

    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Systems;

    res : Poly_Sys(p'range);
    cnt : integer32 := res'first-1;

  begin
    for i in p'range loop
      if p(i) /= Null_Poly then
        if Degree(p(i)) = 1 then
          cnt := cnt + 1;
          res(cnt) := p(i);
        end if;
      end if;
    end loop;
    return res(res'first..cnt);
  end Linear_Equations;

  procedure Coefficients
             ( p : in Standard_Complex_Poly_Systems.Poly_Sys; 
               A : out Standard_Complex_Matrices.Matrix;
               b : out Standard_Complex_Vectors.Vector ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    h : Standard_Complex_Vectors.Vector(0..integer32(n));

  begin
    for i in p'range loop
      h := Planes_and_Polynomials.Polynomial(p(i));
      b(i) := -h(0);
      for j in A'range(2) loop
        A(i,j) := h(j);
      end loop;
    end loop;
  end Coefficients;

  procedure Coefficients
             ( p : in DoblDobl_Complex_Poly_Systems.Poly_Sys; 
               A : out DoblDobl_Complex_Matrices.Matrix;
               b : out DoblDobl_Complex_Vectors.Vector ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    h : DoblDobl_Complex_Vectors.Vector(0..integer32(n));

  begin
    for i in p'range loop
      h := Planes_and_Polynomials.Polynomial(p(i));
      b(i) := -h(0);
      for j in A'range(2) loop
        A(i,j) := h(j);
      end loop;
    end loop;
  end Coefficients;

  procedure Coefficients
             ( p : in QuadDobl_Complex_Poly_Systems.Poly_Sys; 
               A : out QuadDobl_Complex_Matrices.Matrix;
               b : out QuadDobl_Complex_Vectors.Vector ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;

    n : constant natural32 := Number_of_Unknowns(p(p'first));
    h : QuadDobl_Complex_Vectors.Vector(0..integer32(n));

  begin
    for i in p'range loop
      h := Planes_and_Polynomials.Polynomial(p(i));
      b(i) := -h(0);
      for j in A'range(2) loop
        A(i,j) := h(j);
      end loop;
    end loop;
  end Coefficients;

  procedure Solve ( A : in out Standard_Complex_Matrices.Matrix;
                    b : in Standard_Complex_Vectors.Vector;
                    x : out Standard_Complex_Vectors.Vector;
                    res : out double_float ) is

    use Standard_Complex_Singular_Values;

    AA : constant Standard_Complex_Matrices.Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    m : constant integer32 := Min0(n+1,p);
    s : Standard_Complex_Vectors.Vector(1..m);
    e : Standard_Complex_Vectors.Vector(1..p);
    u : Standard_Complex_Matrices.Matrix(1..n,1..n);
    v : Standard_Complex_Matrices.Matrix(1..p,1..p);
    info : integer32;
    r : Standard_Complex_Vectors.Vector(b'range);
    use Standard_Complex_Vectors;
    use Standard_Complex_Matrices;

  begin
    SVD(A,n,p,s,e,u,v,11,info);
    x := Solve(u,v,s,b);
    r := b - AA*x;
   -- put_line("the singular values : "); put_line(s);
   -- put_line("The residual :"); put_line(r);
    res := Standard_Complex_Norms_Equals.Norm2(r);
  end Solve;

  procedure Solve ( A : in out DoblDobl_Complex_Matrices.Matrix;
                    b : in DoblDobl_Complex_Vectors.Vector;
                    x : out DoblDobl_Complex_Vectors.Vector;
                    res : out double_double ) is


    use DoblDobl_Complex_Singular_Values;

    AA : constant DoblDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    m : constant integer32 := Min0(n+1,p);
    s : DoblDobl_Complex_Vectors.Vector(1..m);
    e : DoblDobl_Complex_Vectors.Vector(1..p);
    u : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : DoblDobl_Complex_Matrices.Matrix(1..p,1..p);
    info : integer32;
    r : DoblDobl_Complex_Vectors.Vector(b'range);
    use DoblDobl_Complex_Vectors;
    use DoblDobl_Complex_Matrices;

  begin
    SVD(A,n,p,s,e,u,v,11,info);
    x := Solve(u,v,s,b);
    r := b - AA*x;
   -- put_line("the singular values : "); put_line(s);
   -- put_line("The residual :"); put_line(r);
    res := DoblDobl_Complex_Vector_Norms.Norm2(r);
  end Solve;

  procedure Solve ( A : in out QuadDobl_Complex_Matrices.Matrix;
                    b : in QuadDobl_Complex_Vectors.Vector;
                    x : out QuadDobl_Complex_Vectors.Vector;
                    res : out quad_double ) is

    use QuadDobl_Complex_Singular_Values;

    AA : constant QuadDobl_Complex_Matrices.Matrix(A'range(1),A'range(2)) := A;
    n : constant integer32 := A'last(1);
    p : constant integer32 := A'last(2);
    m : constant integer32 := Min0(n+1,p);
    s : QuadDobl_Complex_Vectors.Vector(1..m);
    e : QuadDobl_Complex_Vectors.Vector(1..p);
    u : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    v : QuadDobl_Complex_Matrices.Matrix(1..p,1..p);
    info : integer32;
    r : QuadDobl_Complex_Vectors.Vector(b'range);
    use QuadDobl_Complex_Vectors;
    use QuadDobl_Complex_Matrices;

  begin
    SVD(A,n,p,s,e,u,v,11,info);
    x := Solve(u,v,s,b);
    r := b - AA*x;
   -- put_line("the singular values : "); put_line(s);
   -- put_line("The residual :"); put_line(r);
    res := QuadDobl_Complex_Vector_Norms.Norm2(r);
  end Solve;

  procedure First_Solution
             ( f : in Standard_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out Standard_Complex_Vectors.Vector;
               res : out double_float ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Systems;

    n : constant integer32 := x'last;

  begin
    fail := true;
    x := (x'range => Create(0.0));
    if not Inconsistent(f) then
      declare
        s : constant Poly_Sys := Linear_Equations(f);
        m : constant integer32 := s'last;
        A : Standard_Complex_Matrices.Matrix(1..m,1..n);
        b : Standard_Complex_Vectors.Vector(1..m);
      begin
        if s'last >= n then
          Coefficients(s,A,b);
         -- put_line("The coefficient matrix : "); put(A,3);
         -- put_line("The righthandside vector : "); put_line(b);
          Solve(A,b,x,res);
          fail := (res > 1.0E-8);
        end if;
      end;
    end if;
  end First_Solution;

  procedure First_Solution
             ( f : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out DoblDobl_Complex_Vectors.Vector;
               res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Systems;

    n : constant integer32 := x'last;
    zero : constant double_double := create(0.0);

  begin
    fail := true;
    x := (x'range => Create(zero));
    if not Inconsistent(f) then
      declare
        s : constant Poly_Sys := Linear_Equations(f);
        m : constant integer32 := s'last;
        A : DoblDobl_Complex_Matrices.Matrix(1..m,1..n);
        b : DoblDobl_Complex_Vectors.Vector(1..m);
      begin
        if s'last >= n then
          Coefficients(s,A,b);
         -- put_line("The coefficient matrix : "); put(A,3);
         -- put_line("The righthandside vector : "); put_line(b);
          Solve(A,b,x,res);
          fail := (res > 1.0E-8);
        end if;
      end;
    end if;
  end First_Solution;

  procedure First_Solution
             ( f : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
               fail : out boolean;
               x : out QuadDobl_Complex_Vectors.Vector;
               res : out quad_double ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Systems;

    n : constant integer32 := x'last;
    zero : constant quad_double := create(0.0);

  begin
    fail := true;
    x := (x'range => Create(zero));
    if not Inconsistent(f) then
      declare
        s : constant Poly_Sys := Linear_Equations(f);
        m : constant integer32 := s'last;
        A : QuadDobl_Complex_Matrices.Matrix(1..m,1..n);
        b : QuadDobl_Complex_Vectors.Vector(1..m);
      begin
        if s'last >= n then
          Coefficients(s,A,b);
         -- put_line("The coefficient matrix : "); put(A,3);
         -- put_line("The righthandside vector : "); put_line(b);
          Solve(A,b,x,res);
          fail := (res > 1.0E-8);
        end if;
      end;
    end if;
  end First_Solution;

  procedure Start_Solution
              ( h : in Standard_Complex_Poly_Systems.Poly_Sys;
                fail : out boolean;
                x : out Standard_Complex_Vectors.Vector;
                res : out double_float ) is

    use Standard_Complex_Numbers;
    use Standard_Complex_Poly_Systems;
    use Standard_Complex_Poly_SysFun;

    n : constant integer32 := x'last;
    h0 : Poly_Sys(h'range) := Eval(h,Create(0.0),n+1);
      -- h0 is system h where t is substituted by zero

  begin
    First_Solution(h0,fail,x,res);
    Clear(h0);
  end Start_Solution;

  procedure Start_Solution
              ( h : in DoblDobl_Complex_Poly_Systems.Poly_Sys;
                fail : out boolean;
                x : out DoblDobl_Complex_Vectors.Vector;
                res : out double_double ) is

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Poly_Systems;
    use DoblDobl_Complex_Poly_SysFun;

    n : constant integer32 := x'last;
    zero : constant double_double := create(0.0);
    h0 : Poly_Sys(h'range) := Eval(h,Create(zero),n+1);
      -- h0 is system h where t is substituted by zero

  begin
    First_Solution(h0,fail,x,res);
    Clear(h0);
  end Start_Solution;

  procedure Start_Solution
              ( h : in QuadDobl_Complex_Poly_Systems.Poly_Sys;
                fail : out boolean;
                x : out QuadDobl_Complex_Vectors.Vector;
                res : out quad_double ) is

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Poly_Systems;
    use QuadDobl_Complex_Poly_SysFun;

    n : constant integer32 := x'last;
    zero : constant quad_double := create(0.0);
    h0 : Poly_Sys(h'range) := Eval(h,Create(zero),n+1);
      -- h0 is system h where t is substituted by zero

  begin
    First_Solution(h0,fail,x,res);
    Clear(h0);
  end Start_Solution;

end Start_Flag_Homotopies;
