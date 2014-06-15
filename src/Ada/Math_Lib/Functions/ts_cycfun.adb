with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Standard_Complex_Numbers;           use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural_VecVecs_io;        use Standard_Natural_VecVecs_io;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Complex_Polynomials;       use Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;    use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;      use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;       use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;     use Standard_Complex_Jaco_Matrices;
with Standard_Speelpenning_Products;     use Standard_Speelpenning_Products;

procedure ts_cycfun is

-- DESCRIPTION :
--   The goal of this test program is to test the performance of the
--   evaluation and differentiation schemes in reverse mode on the
--   benchmark problem of cyclic n-roots.
--   Run the performance test on n = 10 and m = 50000 to see the improvement.

  function Number_of_Monomials ( n,i : integer32 ) return integer32 is

  -- DESCRIPTION :
  --   Returns the number of monomials in the i-th polynomial of
  --   the cyclic n-roots problem.

  begin
    if i = n
     then return 2;
     else return n;
    end if;
  end Number_of_Monomials;

  function Support_of_Cyclic ( n,i : integer32 )
             return Standard_Natural_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the support of the i-th polynomial in the
  --   cyclic n-roots system.

  -- REQUIRED : 0 < i < n

    dim : constant integer32 := Number_of_Monomials(n,i);
    res : Standard_Natural_VecVecs.VecVec(1..dim);
    deg : Standard_Natural_Vectors.Vector(1..n) := (1..n => 0);
    tmp : natural32;

  begin
    for k in 1..i loop
      deg(k) := 1;
    end loop;
    res(1) := new Standard_Natural_Vectors.Vector'(deg);
    if dim = 2 then
      deg := (1..n => 0);
      res(2) := new Standard_Natural_Vectors.Vector'(deg);
    else
      for k in 2..res'last loop
        tmp := deg(n);
        for j in reverse 1..n-1 loop
          deg(j+1) := deg(j);
        end loop;
        deg(1) := tmp;
        res(k) := new Standard_Natural_Vectors.Vector'(deg);
      end loop;
    end if;
    return res;
  end Support_of_Cyclic;

  function Supports_of_Cyclic ( n : integer32 )
             return Standard_Natural_VecVecs.Array_of_VecVecs is

  -- DESCRIPTION :
  --   Returns the supports of the cyclic n-roots problem.

    res : Standard_Natural_VecVecs.Array_of_VecVecs(1..n);

  begin
    for i in 1..n loop
      declare
        s : constant Standard_Natural_VecVecs.VecVec
          := Support_of_Cyclic(n,i);
      begin
        res(i) := new Standard_Natural_VecVecs.VecVec'(s);
      end;
    end loop;
    return res;
  end Supports_of_Cyclic;

  function Cyclic_Polynomial
             ( s : Standard_Natural_VecVecs.VecVec ) return Poly is

  -- DESCRIPTION :
  --   Returns the cyclic polynomial defined by the support in s.

    res : Poly := Null_Poly;
    t : Term;

  begin
    t.cf := Create(1.0);
    t.dg := new Standard_Natural_Vectors.Vector'(s(s'first).all);
    Add(res,t);
    if s'last = s'first+1 then
      t.dg.all := s(s'last).all;
      Sub(res,t);
    else
      for i in s'first+1..s'last loop
        t.dg.all := s(i).all;
        Add(res,t);
      end loop;
    end if;
    Clear(t);
    return res;
  end Cyclic_Polynomial;

  function Cyclic_Polynomial_System
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Poly_Sys is

  -- DESCRIPTION :
  --   Returns the polynomial system defined by the supports in s.

    res : Poly_Sys(s'range);

  begin
    for i in s'range loop
      res(i) := Cyclic_Polynomial(s(i).all);
    end loop;
    return res;
  end Cyclic_Polynomial_System;

  function Indexed_Support
             ( s : Standard_Natural_VecVecs.VecVec )
             return Standard_Integer_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns the indexed supports of s.

    res : Standard_Integer_VecVecs.VecVec(s'range);

  begin
    for i in s'range loop
      declare
        idx : constant Standard_Integer_Vectors.Vector
            := Standard_Speelpenning_Products.Nonzero_Indices(s(i).all);
      begin
        res(i) := new Standard_Integer_Vectors.Vector'(idx);
      end;
    end loop;
    return res;
  end Indexed_Support;

  function Indexed_Supports
             ( s : Standard_Natural_VecVecs.Array_of_VecVecs )
             return Standard_Integer_VecVecs.Array_of_VecVecs is

  -- DESCRIPTION :
  --   Returns the indexed supports of s.

    res : Standard_Integer_VecVecs.Array_of_VecVecs(s'range);

  begin
    for i in s'range loop
      declare
        idx : constant Standard_Integer_VecVecs.VecVec
            := Indexed_Support(s(i).all);
      begin
        res(i) := new Standard_Integer_VecVecs.VecVec'(idx);
      end;
    end loop;
    return res;
  end Indexed_Supports;

  function Gradient_of_Cyclic
             ( s : Standard_Natural_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector ) 
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its support in s at x.

    res,y : Standard_Complex_Vectors.Vector(0..x'last);

  begin
    res := Reverse_Speel(s(s'first).all,x);
    if s'first+1 = s'last then
      y := Reverse_Speel(s(s'last).all,x);
      Standard_Complex_Vectors.Sub(res,y);
    else
      for i in s'first+1..s'last loop
        y := Reverse_Speel(s(i).all,x);
        Standard_Complex_Vectors.Add(res,y);
      end loop;
    end if;
    return res;
  end Gradient_of_Cyclic;

  function Indexed_Gradient_of_Cyclic
             ( s : Standard_Integer_VecVecs.VecVec;
               x : Standard_Complex_Vectors.Vector ) 
             return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its indexed support in s at x.

    res,y : Standard_Complex_Vectors.Vector(0..x'last);

  begin
    res := Indexed_Reverse_Speel(s(s'first).all,x);
    if s'first+1 = s'last then
      y := Indexed_Reverse_Speel(s(s'last).all,x);
      Standard_Complex_Vectors.Sub(res,y);
    else
      for i in s'first+1..s'last loop
        y := Indexed_Reverse_Speel(s(i).all,x);
        Standard_Complex_Vectors.Add(res,y);
      end loop;
    end if;
    return res;
  end Indexed_Gradient_of_Cyclic;

  procedure Evaluate_and_Differentiate
              ( s : in Standard_Natural_VecVecs.Array_of_VecVecs;
                x : in Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x.
  --   The function evaluations are in y and the Jacobian matrix in A.

    z : Standard_Complex_Vectors.Vector(0..x'last);

  begin
    for i in s'range loop
      z := Gradient_of_Cyclic(s(i).all,x);
      y(i) := z(0);
      for j in 1..z'last loop
        A(i,j) := z(j);
      end loop;
    end loop;
  end Evaluate_and_Differentiate;

  procedure Indexed_Evaluate_and_Differentiate
              ( s : in Standard_Integer_VecVecs.Array_of_VecVecs;
                x : in Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s.
  --   The function evaluations are in y and the Jacobian matrix in A.

    z : Standard_Complex_Vectors.Vector(0..x'last);

  begin
    for i in s'range loop
      z := Indexed_Gradient_of_Cyclic(s(i).all,x);
      y(i) := z(0);
      for j in 1..z'last loop
        A(i,j) := z(j);
      end loop;
    end loop;
  end Indexed_Evaluate_and_Differentiate;

  procedure Evaluate_and_Differentiate
              ( f : in Eval_Poly_Sys; jf : in Eval_Jaco_Mat;
                x : in Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x.
  --   The function evaluations are in y and the Jacobian matrix in A.

  begin
    y := Eval(f,x);
    A := Eval(jf,x);
  end Evaluate_and_Differentiate;

  procedure Random_Point_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point.

    x : constant Standard_Complex_Vectors.Vector
      := Standard_Random_Vectors.Random_Vector(1,n);
    y,y2 : Standard_Complex_Vectors.Vector(0..n);
    p,q : Poly;
    z : Complex_Number;

  begin
    for i in 1..n loop
      declare
        s : constant Standard_Natural_VecVecs.VecVec
          := Support_of_Cyclic(n,i);    
        idx : constant Standard_Integer_VecVecs.VecVec := Indexed_Support(s);
      begin
        put("Support "); put(i,1); put(" of cyclic "); put(n,1);
        put_line("-roots :");
        put(s);
        p := Cyclic_Polynomial(s);
        put_line("The polynomial : "); put(p); new_line;
        y := Gradient_of_Cyclic(s,x);
        put_line("Its gradient at a random point : "); put_line(y);
        y2 := Indexed_Gradient_of_Cyclic(idx,x);
        put_line("Computed via its indexing : "); put_line(y2);
        z := Eval(p,x);
        put_line("The function value and its derivatives at the point :");
        put(z); new_line;
        for i in 1..n loop
          q := Diff(p,i);
          z := Eval(q,x);
          put(z); new_line;
          Clear(q);
        end loop;
        Clear(p);
      end;
    end loop;
  end Random_Point_Test;

  procedure Performance_Test ( n : integer32 ) is

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem.

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    idx : Standard_Integer_VecVecs.Array_of_VecVecs(1..n)
        := Indexed_Supports(s);
    p : Poly_Sys(1..n) := Cyclic_Polynomial_System(s);
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    x : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    y : Standard_Complex_Vectors.Vector(1..n);
    A : Standard_Complex_Matrices.Matrix(1..n,1..n);
    timer : Timing_Widget;
    m : integer32 := 0;

  begin
    new_line;
    put("Give the number of evaluations and differentiations : ");
    get(m);
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      Evaluate_and_Differentiate(s,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"in reverse mode");
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      Indexed_Evaluate_and_Differentiate(idx,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"with precomputed indexes");
    tstart(timer);
    for i in 1..m loop
      x := Standard_Random_Vectors.Random_Vector(1,n);
      Evaluate_and_Differentiate(f,jf,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"with Horner schemes");
  end Performance_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension and generates then the
  --   support for the cyclic n-roots problem.

    n : integer32 := 0;
    ans : character;

  begin
    put("Give the dimension : "); get(n);
    declare
      s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
        := Supports_of_Cyclic(n);
      p : Poly_Sys(1..n) := Cyclic_Polynomial_System(s);
    begin
      put("The cyclic "); put(n,1); put_line("-roots system :"); put(p);
      put_line("The supports : ");
      for i in s'range loop
        put(s(i).all); new_line;
      end loop;
    end;
    put("Do a random point test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Random_Point_Test(n);
    end if;
    put("Do a performance test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Performance_Test(n);
    end if;
  end Main;

begin
  Main;
end ts_cycfun;
