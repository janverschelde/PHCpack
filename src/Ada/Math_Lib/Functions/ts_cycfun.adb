with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Standard_Natural_Numbers;           use Standard_Natural_Numbers;
with Standard_Integer_Numbers;           use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;        use Standard_Integer_Numbers_io;
with Double_Double_Numbers;              use Double_Double_Numbers;
with Quad_Double_Numbers;                use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;        use Standard_Complex_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;        use DoblDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;        use QuadDobl_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_VecVecs;
with Standard_Natural_VecVecs_io;        use Standard_Natural_VecVecs_io;
with Standard_Integer_Vectors;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;        use Standard_Complex_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;        use DoblDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;        use QuadDobl_Complex_Vectors_io;
with Standard_Random_Vectors;
with DoblDobl_Random_Vectors;
with QuadDobl_Random_Vectors;
with Standard_Complex_Matrices;
with DoblDobl_Complex_Matrices;
with QuadDobl_Complex_Matrices;
with Standard_Complex_Polynomials;
with DoblDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials;
with Standard_Complex_Polynomials_io;    use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials_io;    use DoblDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Polynomials_io;    use QuadDobl_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;   use Standard_Complex_Poly_Systems_io;
with DoblDobl_Complex_Poly_Systems_io;   use DoblDobl_Complex_Poly_Systems_io;
with QuadDobl_Complex_Poly_Systems_io;   use QuadDobl_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Standard_Speelpenning_Products;     use Standard_Speelpenning_Products;
with DoblDobl_Speelpenning_Products;     use DoblDobl_Speelpenning_Products;
with QuadDobl_Speelpenning_Products;     use QuadDobl_Speelpenning_Products;
with Cyclic_Roots_System;                use Cyclic_Roots_System;

procedure ts_cycfun is

-- DESCRIPTION :
--   The goal of this test program is to test the performance of the
--   evaluation and differentiation schemes in reverse mode on the
--   benchmark problem of cyclic n-roots.
--   Run the performance test on n = 10 and m = 50000 to see the improvement.

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
  --   defined by its support in s at x, in standard double precision.

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

  function Gradient_of_Cyclic
             ( s : Standard_Natural_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector ) 
             return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its support in s at x, in double double precision.

    res,y : DoblDobl_Complex_Vectors.Vector(0..x'last);

  begin
    res := Reverse_Speel(s(s'first).all,x);
    if s'first+1 = s'last then
      y := Reverse_Speel(s(s'last).all,x);
      DoblDobl_Complex_Vectors.Sub(res,y);
    else
      for i in s'first+1..s'last loop
        y := Reverse_Speel(s(i).all,x);
        DoblDobl_Complex_Vectors.Add(res,y);
      end loop;
    end if;
    return res;
  end Gradient_of_Cyclic;

  function Gradient_of_Cyclic
             ( s : Standard_Natural_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector ) 
             return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its support in s at x, in quad double precision.

    res,y : QuadDobl_Complex_Vectors.Vector(0..x'last);

  begin
    res := Reverse_Speel(s(s'first).all,x);
    if s'first+1 = s'last then
      y := Reverse_Speel(s(s'last).all,x);
      QuadDobl_Complex_Vectors.Sub(res,y);
    else
      for i in s'first+1..s'last loop
        y := Reverse_Speel(s(i).all,x);
        QuadDobl_Complex_Vectors.Add(res,y);
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
  --   defined by its indexed support in s at x, in standard double precision.

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

  function Indexed_Gradient_of_Cyclic
             ( s : Standard_Integer_VecVecs.VecVec;
               x : DoblDobl_Complex_Vectors.Vector ) 
             return DoblDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its indexed support in s at x, in double double precision.

    res,y : DoblDobl_Complex_Vectors.Vector(0..x'last);

  begin
    res := Indexed_Reverse_Speel(s(s'first).all,x);
    if s'first+1 = s'last then
      y := Indexed_Reverse_Speel(s(s'last).all,x);
      DoblDobl_Complex_Vectors.Sub(res,y);
    else
      for i in s'first+1..s'last loop
        y := Indexed_Reverse_Speel(s(i).all,x);
        DoblDobl_Complex_Vectors.Add(res,y);
      end loop;
    end if;
    return res;
  end Indexed_Gradient_of_Cyclic;

  function Indexed_Gradient_of_Cyclic
             ( s : Standard_Integer_VecVecs.VecVec;
               x : QuadDobl_Complex_Vectors.Vector ) 
             return QuadDobl_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the gradient vector of the cyclic n-roots polynomial
  --   defined by its indexed support in s at x, in quad double precision.

    res,y : QuadDobl_Complex_Vectors.Vector(0..x'last);

  begin
    res := Indexed_Reverse_Speel(s(s'first).all,x);
    if s'first+1 = s'last then
      y := Indexed_Reverse_Speel(s(s'last).all,x);
      QuadDobl_Complex_Vectors.Sub(res,y);
    else
      for i in s'first+1..s'last loop
        y := Indexed_Reverse_Speel(s(i).all,x);
        QuadDobl_Complex_Vectors.Add(res,y);
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
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in standard double precision.
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

  procedure Evaluate_and_Differentiate
              ( s : in Standard_Natural_VecVecs.Array_of_VecVecs;
                x : in DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in double double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    z : DoblDobl_Complex_Vectors.Vector(0..x'last);

  begin
    for i in s'range loop
      z := Gradient_of_Cyclic(s(i).all,x);
      y(i) := z(0);
      for j in 1..z'last loop
        A(i,j) := z(j);
      end loop;
    end loop;
  end Evaluate_and_Differentiate;

  procedure Evaluate_and_Differentiate
              ( s : in Standard_Natural_VecVecs.Array_of_VecVecs;
                x : in QuadDobl_Complex_Vectors.Vector;
                y : out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in quad double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    z : QuadDobl_Complex_Vectors.Vector(0..x'last);

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
                z : in out Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s, in standard double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

  -- REQUIRED :
  --   The vector z is used as work space and is of range 0..x'last.

  -- z : Standard_Complex_Vectors.Vector(0..x'last);

  begin
    for i in s'range loop
      z := Indexed_Gradient_of_Cyclic(s(i).all,x);
      y(i) := z(0);
      for j in 1..z'last loop
        A(i,j) := z(j);
      end loop;
    end loop;
  end Indexed_Evaluate_and_Differentiate;

  procedure Indexed_Evaluate_and_Differentiate
              ( s : in Standard_Integer_VecVecs.Array_of_VecVecs;
                x : in DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s, in double double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    z : DoblDobl_Complex_Vectors.Vector(0..x'last);

  begin
    for i in s'range loop
      z := Indexed_Gradient_of_Cyclic(s(i).all,x);
      y(i) := z(0);
      for j in 1..z'last loop
        A(i,j) := z(j);
      end loop;
    end loop;
  end Indexed_Evaluate_and_Differentiate;

  procedure Indexed_Evaluate_and_Differentiate
              ( s : in Standard_Integer_VecVecs.Array_of_VecVecs;
                x : in QuadDobl_Complex_Vectors.Vector;
                y : out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   using the indexed supports in s, in quad double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    z : QuadDobl_Complex_Vectors.Vector(0..x'last);

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
              ( f : in Standard_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in Standard_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in Standard_Complex_Vectors.Vector;
                y : out Standard_Complex_Vectors.Vector;
                A : out Standard_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in standard double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;

  begin
    y := Eval(f,x);
    A := Eval(jf,x);
  end Evaluate_and_Differentiate;

  procedure Evaluate_and_Differentiate
              ( f : in DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in DoblDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in DoblDobl_Complex_Vectors.Vector;
                y : out DoblDobl_Complex_Vectors.Vector;
                A : out DoblDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in double double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

  begin
    y := Eval(f,x);
    A := Eval(jf,x);
  end Evaluate_and_Differentiate;

  procedure Evaluate_and_Differentiate
              ( f : in QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys;
                jf : in QuadDobl_Complex_Jaco_Matrices.Eval_Jaco_Mat;
                x : in QuadDobl_Complex_Vectors.Vector;
                y : out QuadDobl_Complex_Vectors.Vector;
                A : out QuadDobl_Complex_Matrices.Matrix ) is

  -- DESCRIPTION :
  --   Evaluates and differentiates the cyclic n-roots system at x,
  --   in quad double precision.
  --   The function evaluations are in y and the Jacobian matrix in A.

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

  begin
    y := Eval(f,x);
    A := Eval(jf,x);
  end Evaluate_and_Differentiate;

  procedure Standard_Random_Point_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point,
  --   in standard double precision.

    use Standard_Complex_Numbers;
    use Standard_Complex_Polynomials;
    use Standard_Complex_Poly_Functions;

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
        p := Standard_Cyclic_Polynomial(s);
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
  end Standard_Random_Point_Test;

  procedure DoblDobl_Random_Point_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point,
  --   in double double precision.

    use DoblDobl_Complex_Numbers;
    use DoblDobl_Complex_Polynomials;
    use DoblDobl_Complex_Poly_Functions;

    x : constant DoblDobl_Complex_Vectors.Vector
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    y,y2 : DoblDobl_Complex_Vectors.Vector(0..n);
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
        p := DoblDobl_Cyclic_Polynomial(s);
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
  end DoblDobl_Random_Point_Test;

  procedure QuadDobl_Random_Point_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Tests the evaluation and differentiation of the polynomials
  --   in the cyclic n-roots problem at a random point,
  --   in double double precision.

    use QuadDobl_Complex_Numbers;
    use QuadDobl_Complex_Polynomials;
    use QuadDobl_Complex_Poly_Functions;

    x : constant QuadDobl_Complex_Vectors.Vector
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    y,y2 : QuadDobl_Complex_Vectors.Vector(0..n);
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
        p := QuadDobl_Cyclic_Polynomial(s);
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
  end QuadDobl_Random_Point_Test;

  procedure Standard_Performance_Test ( n : integer32 ) is

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem,
  --   in standard double precision.

    use Standard_Complex_Poly_SysFun;
    use Standard_Complex_Jaco_Matrices;

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    idx : Standard_Integer_VecVecs.Array_of_VecVecs(1..n)
        := Indexed_Supports(s);
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..n)
      := Standard_Cyclic_Polynomial_System(s);
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    x : Standard_Complex_Vectors.Vector(1..n)
      := Standard_Random_Vectors.Random_Vector(1,n);
    y : Standard_Complex_Vectors.Vector(1..n);
    z : Standard_Complex_Vectors.Vector(0..n);
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
      Indexed_Evaluate_and_Differentiate(idx,x,z,y,A);
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
  end Standard_Performance_Test;

  procedure DoblDobl_Performance_Test ( n : integer32 ) is

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem,
  --   in double double precision.

    use DoblDobl_Complex_Poly_SysFun;
    use DoblDobl_Complex_Jaco_Matrices;

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    idx : Standard_Integer_VecVecs.Array_of_VecVecs(1..n)
        := Indexed_Supports(s);
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..n)
      := DoblDobl_Cyclic_Polynomial_System(s);
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    x : DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    y : DoblDobl_Complex_Vectors.Vector(1..n);
    A : DoblDobl_Complex_Matrices.Matrix(1..n,1..n);
    timer : Timing_Widget;
    m : integer32 := 0;

  begin
    new_line;
    put("Give the number of evaluations and differentiations : ");
    get(m);
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      Evaluate_and_Differentiate(s,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"in reverse mode");
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      Indexed_Evaluate_and_Differentiate(idx,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"with precomputed indexes");
    tstart(timer);
    for i in 1..m loop
      x := DoblDobl_Random_Vectors.Random_Vector(1,n);
      Evaluate_and_Differentiate(f,jf,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"with Horner schemes");
  end DoblDobl_Performance_Test;

  procedure QuadDobl_Performance_Test ( n : integer32 ) is

  -- DESCRIPTION :
  --   Does a performance test on the cyclic n-roots problem,
  --   in quad double precision.

    use QuadDobl_Complex_Poly_SysFun;
    use QuadDobl_Complex_Jaco_Matrices;

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    idx : Standard_Integer_VecVecs.Array_of_VecVecs(1..n)
        := Indexed_Supports(s);
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..n)
      := QuadDobl_Cyclic_Polynomial_System(s);
    f : Eval_Poly_Sys(1..n) := Create(p);
    jm : Jaco_Mat(1..n,1..n) := Create(p);
    jf : Eval_Jaco_Mat(1..n,1..n) := Create(jm);
    x : QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    y : QuadDobl_Complex_Vectors.Vector(1..n);
    A : QuadDobl_Complex_Matrices.Matrix(1..n,1..n);
    timer : Timing_Widget;
    m : integer32 := 0;

  begin
    new_line;
    put("Give the number of evaluations and differentiations : ");
    get(m);
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      Evaluate_and_Differentiate(s,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"in reverse mode");
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      Indexed_Evaluate_and_Differentiate(idx,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"with precomputed indexes");
    tstart(timer);
    for i in 1..m loop
      x := QuadDobl_Random_Vectors.Random_Vector(1,n);
      Evaluate_and_Differentiate(f,jf,x,y,A);
    end loop;
    tstop(timer);
    print_times(standard_output,timer,"with Horner schemes");
  end QuadDobl_Performance_Test;

  procedure Standard_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates the support for the cyclic n-roots system and
  --   the coefficients in standard double precision.

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    p : Standard_Complex_Poly_Systems.Poly_Sys(1..n)
      := Standard_Cyclic_Polynomial_System(s);
    ans : character;

  begin
    put("The cyclic "); put(n,1); put_line("-roots system :"); put(p);
    put_line("The supports : ");
    for i in s'range loop
      put(s(i).all); new_line;
    end loop;
    put("Do a random point test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Random_Point_Test(n);
    end if;
    put("Do a performance test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then Standard_Performance_Test(n);
    end if;
  end Standard_Test;

  procedure DoblDobl_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates the support for the cyclic n-roots system and
  --   the coefficients in double double precision.

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..n)
      := DoblDobl_Cyclic_Polynomial_System(s);
    ans : character;

  begin
    put("The cyclic "); put(n,1); put_line("-roots system :"); put(p);
    put_line("The supports : ");
    for i in s'range loop
      put(s(i).all); new_line;
    end loop;
    put("Do a random point test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Random_Point_Test(n);
    end if;
    put("Do a performance test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then DoblDobl_Performance_Test(n);
    end if;
  end DoblDobl_Test;

  procedure QuadDobl_Test ( n : in integer32 ) is

  -- DESCRIPTION :
  --   Generates the support for the cyclic n-roots system and
  --   the coefficients in quad double precision.

    s : Standard_Natural_VecVecs.Array_of_VecVecs(1..n)
      := Supports_of_Cyclic(n);
    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..n)
      := QuadDobl_Cyclic_Polynomial_System(s);
    ans : character;

  begin
    put("The cyclic "); put(n,1); put_line("-roots system :"); put(p);
    put_line("The supports : ");
    for i in s'range loop
      put(s(i).all); new_line;
    end loop;
    put("Do a random point test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Random_Point_Test(n);
    end if;
    put("Do a performance test ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then QuadDobl_Performance_Test(n);
    end if;
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimension and the precision.
  --
    n : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Testing evaluation and differentiation of cyclic n system.");
    new_line;
    put("Give the dimension : "); get(n);
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. standard double precision; or");
    put_line("  1. double double precision; or");
    put_line("  2. quad double precision.");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    case ans is
      when '0' => Standard_Test(n);
      when '1' => DoblDobl_Test(n);
      when '2' => QuadDobl_Test(n);
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_cycfun;
