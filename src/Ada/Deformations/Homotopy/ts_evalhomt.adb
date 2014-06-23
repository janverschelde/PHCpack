with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Random_Numbers;
with Standard_Integer_Vectors;
with Standard_Integer_Vectors_io;       use Standard_Integer_Vectors_io;
with Standard_Integer_VecVecs;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Random_Vectors;
with Standard_Complex_Matrices;         use Standard_Complex_Matrices;
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Symbol_Table;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
with Standard_Complex_Poly_Randomizers;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;    use Standard_Complex_Jaco_Matrices;
with Standard_Homotopy;
with Standard_Coefficient_Homotopy;     use Standard_Coefficient_Homotopy;

procedure ts_evalhomt is

-- DESCRIPTION :
--   Development of the evaluation of (1-t)*f + t*g
--   where monomials of f and g are shared
--   as in a coefficient-parameter homotopy.

  function Eval ( p,q : Poly; x : Standard_Complex_Vectors.Vector;
                  t : double_float ) return Complex_Number is

    res : Complex_Number;
    px : constant Complex_Number := Eval(p,x);
    qx : constant Complex_Number := Eval(q,x);

  begin
    res := (1.0-t)*px + t*qx;
    return res;
  end Eval;

  function Eval ( p,q : Poly_Sys;
                  x : Standard_Complex_Vectors.Vector;
                  t : double_float )
                return Standard_Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the value of (1-t)*p + t*q, evaluated at x.

    res : Standard_Complex_Vectors.Vector(p'range);
    px : constant Standard_Complex_Vectors.Vector := Eval(p,x);
    qx : constant Standard_Complex_Vectors.Vector := Eval(q,x);

  begin
    for i in p'range loop
      res(i) := (1.0-t)*px(i) + t*qx(i);
    end loop;
    return res;
  end Eval;

  function Diff ( p,q : Poly_Sys;
                  x : Standard_Complex_Vectors.Vector;
                  t : double_float ) return Matrix is

  -- DESCRIPTION :
  --   Returns the Jacobian matrix of (1-t)*p + t*q, evaluated at x.

    res : Matrix(p'range,x'range);
    pjm : Jaco_Mat(p'range,x'range) := Create(p);
    qjm : Jaco_Mat(p'range,x'range) := Create(q);

  begin
    for i in res'range(1) loop
      for j in res'range(2) loop
        res(i,j) := (1.0-t)*Eval(pjm(i,j),x) + t*Eval(qjm(i,j),x);
      end loop;
    end loop;
    return res;
  end Diff;

  procedure Eval ( n : in natural32; p,q,h : in Poly;
                   cp,cq,ch : in Standard_Complex_Vectors.Vector;
                   ip,iq : in Standard_Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables;
  --   p        first polynomial in several variables;
  --   q        second polynomial in several variables;
  --   h        labeled structure of the sum of p and q;
  --   cp       coefficient vector of p;
  --   cq       coefficient vector of q;
  --   ch       coefficient vector of h;
  --   ip       index of the labels of coefficients of p in h:
  --            ip(k) locates in ch the k-th coefficient of p,
  --            in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq       index of the labels of coefficients of q in h:
  --            ip(k) locates in ch the k-th coefficient of q,
  --            in particular: ch(iq(k)) := cq(k) sets h to q.

    t : constant double_float := abs(Standard_Random_Numbers.Random);
    x : constant Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    y : constant Complex_Number := Eval(p,q,x,t);
    f : Eval_Coeff_Poly := Create(h);
    c : constant Standard_Complex_Vectors.Vector
      := Evaluated_Coefficients(ch'last,cp,cq,ip,iq,t);
    z : constant Complex_Number := Eval(f,c,x);

  begin
    put("A random t : "); put(t); new_line;
    put_line("A random point : "); put_line(x);
    put("-> y = "); put(y); new_line;
    put("-> z = "); put(z); new_line;
  end Eval;

  procedure Write_Elements ( A,B : in Matrix ) is

  -- DESCRIPTION :
  --   Writes the elements of the matrices A and B to screen.

  begin
    for i in A'range(1) loop
      for j in A'range(2) loop
        put("A("); put(i,1); put(","); put(j,1); put(") = ");
        put(A(i,j)); new_line;
        put("B("); put(i,1); put(","); put(j,1); put(") = ");
        put(B(i,j)); new_line;
      end loop;
    end loop;
  end Write_Elements;

  procedure Eval ( n : in natural32; p,q,h : in Poly_Sys;
                   cp,cq,ch : in Standard_Complex_VecVecs.VecVec;
                   ip,iq : in Standard_Integer_VecVecs.VecVec ) is

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables;
  --   p        first polynomial in several variables;
  --   q        second polynomial in several variables;
  --   h        labeled structure of the sum of p and q;
  --   cp       coefficient vector of p;
  --   cq       coefficient vector of q;
  --   ch       coefficient vector of h;
  --   ip       index of the labels of coefficients of p in h:
  --            ip(k) locates in ch the k-th coefficient of p,
  --            in particular: ch(ip(k)) := cp(k) sets h to p;
  --   iq       index of the labels of coefficients of q in h:
  --            ip(k) locates in ch the k-th coefficient of q,
  --            in particular: ch(iq(k)) := cq(k) sets h to q.

    t : constant double_float := abs(Standard_Random_Numbers.Random);
    x : constant Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    y : constant Standard_Complex_Vectors.Vector := Eval(p,q,x,t);
    f : Eval_Coeff_Poly_Sys(p'range) := Create(h);
    c : Standard_Complex_VecVecs.VecVec(f'range);
    z : Standard_Complex_Vectors.Vector(f'range);
    m : Mult_Factors(p'range,x'range);
    jf : Eval_Coeff_Jaco_Mat(p'range,x'range);
    A : Matrix(p'range,x'range) := Diff(p,q,x,t);
    B : Matrix(p'range,x'range);

  begin
    put("A random t : "); put(t); new_line;
    put_line("A random point : "); put_line(x);
    put_line("-> y = "); put_line(y);
    for i in c'range loop
      declare
        cff : Standard_Complex_Vectors.Vector(ch(i)'range);
      begin
        cff := (cff'range => Create(0.0));
        c(i) := new Standard_Complex_Vectors.Vector'(cff);
      end;
    end loop;
    Evaluated_Coefficients(c,cp,cq,ip,iq,t);
    z := Eval(f,c,x);
    put_line("-> z = "); put_line(z);
    Standard_Complex_Jaco_Matrices.Create(h,jf,m);
    B := Eval(jf,m,c,x);
    Write_Elements(A,B);
  end Eval;

  procedure Encapsulated_Eval
              ( n : in natural32; p,q : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables;
  --   p        first polynomial in several variables;
  --   q        second polynomial in several variables.

    t : constant double_float := abs(Standard_Random_Numbers.Random);
    ct : constant Complex_Number := Create(t);
    x : constant Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    y : constant Standard_Complex_Vectors.Vector := Eval(p,q,x,t);
    z : constant Standard_Complex_Vectors.Vector
      := Standard_Coefficient_Homotopy.Eval(x,ct);
    A : Matrix(p'range,x'range) := Diff(p,q,x,t);
    B : Matrix(p'range,x'range)
      := Standard_Coefficient_Homotopy.Diff(x,ct);

  begin
    put("A random t : "); put(ct); new_line;
    put_line("A random point : "); put_line(x);
    put_line("-> y = "); put_line(y);
    put_line("-> z = "); put_line(z);
    Write_Elements(A,B);
  end Encapsulated_Eval;

  procedure Compared_Encapsulated_Eval ( n : in natural32 ) is

  -- DESCRIPTION :
  --   Evaluates (1-t)*p + t*q at some random t in [0,1]
  --   and at a random vector.

  -- ON ENTRY :
  --   n        number of variables.

    t : constant double_float := abs(Standard_Random_Numbers.Random);
    ct : constant Complex_Number := Create(t);
    x : constant Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    y : constant Standard_Complex_Vectors.Vector
      := Standard_Homotopy.Eval(x,ct);
    z : constant Standard_Complex_Vectors.Vector
      := Standard_Coefficient_Homotopy.Eval(x,ct);
    A : Matrix(x'range,x'range) := Standard_Homotopy.Diff(x,ct);
    B : Matrix(x'range,x'range)
      := Standard_Coefficient_Homotopy.Diff(x,ct);

  begin
    put("A random t : "); put(ct); new_line;
    put_line("A random point : "); put_line(x);
    put_line("-> y = "); put_line(y);
    put_line("-> z = "); put_line(z);
    Write_Elements(A,B);
  end Compared_Encapsulated_Eval;

  procedure Test_Evaluation ( n : in natural32; p,q : in Poly ) is

  -- DESCRIPTION :
  --   Tests the evaluation of (1-t)*p + t*q for two polynomials
  --   p and q in n variables.

    cp : constant Standard_Complex_Vectors.Vector := Coefficients(p);
    cq : constant Standard_Complex_Vectors.Vector := Coefficients(q);
    lp : Poly := Labeled_Coefficients(p,true);
    lq : Poly := Labeled_Coefficients(q,false);
    lh : Poly := lp + lq;

  begin
    put_line("The coefficients of p : "); put_line(cp);
    put_line("The coefficients of q : "); put_line(cq);
    put_line("-> labeled p : "); put(lp); new_line;
    put_line("-> labeled q : "); put(lq); new_line;
    put_line("-> labeled h : "); put(lh); new_line;
    declare
      ch : constant Standard_Complex_Vectors.Vector := Coefficients(lh);
      ip : constant Standard_Integer_Vectors.Vector
         := Index_of_Labels(ch,true);
      iq : constant Standard_Integer_Vectors.Vector
         := Index_of_Labels(ch,false);
    begin
      put_line("The coefficients of labeled h :"); put_line(ch);
      put("indices for p : "); put(ip); new_line;
      put("indices for q : "); put(iq); new_line;
      Eval(n,p,q,lh,cp,cq,ch,ip,iq);
    end;
  end Test_Evaluation;

  procedure Test_System_Evaluation 
              ( n : in natural32; p,q : in Poly_Sys ) is

  -- DESCRIPTION :
  --   Tests the evaluation of (1-t)*p + t*q for two polynomial
  --   systems p and q in n variables.

    cp,cq,ch : Standard_Complex_VecVecs.VecVec(p'range);
    lp,lq,lh : Poly_Sys(p'range);
    ip,iq : Standard_Integer_VecVecs.VecVec(p'range);

  begin
    cp := Coefficients(p);
    cq := Coefficients(q);
    lp := Labeled_Coefficients(p,true);
    lq := Labeled_Coefficients(q,false);
    lh := lp + lq;
    ch := Coefficients(lh);
    ip := Index_of_Labels(ch,true);
    iq := Index_of_Labels(ch,false);
    Eval(n,p,q,lh,cp,cq,ch,ip,iq);
  end Test_System_Evaluation;

  procedure Interactive_Test is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables
  --   and then for two polynomials in that number of variables,
  --   before running the evaluation test.

    n : natural32 := 0;
    p,q : Poly;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    put_line("Reading the first polynomial p ..."); get(p);
    put("-> p = "); put(p); new_line;
    put_line("Reading the second polynomial q ..."); get(q);
    put("-> q = "); put(q); new_line;
    Test_Evaluation(n,p,q);
  end Interactive_Test;

  procedure Random_Test is

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated polynomials.

    n,m,d : natural32 := 0;
    p,q : Poly;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    p := Random_Sparse_Poly(n,d,m,0);
    put("-> p = "); put(p); new_line;
    q := Random_Sparse_Poly(n,d,m,0);
    put("-> q = "); put(q); new_line;
    Test_Evaluation(n,p,q);
  end Random_Test;

  procedure Random_Systems
              ( n : in integer32; p,q : out Poly_Sys ) is

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q.

    m,d : natural32 := 0;

  begin
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    for i in p'range loop
      p(i) := Random_Sparse_Poly(natural32(n),d,m,0);
      q(i) := Random_Sparse_Poly(natural32(n),d,m,0);
    end loop;
    put_line("-> p = "); put(p);
    put_line("-> q = "); put(q);
  end Random_Systems;

  procedure Random_Coefficient_Systems
              ( n : in integer32; p,q : out Poly_Sys ) is

  -- DESCRIPTION :
  --   Generates two n-dimensional systems p and q,
  --   with the same supports.

    m,d : natural32 := 0;

  begin
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    for i in p'range loop
      p(i) := Random_Sparse_Poly(natural32(n),d,m,0);
      q(i) := Standard_Complex_Poly_Randomizers.Complex_Randomize1(p(i));
    end loop;
    put_line("-> p = "); put(p);
    put_line("-> q = "); put(q);
  end Random_Coefficient_Systems;

  procedure Random_System_Test is

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems.

    n : natural32 := 0;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : Poly_Sys(1..integer32(n));
    begin
      Random_Systems(integer32(n),p,q);
      Test_System_Evaluation(n,p,q);
    end;
  end Random_System_Test;

  procedure Random_Encapsulation_Test is

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   using the encapsulated interface.

    n : natural32 := 0;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : Poly_Sys(1..integer32(n));
      gamma : Complex_Number := Create(1.0);
    begin
      Random_Systems(integer32(n),p,q);
      Standard_Coefficient_Homotopy.Create(p,q,1,gamma);
      Encapsulated_Eval(n,p,q);
    end;
  end Random_Encapsulation_Test;

  procedure Compared_Encapsulation_Test is

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the Standard_Homotopy package.

    n : natural32 := 0;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : Poly_Sys(1..integer32(n));
      gamma : Complex_Number := Standard_Random_Numbers.Random1;
    begin
      Random_Systems(integer32(n),p,q);
      Standard_Coefficient_Homotopy.Create(p,q,2,gamma);
      Standard_Homotopy.Create(q,p,2,gamma);
      Compared_Encapsulated_Eval(n);
    end;
  end Compared_Encapsulation_Test;

  procedure Standard_Homotopy_Performance ( n,m : natural32 ) is

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the homotopy with timings taken.

    t : double_float;
    ct : Complex_Number;
    x,y : Standard_Complex_Vectors.Vector(1..integer32(n));
    A : Matrix(x'range,x'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for i in 1..m loop
      t := abs(Standard_Random_Numbers.Random);
      ct := Create(t);
      x := Standard_Random_Vectors.Random_Vector(1,integer32(n));
      y := Standard_Homotopy.Eval(x,ct);
      A := Standard_Homotopy.Diff(x,ct);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"eval & diff of standard homotopy");
  end Standard_Homotopy_Performance;

  procedure Standard_Coefficient_Homotopy_Performance ( n,m : natural32 ) is

  -- DESCRIPTION :
  --   Generates m random values for x and t
  --   for evaluation in the coefficient homotopy with timings taken.

    t : double_float;
    ct : Complex_Number;
    x,y : Standard_Complex_Vectors.Vector(1..integer32(n));
    A : Matrix(x'range,x'range);
    timer : Timing_Widget;

  begin
    tstart(timer);
    for i in 1..m loop
      t := abs(Standard_Random_Numbers.Random);
      ct := Create(t);
      x := Standard_Random_Vectors.Random_Vector(1,integer32(n));
      y := Standard_Coefficient_Homotopy.Eval(x,ct);
      A := Standard_Coefficient_Homotopy.Diff(x,ct);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"eval & diff of coefficient homotopy");
  end Standard_Coefficient_Homotopy_Performance;

  procedure Performance_Test is

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems,
  --   comparing with the Standard_Homotopy package.

    n,m : natural32 := 0;

  begin
    put("Give the number of evaluations : "); get(m);
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    declare
      p,q : Poly_Sys(1..integer32(n));
      gamma : Complex_Number := Standard_Random_Numbers.Random1;
    begin
      Random_Coefficient_Systems(integer32(n),p,q);
      Standard_Coefficient_Homotopy.Create(p,q,2,gamma);
      Standard_Homotopy.Create(q,p,2,gamma);
    end;
    Standard_Homotopy_Performance(n,m);
    Standard_Coefficient_Homotopy_Performance(n,m);
  end Performance_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Tests (1-t)*p + t*q for given or random polynomials.

    ans : character;

  begin
    new_line;
    put_line("Evaluating a homotopy pair of polynomials ...");
    put_line("  1. interactive test on user given polynomials;");
    put_line("  2. test on randomly generated polynomials;");
    put_line("  3. test on randomly generated polynomial systems;");
    put_line("  4. test on encapsulation for random systems;");
    put_line("  5. compare with homotopy package for random systems;");
    put_line("  6. performance test on evaluations of homotopies.");
    put("Type 1, 2, 3, 4, 5, or 6 to choose a test : ");
    Ask_Alternative(ans,"123456");
    new_line;
    case ans is
      when '1' => Interactive_Test;
      when '2' => Random_Test;
      when '3' => Random_System_Test;
      when '4' => Random_Encapsulation_Test;
      when '5' => Compared_Encapsulation_Test;
      when '6' => Performance_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_evalhomt;
