with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
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
with Standard_Complex_Polynomials;      use Standard_Complex_Polynomials;
with Symbol_Table;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Random_Polynomials;       use Standard_Random_Polynomials;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Complex_Poly_Systems;     use Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_Systems_io;  use Standard_Complex_Poly_Systems_io;
with Standard_Complex_Poly_SysFun;      use Standard_Complex_Poly_SysFun;
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

    res : Standard_Complex_Vectors.Vector(p'range);
    px : constant Standard_Complex_Vectors.Vector := Eval(p,x);
    qx : constant Standard_Complex_Vectors.Vector := Eval(q,x);

  begin
    for i in p'range loop
      res(i) := (1.0-t)*px(i) + t*qx(i);
    end loop;
    return res;
  end Eval;

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
  end Eval;

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

  procedure Random_System_Test is

  -- DESCRIPTION :
  --   Performs the evaluation test on randomly generated systems.

    n,m,d : natural32 := 0;
    p,q : Poly;

  begin
    put("Give the number variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the number of monomials : "); get(m);
    put("Give upper bound on degree : "); get(d);
    declare
      p,q : Poly_Sys(1..integer32(n));
    begin
      for i in p'range loop
        p(i) := Random_Sparse_Poly(n,d,m,0);
        q(i) := Random_Sparse_Poly(n,d,m,0);
      end loop;
      put_line("-> p = "); put(p);
      put_line("-> q = "); put(q);
      Test_System_Evaluation(n,p,q);
    end;
  end Random_System_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Tests (1-t)*p + t*q for given or random polynomials.

    ans : character;

  begin
    put_line("Evaluating a homotopy pair of polynomials ...");
    put_line("  1. interactive test on user given polynomials;");
    put_line("  2. test on randomly generated polynomials;");
    put_line("  3. test on randomly generated polynomial systems.");
    put("Type 1, 2 or 3 to choose a test : ");
    Ask_Alternative(ans,"123");
    case ans is
      when '1' => Interactive_Test;
      when '2' => Random_Test;
      when '3' => Random_System_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_evalhomt;
