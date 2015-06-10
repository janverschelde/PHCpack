with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Mathematical_Functions;   use Standard_Mathematical_Functions;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;       use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;      use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Norms_Equals;     use Standard_Complex_Norms_Equals;
with Standard_Random_Numbers;
with Standard_Random_Vectors;           use Standard_Random_Vectors;
with Standard_Complex_Matrices;
with Standard_Speelpenning_Products;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Double_Double_Numbers_io;          use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Mathematical_Functions;   use DoblDobl_Mathematical_Functions;
with Double_Double_Vectors;
with Double_Double_Vectors_io;          use Double_Double_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Vector_Norms;     use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with DoblDobl_Complex_Matrices;
with DoblDobl_Speelpenning_Products;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Quad_Double_Numbers_io;            use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Mathematical_Functions;   use QuadDobl_Mathematical_Functions;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;            use Quad_Double_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Vector_Norms;     use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with QuadDobl_Complex_Matrices;
with QuadDobl_Speelpenning_Products;
with Multprec_Floating_Numbers;         use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;      use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;       use Multprec_Complex_Numbers_io;
with Multprec_Mathematical_Functions;   use Multprec_Mathematical_Functions;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;      use Multprec_Floating_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;       use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;     use Multprec_Complex_Norms_Equals;
with Multprec_Random_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Complex_Matrices;
with Multprec_Speelpenning_Products;
with Standard_Monomial_Evaluations;
with DoblDobl_Monomial_Evaluations;
with QuadDobl_Monomial_Evaluations;
with Multprec_Monomial_Evaluations;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Polynomials_io;   use QuadDobl_Complex_Polynomials_io;
with QuadDobl_Complex_Poly_Functions;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Polynomials_io;   use Multprec_Complex_Polynomials_io;
with Multprec_Complex_Poly_Functions;
with Coefficient_Supported_Polynomials; use Coefficient_Supported_Polynomials;
with Standard_Gradient_Evaluations;
with DoblDobl_Gradient_Evaluations;
with QuadDobl_Gradient_Evaluations;
with Multprec_Gradient_Evaluations;
with Standard_Complex_Poly_Systems;
with Standard_Complex_Poly_SysFun;
with Standard_Complex_Jaco_Matrices;
with DoblDobl_Complex_Poly_Systems;
with DoblDobl_Complex_Poly_SysFun;
with DoblDobl_Complex_Jaco_Matrices;
with QuadDobl_Complex_Poly_Systems;
with QuadDobl_Complex_Poly_SysFun;
with QuadDobl_Complex_Jaco_Matrices;
with Multprec_Complex_Poly_Systems;
with Multprec_Complex_Poly_SysFun;
with Multprec_Complex_Jaco_Matrices;
with Cyclic_Roots_System;               use Cyclic_Roots_System;
with Random_Conditioned_Matrices;       use Random_Conditioned_Matrices;
with Random_Conditioned_Evaluations;    use Random_Conditioned_Evaluations;
with Varbprec_Complex_Linear_Solvers;   use Varbprec_Complex_Linear_Solvers;
with Varbprec_Polynomial_Evaluations;   use Varbprec_Polynomial_Evaluations;
with Varbprec_Gradient_Evaluations;     use Varbprec_Gradient_Evaluations;

procedure ts_vmpdiff is

-- DESCRIPTION :
--   Elaboration of gradient computation with condition numbers
--   at variable levels of precision.
--   This test program is an extension of ts_speel and ts_vmpeval.

  function Random_Exponent
             ( n,d : integer32 ) return Standard_Natural_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns an exponent vector in n variables
  --   with entries ranging between 0 and d.

    res : Standard_Natural_Vectors.Vector(1..n);

  begin
    for i in 1..n loop
      res(i) := natural32(Standard_Random_Numbers.Random(0,d));
    end loop;
    return res;
  end Random_Exponent;

  function Random_Exponents
             ( n,d,m : integer32 ) return Standard_Natural_VecVecs.VecVec is

  -- DESCRIPTION :
  --   Returns m random exponents in n variables with entries
  --   ranging between 0 and d.

    res : Standard_Natural_VecVecs.VecVec(1..m);
    exp : Standard_Natural_Vectors.Vector(1..n);

  begin
    for i in 1..m loop
      exp := Random_Exponent(n,d);
      res(i) := new Standard_Natural_Vectors.Vector'(exp);
    end loop;
    return res;
  end Random_Exponents;

  procedure Standard_Evaluate_Gradient ( n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Randomly generates a sequence of m monomials in n variables
  --   of degrees at most d, using standard arithmetic.

    e : Standard_Natural_VecVecs.VecVec(1..m) := Random_Exponents(n,d,m);
    c : constant Standard_Complex_Vectors.Vector(1..m) := Random_Vector(1,m);
    f,b : Standard_Natural_VecVecs.VecVec(1..m);
    wrk : Standard_Complex_VecVecs.VecVec(1..m);
    x : constant Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    z,fgz : Standard_Complex_Vectors.Vector(0..n);
    numfz,denfz,fzrco,maxng,mindg,gzrco : double_float;
    p : Standard_Complex_Polynomials.Poly := Create_Standard_Polynomial(c,e);
    fz : Standard_Complex_Numbers.Complex_Number;
    absfz,denrco,rco : double_float;

    use Standard_Gradient_Evaluations;

  begin
    Split_Common_Factors(e,f,b);
    put_line("The exponents, with splitted factors : ");
    for i in 1..m loop
      put(e(i).all);
      put(" = "); put(f(i).all);
      put(" + "); put(b(i).all); new_line;
    end loop;
    z := Gradient_of_Polynomial(f,b,c,x);
    for i in b'range loop
      wrk(i) := new Standard_Complex_Vectors.Vector(0..n);
    end loop;
    Gradient_with_Inverse_Condition
      (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,maxng,mindg,gzrco);
    Inverse_Condition_Number(p,x,fz,absfz,denrco,rco);
    put_line("Function value with three different procedures : ");
    put("f(z) : "); put(z(0)); new_line;
    put("f(z) : "); put(fgz(0)); new_line;
    put("f(z) : "); put(fz); new_line;
    put_line("Condition number of evaluation by different procedures :");
    put("inverse condition number : "); put(fzrco,3); new_line;
    put("inverse condition number : "); put(rco,3); new_line;
    put_line("Condition of gradient :");
    put("  max of numerator : "); put(maxng,3); new_line;
    put("  min of denominator : "); put(mindg,3); new_line;
    put("  inverse condition number : "); put(gzrco,3); new_line;
    Standard_Natural_VecVecs.Clear(e);
    Standard_Natural_VecVecs.Clear(f);
    Standard_Natural_VecVecs.Clear(b);
    Standard_Complex_VecVecs.Clear(wrk);
  end Standard_Evaluate_Gradient;

  procedure DoblDobl_Evaluate_Gradient ( n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Randomly generates a sequence of m monomials in n variables
  --   of degrees at most d, using double double arithmetic.

    e : Standard_Natural_VecVecs.VecVec(1..m) := Random_Exponents(n,d,m);
    c : constant DoblDobl_Complex_Vectors.Vector(1..m)
      := DoblDobl_Random_Vectors.Random_Vector(1,m);
    f,b : Standard_Natural_VecVecs.VecVec(1..m);
    wrk : DoblDobl_Complex_VecVecs.VecVec(1..m);
    x : constant DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    z,fgz : DoblDobl_Complex_Vectors.Vector(0..n);
    numfz,denfz,fzrco,maxng,mindg,gzrco : double_double;
    p : DoblDobl_Complex_Polynomials.Poly := Create_DoblDobl_Polynomial(c,e);
    fz : DoblDobl_Complex_Numbers.Complex_Number;
    absfz,denrco,rco : double_double;

    use DoblDobl_Gradient_Evaluations;

  begin
    Split_Common_Factors(e,f,b);
    put_line("The exponents, with splitted factors : ");
    for i in 1..m loop
      put(e(i).all);
      put(" = "); put(f(i).all);
      put(" + "); put(b(i).all); new_line;
    end loop;
    z := Gradient_of_Polynomial(f,b,c,x);
    for i in b'range loop
      wrk(i) := new DoblDobl_Complex_Vectors.Vector(0..n);
    end loop;
    Gradient_with_Inverse_Condition
      (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,maxng,mindg,gzrco);
    Inverse_Condition_Number(p,x,fz,absfz,denrco,rco);
    put_line("Function value with three different procedures : ");
    put("f(z) : "); put(z(0)); new_line;
    put("f(z) : "); put(fgz(0)); new_line;
    put("f(z) : "); put(fz); new_line;
    put_line("Condition number of evaluation by different procedures :");
    put("inverse condition number : "); put(fzrco,3); new_line;
    put("inverse condition number : "); put(rco,3); new_line;
    put_line("Condition of gradient :");
    put("  max of numerator : "); put(maxng,3); new_line;
    put("  min of denominator : "); put(mindg,3); new_line;
    put("  inverse condition number : "); put(gzrco,3); new_line;
    Standard_Natural_VecVecs.Clear(b);
    Standard_Natural_VecVecs.Clear(f);
    Standard_Natural_VecVecs.Clear(e);
    DoblDobl_Complex_VecVecs.Clear(wrk);
  end DoblDobl_Evaluate_Gradient;

  procedure QuadDobl_Evaluate_Gradient ( n,d,m : in integer32 ) is

  -- DESCRIPTION :
  --   Randomly generates a sequence of m monomials in n variables
  --   of degrees at most d, using quad double arithmetic.

    e : Standard_Natural_VecVecs.VecVec(1..m) := Random_Exponents(n,d,m);
    c : constant QuadDobl_Complex_Vectors.Vector(1..m)
      := QuadDobl_Random_Vectors.Random_Vector(1,m);
    f,b : Standard_Natural_VecVecs.VecVec(1..m);
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..m);
    x : constant QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    z,fgz : QuadDobl_Complex_Vectors.Vector(0..n);
    numfz,denfz,fzrco,maxng,mindg,gzrco : quad_double;
    p : QuadDobl_Complex_Polynomials.Poly := Create_QuadDobl_Polynomial(c,e);
    fz : QuadDobl_Complex_Numbers.Complex_Number;
    absfz,denrco,rco : quad_double;

    use QuadDobl_Gradient_Evaluations;

  begin
    Split_Common_Factors(e,f,b);
    put_line("The exponents, with splitted factors : ");
    for i in 1..m loop
      put(e(i).all);
      put(" = "); put(f(i).all);
      put(" + "); put(b(i).all); new_line;
    end loop;
    z := Gradient_of_Polynomial(f,b,c,x);
    for i in b'range loop
      wrk(i) := new QuadDobl_Complex_Vectors.Vector(0..n);
    end loop;
    Gradient_with_Inverse_Condition
      (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,maxng,mindg,gzrco);
    Inverse_Condition_Number(p,x,fz,absfz,denrco,rco);
    put_line("Function value with three different procedures : ");
    put("f(z) : "); put(z(0)); new_line;
    put("f(z) : "); put(fgz(0)); new_line;
    put("f(z) : "); put(fz); new_line;
    put_line("Condition number of evaluation by different procedures :");
    put("inverse condition number : "); put(fzrco,3); new_line;
    put("inverse condition number : "); put(rco,3); new_line;
    put_line("Condition of gradient :");
    put("  max of numerator : "); put(maxng,3); new_line;
    put("  min of denominator : "); put(mindg,3); new_line;
    put("  inverse condition number : "); put(gzrco,3); new_line;
    Standard_Natural_VecVecs.Clear(e);
    Standard_Natural_VecVecs.Clear(f);
    Standard_Natural_VecVecs.Clear(b);
    QuadDobl_Complex_VecVecs.Clear(wrk);
  end QuadDobl_Evaluate_Gradient;

  procedure Multprec_Evaluate_Gradient
              ( n,d,m : in integer32; size : in natural32 ) is

  -- DESCRIPTION :
  --   Randomly generates a sequence of m monomials in n variables
  --   of degrees at most d, with size long multiprecision numbers.

    e : Standard_Natural_VecVecs.VecVec(1..m) := Random_Exponents(n,d,m);
    c : Multprec_Complex_Vectors.Vector(1..m)
      := Multprec_Random_Vectors.Random_Vector(1,m,size);
    f,b : Standard_Natural_VecVecs.VecVec(1..m);
    wrk : Multprec_Complex_VecVecs.VecVec(1..m);
    x : Multprec_Complex_Vectors.Vector(1..n)
      := Multprec_Random_Vectors.Random_Vector(1,n,size);
    z,fgz : Multprec_Complex_Vectors.Vector(0..n);
    numfz,denfz,fzrco,maxng,mindg,gzrco : Floating_Number;
    p : Multprec_Complex_Polynomials.Poly := Create_Multprec_Polynomial(c,e);
    fz : Multprec_Complex_Numbers.Complex_Number;
    absfz,denrco,rco : Floating_Number;

    use Multprec_Gradient_Evaluations;

  begin
    Split_Common_Factors(e,f,b);
    put_line("The exponents, with splitted factors : ");
    for i in 1..m loop
      put(e(i).all);
      put(" = "); put(f(i).all);
      put(" + "); put(b(i).all); new_line;
    end loop;
    z := Gradient_of_Polynomial(f,b,c,x);
    for i in b'range loop
      wrk(i) := new Multprec_Complex_Vectors.Vector(0..n);
    end loop;
    Gradient_with_Inverse_Condition
      (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,maxng,mindg,gzrco);
    Inverse_Condition_Number(p,x,fz,absfz,denrco,rco);
    put_line("Function value with three different procedures : ");
    put("f(z) : "); put(z(0)); new_line;
    put("f(z) : "); put(fgz(0)); new_line;
    put("f(z) : "); put(fz); new_line;
    put_line("Condition number of evaluation by different procedures :");
    put("inverse condition number : "); put(fzrco,3); new_line;
    put("inverse condition number : "); put(rco,3); new_line;
    put_line("Condition of gradient :");
    put("  max of numerator : "); put(maxng,3); new_line;
    put("  min of denominator : "); put(mindg,3); new_line;
    put("  inverse condition number : "); put(gzrco,3); new_line;
    Standard_Natural_VecVecs.Clear(e);
    Standard_Natural_VecVecs.Clear(f);
    Standard_Natural_VecVecs.Clear(b);
    Multprec_Complex_VecVecs.Clear(wrk);
  end Multprec_Evaluate_Gradient;

  procedure Evaluate_Gradient is

  -- DESCRIPTION :
  --   Interactive test for the development of the gradient evaluation of
  --   monomials at some random vectors, with and without power table.

    n,d,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put_line("Reading parameters for a random exponent sequence ...");
    put("  give the number of variables, n = "); get(n);
    put("  give the largest degree, d = "); get(d);
    put("  give number of monomials, m = "); get(m);
    put("double, double double, quad double, or multiprecision? (s/d/q/m) ");
    Ask_Alternative(ans,"sdqm");
    case ans is
      when 's' => Standard_Evaluate_Gradient(n,d,m);
      when 'd' => DoblDobl_Evaluate_Gradient(n,d,m);
      when 'q' => QuadDobl_Evaluate_Gradient(n,d,m);
      when 'm' =>
        declare
          deci,size : natural32 := 0;
        begin
          new_line;
          put("Give the number of decimal places : "); get(deci);
          size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
          Multprec_Evaluate_Gradient(n,d,m,size);
        end;
      when others => null;
    end case;
  end Evaluate_Gradient;

  procedure Standard_Performance_Test ( n,m : integer32 ) is

  -- DESCRIPTION :
  --   Generates the next to last equation of the cyclic n-roots problem
  --   with standard complex coefficients and performs m evaluations
  --   of the gradient, once in the reverse mode and once straightforward.

    s : Standard_Natural_VecVecs.VecVec(1..n) := Support_of_Cyclic(n,n-1);
    f,b : Standard_Natural_VecVecs.VecVec(1..n);
    p : Standard_Complex_Polynomials.Poly := Standard_Cyclic_Polynomial(s);
    c : Standard_Complex_Vectors.Vector(1..n)
      := (1..n => Standard_Complex_Numbers.Create(1.0));
    x : Standard_Complex_Vectors.Vector(1..n) := Random_Vector(1,n);
    z : Standard_Complex_Vectors.Vector(0..n);
    g : Standard_Complex_Poly_Systems.Poly_Sys(0..n);
    eg : Standard_Complex_Poly_SysFun.Eval_Poly_Sys(0..n);
    wrk : Standard_Complex_VecVecs.VecVec(1..n);
    timer : Timing_Widget;

    use Standard_Gradient_Evaluations;

  begin
    put_line("The polynomial : "); put(p);
    Split_Common_Factors(s,f,b);
    tstart(timer);
    for i in 1..m loop
      z := Gradient_of_Polynomial(f,b,c,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient function in reverse mode");
    for i in b'range loop
      wrk(i) := new Standard_Complex_Vectors.Vector(0..n);
    end loop;
    tstart(timer);
    for i in 1..m loop
      Gradient_of_Polynomial(f,b,c,x,wrk,z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient procedure in reverse mode");
    for k in 1..n loop
      g(k) := Standard_Complex_Polynomials.Diff(p,k);
    end loop;
    tstart(timer);
    for i in 1..m loop
      z := Standard_Complex_Poly_SysFun.Eval(g,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient plain mode");
    eg := Standard_Complex_Poly_SysFun.Create(g);
    tstart(timer);
    for i in 1..m loop
      z := Standard_Complex_Poly_SysFun.Eval(eg,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient with Horner");
    Standard_Natural_VecVecs.Clear(s);
    Standard_Complex_Polynomials.Clear(p);
    Standard_Complex_Poly_Systems.Clear(g);
    Standard_Complex_Poly_SysFun.Clear(eg);
  end Standard_Performance_Test;

  procedure DoblDobl_Performance_Test ( n,m : integer32 ) is

  -- DESCRIPTION :
  --   Generates the next to last equation of the cyclic n-roots problem
  --   with standard complex coefficients and performs m evaluations
  --   of the gradient, once in the reverse mode and once straightforward.

    s : Standard_Natural_VecVecs.VecVec(1..n) := Support_of_Cyclic(n,n-1);
    f,b : Standard_Natural_VecVecs.VecVec(1..n);
    p : DoblDobl_Complex_Polynomials.Poly := DoblDobl_Cyclic_Polynomial(s);
    c : DoblDobl_Complex_Vectors.Vector(1..n)
      := (1..n => DoblDobl_Complex_Numbers.Create(integer(1)));
    x : DoblDobl_Complex_Vectors.Vector(1..n)
      := DoblDobl_Random_Vectors.Random_Vector(1,n);
    z : DoblDobl_Complex_Vectors.Vector(0..n);
    g : DoblDobl_Complex_Poly_Systems.Poly_Sys(0..n);
    eg : DoblDobl_Complex_Poly_SysFun.Eval_Poly_Sys(0..n);
    wrk : DoblDobl_Complex_VecVecs.VecVec(1..n);
    timer : Timing_Widget;

    use DoblDobl_Gradient_Evaluations;

  begin
    put_line("The polynomial : "); put(p);
    Split_Common_Factors(s,f,b);
    tstart(timer);
    for i in 1..m loop
      z := Gradient_of_Polynomial(f,b,c,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient function in reverse mode");
    for i in b'range loop
      wrk(i) := new DoblDobl_Complex_Vectors.Vector(0..n);
    end loop;
    tstart(timer);
    for i in 1..m loop
      Gradient_of_Polynomial(f,b,c,x,wrk,z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient procedure in reverse mode");
    for k in 1..n loop
      g(k) := DoblDobl_Complex_Polynomials.Diff(p,k);
    end loop;
    tstart(timer);
    for i in 1..m loop
      z := DoblDobl_Complex_Poly_SysFun.Eval(g,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient plain mode");
    eg := DoblDobl_Complex_Poly_SysFun.Create(g);
    tstart(timer);
    for i in 1..m loop
      z := DoblDobl_Complex_Poly_SysFun.Eval(eg,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient with Horner");
    Standard_Natural_VecVecs.Clear(s);
    DoblDobl_Complex_Polynomials.Clear(p);
    DoblDobl_Complex_Poly_Systems.Clear(g);
    DoblDobl_Complex_Poly_SysFun.Clear(eg);
  end DoblDobl_Performance_Test;

  procedure QuadDobl_Performance_Test ( n,m : integer32 ) is

  -- DESCRIPTION :
  --   Generates the next to last equation of the cyclic n-roots problem
  --   with standard complex coefficients and performs m evaluations
  --   of the gradient, once in the reverse mode and once straightforward.

    s : Standard_Natural_VecVecs.VecVec(1..n) := Support_of_Cyclic(n,n-1);
    f,b : Standard_Natural_VecVecs.VecVec(1..n);
    p : QuadDobl_Complex_Polynomials.Poly := QuadDobl_Cyclic_Polynomial(s);
    c : QuadDobl_Complex_Vectors.Vector(1..n)
      := (1..n => QuadDobl_Complex_Numbers.Create(integer(1)));
    x : QuadDobl_Complex_Vectors.Vector(1..n)
      := QuadDobl_Random_Vectors.Random_Vector(1,n);
    z : QuadDobl_Complex_Vectors.Vector(0..n);
    g : QuadDobl_Complex_Poly_Systems.Poly_Sys(0..n);
    eg : QuadDobl_Complex_Poly_SysFun.Eval_Poly_Sys(0..n);
    wrk : QuadDobl_Complex_VecVecs.VecVec(1..n);
    timer : Timing_Widget;

    use QuadDobl_Gradient_Evaluations;

  begin
    put_line("The polynomial : "); put(p);
    Split_Common_Factors(s,f,b);
    tstart(timer);
    for i in 1..m loop
      z := Gradient_of_Polynomial(f,b,c,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient function in reverse mode");
    for i in b'range loop
      wrk(i) := new QuadDobl_Complex_Vectors.Vector(0..n);
    end loop;
    tstart(timer);
    for i in 1..m loop
      Gradient_of_Polynomial(f,b,c,x,wrk,z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient procedure in reverse mode");
    for k in 1..n loop
      g(k) := QuadDobl_Complex_Polynomials.Diff(p,k);
    end loop;
    tstart(timer);
    for i in 1..m loop
      z := QuadDobl_Complex_Poly_SysFun.Eval(g,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient plain mode");
    eg := QuadDobl_Complex_Poly_SysFun.Create(g);
    tstart(timer);
    for i in 1..m loop
      z := QuadDobl_Complex_Poly_SysFun.Eval(eg,x);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient with Horner");
    Standard_Natural_VecVecs.Clear(s);
    QuadDobl_Complex_Polynomials.Clear(p);
    QuadDobl_Complex_Poly_Systems.Clear(g);
    QuadDobl_Complex_Poly_SysFun.Clear(eg);
  end QuadDobl_Performance_Test;

  procedure Multprec_Performance_Test
             ( n,m : integer32; size : natural32 ) is

  -- DESCRIPTION :
  --   Generates the next to last equation of the cyclic n-roots problem
  --   with standard complex coefficients and performs m evaluations
  --   of the gradient, once in the reverse mode and once straightforward.

    s : Standard_Natural_VecVecs.VecVec(1..n) := Support_of_Cyclic(n,n-1);
    f,b : Standard_Natural_VecVecs.VecVec(1..n);
    p : Multprec_Complex_Polynomials.Poly := Multprec_Cyclic_Polynomial(s);
    c : Multprec_Complex_Vectors.Vector(1..n)
      := (1..n => Multprec_Complex_Numbers.Create(integer(1)));
    x : Multprec_Complex_Vectors.Vector(1..n)
      := Multprec_Random_Vectors.Random_Vector(1,n,size);
    z : Multprec_Complex_Vectors.Vector(0..n);
    g : Multprec_Complex_Poly_Systems.Poly_Sys(0..n);
    eg : Multprec_Complex_Poly_SysFun.Eval_Poly_Sys(0..n);
    wrk : Multprec_Complex_VecVecs.VecVec(1..n);
    timer : Timing_Widget;

    use Multprec_Gradient_Evaluations;

  begin
    put_line("The polynomial : "); put(p);
    Split_Common_Factors(s,f,b);
    tstart(timer);
    for i in 1..m loop
      z := Gradient_of_Polynomial(f,b,c,x);
      Multprec_Complex_Vectors.Clear(z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient function in reverse mode");
    for i in b'range loop
      wrk(i) := new Multprec_Complex_Vectors.Vector(0..n);
    end loop;
    tstart(timer);
    for i in 1..m loop
      Gradient_of_Polynomial(f,b,c,x,wrk,z);
      Multprec_Complex_Vectors.Clear(z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient procedure in reverse mode");
    for k in 1..n loop
      g(k) := Multprec_Complex_Polynomials.Diff(p,k);
    end loop;
    tstart(timer);
    for i in 1..m loop
      z := Multprec_Complex_Poly_SysFun.Eval(g,x);
      Multprec_Complex_Vectors.Clear(z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient plain mode");
    eg := Multprec_Complex_Poly_SysFun.Create(g);
    tstart(timer);
    for i in 1..m loop
      z := Multprec_Complex_Poly_SysFun.Eval(eg,x);
      Multprec_Complex_Vectors.Clear(z);
    end loop;
    tstop(timer);
    new_line;
    print_times(standard_output,timer,"gradient with Horner");
    Standard_Natural_VecVecs.Clear(s);
    Multprec_Complex_Polynomials.Clear(p);
    Multprec_Complex_Poly_Systems.Clear(g);
    Multprec_Complex_Poly_SysFun.Clear(eg);
  end Multprec_Performance_Test;

  procedure Performance_Test is

  -- DESCRIPTION :
  --   We take the next to last equation of the cyclic n-roots problem
  --   and compute its gradient once evaluating seperate derivatives
  --   and once via the reverse mode, for various level of precision.

    n,m : integer32 := 0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the frequency : "); get(m);

    put("double, double double, quad double, or multiprecision? (s/d/q/m) ");
    Ask_Alternative(ans,"sdqm");
    case ans is
      when 's' => Standard_Performance_Test(n,m);
      when 'd' => DoblDobl_Performance_Test(n,m);
      when 'q' => QuadDobl_Performance_Test(n,m);
      when 'm' =>
        declare
          deci,size : natural32 := 0;
        begin
          new_line;
          put("Give the number of decimal places : "); get(deci);
          size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
          Multprec_Performance_Test(n,m,size);
        end;
      when others => null;
    end case;
  end Performance_Test;

  procedure Compare ( fz : in Standard_Complex_Numbers.Complex_Number;
                      gz : in Standard_Complex_Vectors.Vector;
                      fgz : in Standard_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the evaluated polynomial in fz with fgz(0)
  --   and the values in gz with fgz(1..fgz'last).
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use Standard_Complex_Numbers;

    dff : Complex_Number;
    val : double_float;
 
  begin
    dff := fz - fgz(0);
    val := AbsVal(dff);
    if val > tol or output then
      put_line("difference in evaluated polynomial :");
      put(fz); new_line;
      put(fgz(0)); new_line;
      put("error :"); put(val,3); new_line;
    end if;
    for i in gz'range loop
      dff := gz(i) - fgz(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference in gradient at "); put(i,1); put_line(" :");
        put(gz(i)); new_line;
        put(fgz(i)); new_line;
        put("error :"); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Compare ( fz : in DoblDobl_Complex_Numbers.Complex_Number;
                      gz : in DoblDobl_Complex_Vectors.Vector;
                      fgz : in DoblDobl_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the evaluated polynomial in fz with fgz(0)
  --   and the values in gz with fgz(1..fgz'last).
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use DoblDobl_Complex_Numbers;

    dff : Complex_Number;
    val : double_double;
 
  begin
    dff := fz - fgz(0);
    val := AbsVal(dff);
    if val > tol or output then
      put_line("difference in evaluated polynomial :");
      put(fz); new_line;
      put(fgz(0)); new_line;
      put("error : "); put(val,3); new_line;
    end if;
    for i in gz'range loop
      dff := gz(i) - fgz(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference in gradient at "); put(i,1); put_line(" :");
        put(gz(i)); new_line;
        put(fgz(i)); new_line;
        put("error : "); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Compare ( fz : in QuadDobl_Complex_Numbers.Complex_Number;
                      gz : in QuadDobl_Complex_Vectors.Vector;
                      fgz : in QuadDobl_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the evaluated polynomial in fz with fgz(0)
  --   and the values in gz with fgz(1..fgz'last).
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use QuadDobl_Complex_Numbers;

    dff : Complex_Number;
    val : quad_double;
 
  begin
    dff := fz - fgz(0);
    val := AbsVal(dff);
    if val > tol or output then
      put_line("difference in evaluated polynomial :");
      put(fz); new_line;
      put(fgz(0)); new_line;
      put("error : "); put(val,3); new_line;
    end if;
    for i in gz'range loop
      dff := gz(i) - fgz(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference in gradient at "); put(i,1); put_line(" :");
        put(gz(i)); new_line;
        put(fgz(i)); new_line;
        put("error : "); put(val,3); new_line;
      end if;
    end loop;
  end Compare;

  procedure Compare ( fz : in Multprec_Complex_Numbers.Complex_Number;
                      gz : in Multprec_Complex_Vectors.Vector;
                      fgz : in Multprec_Complex_Vectors.Vector;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the evaluated polynomial in fz with fgz(0)
  --   and the values in gz with fgz(1..fgz'last).
  --   If the differences between the corresponding values is larger
  --   than the tolerance tol, then an error message is written.
  --   If the flag output is true, then all differences are written.

    use Multprec_Complex_Numbers;

    dff : Complex_Number;
    val : Floating_Number;
 
  begin
    dff := fz - fgz(0);
    val := AbsVal(dff);
    if val > tol or output then
      put_line("difference in evaluated polynomial :");
      put(fz); new_line;
      put(fgz(0)); new_line;
      put("error : "); put(val,3); new_line;
    end if;
    Clear(dff); Clear(val);
    for i in gz'range loop
      dff := gz(i) - fgz(i);
      val := AbsVal(dff);
      if val > tol or output then
        put("difference in gradient at "); put(i,1); put_line(" :");
        put(gz(i)); new_line;
        put(fgz(i)); new_line;
        put("error : "); put(val,3); new_line;
      end if;
      Clear(dff); Clear(val);
    end loop;
  end Compare;

  procedure Standard_Conditioning_Test
              ( n,d,m,c : natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Generates a polynomial, a point, and its gradient of
  --   prescribed conditioning number, using double arithmetic.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : Standard_Complex_Polynomials.Poly;
    x : Standard_Complex_Vectors.Vector(1..integer32(n));
    g : Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    fgz : Standard_Complex_Vectors.Vector(0..x'last);
    fz : Standard_Complex_Numbers.Complex_Number;
    nt : integer32;
    numfz,denfz,fzrco,gznrc,gzdrc,gzrco : double_float;

    use Standard_Gradient_Evaluations;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,cffsz,pntsz,close,g,p,x);
    fz := Standard_Complex_Poly_Functions.Eval(p,x);
    nt := integer32(Standard_Complex_Polynomials.Number_of_Terms(p));
    declare
      c : Standard_Complex_Vectors.Vector(1..nt);
      wrk : Standard_Complex_VecVecs.VecVec(1..nt);
      e,f,b : Standard_Natural_VecVecs.VecVec(1..nt);
    begin
      Coefficients_and_Supports(p,c,e);
      Split_Common_Factors(e,f,b);
      for i in b'range loop
        wrk(i) := new Standard_Complex_Vectors.Vector(0..integer32(n));
      end loop;
      Gradient_with_Inverse_Condition
        (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,gznrc,gzdrc,gzrco);
      put_line("Comparing evaluated values ...");
      Compare(fz,g,fgz,1.0E-8,true);
      put_line("Condition number of evaluation :");
      put("    numerator : "); put(numfz,3); new_line;
      put("  denominator : "); put(denfz,3); new_line;
      put("  inverse condition number : "); put(fzrco,3); new_line;
      put_line("Condition of gradient :");
      put("    numerator : "); put(gznrc,3); new_line;
      put("  denominator : "); put(gzdrc,3); new_line;
      put("  inverse condition number : "); put(gzrco,3); new_line;
    end;
  end Standard_Conditioning_Test;

  procedure DoblDobl_Conditioning_Test
              ( n,d,m,c : natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Generates a polynomial, a point, and its gradient of
  --   prescribed conditioning number, using double double arithmetic.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : DoblDobl_Complex_Polynomials.Poly;
    x : DoblDobl_Complex_Vectors.Vector(1..integer32(n));
    g : DoblDobl_Complex_Vectors.Vector(1..integer32(n))
      := DoblDobl_Random_Vectors.Random_Vector(1,integer32(n));
    fz : DoblDobl_Complex_Numbers.Complex_Number;
    fgz : DoblDobl_Complex_Vectors.Vector(0..x'last);
    nt : integer32;
    numfz,denfz,fzrco,gznrc,gzdrc,gzrco : double_double;

    use DoblDobl_Gradient_Evaluations;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,cffsz,pntsz,close,g,p,x);
    fz := DoblDobl_Complex_Poly_Functions.Eval(p,x);
    nt := integer32(DoblDobl_Complex_Polynomials.Number_of_Terms(p));
    declare
      c : DoblDobl_Complex_Vectors.Vector(1..nt);
      wrk : DoblDobl_Complex_VecVecs.VecVec(1..nt);
      e,f,b : Standard_Natural_VecVecs.VecVec(1..nt);
    begin
      Coefficients_and_Supports(p,c,e);
      Split_Common_Factors(e,f,b);
      for i in b'range loop
        wrk(i) := new DoblDobl_Complex_Vectors.Vector(0..integer32(n));
      end loop;
      Gradient_with_Inverse_Condition
        (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,gznrc,gzdrc,gzrco);
      put_line("Comparing evaluated values ...");
      Compare(fz,g,fgz,1.0E-8,true);
      put_line("Condition number of evaluation :");
      put("    numerator : "); put(numfz,3); new_line;
      put("  denominator : "); put(denfz,3); new_line;
      put("  inverse condition number : "); put(fzrco,3); new_line;
      put_line("Condition of gradient :");
      put("    numerator : "); put(gznrc,3); new_line;
      put("  denominator : "); put(gzdrc,3); new_line;
      put("  inverse condition number : "); put(gzrco,3); new_line;
    end;
  end DoblDobl_Conditioning_Test;

  procedure QuadDobl_Conditioning_Test
              ( n,d,m,c : natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Generates a polynomial, a point, and its gradient of
  --   prescribed conditioning number, using quad double arithmetic.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : QuadDobl_Complex_Polynomials.Poly;
    x : QuadDobl_Complex_Vectors.Vector(1..integer32(n));
    g : QuadDobl_Complex_Vectors.Vector(1..integer32(n))
      := QuadDobl_Random_Vectors.Random_Vector(1,integer32(n));
    fz : QuadDobl_Complex_Numbers.Complex_Number;
    fgz : QuadDobl_Complex_Vectors.Vector(0..x'last);
    nt : integer32;
    numfz,denfz,fzrco,gznrc,gzdrc,gzrco : quad_double;

    use QuadDobl_Gradient_Evaluations;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,cffsz,pntsz,close,g,p,x);
    fz := QuadDobl_Complex_Poly_Functions.Eval(p,x);
    nt := integer32(QuadDobl_Complex_Polynomials.Number_of_Terms(p));
    declare
      c : QuadDobl_Complex_Vectors.Vector(1..nt);
      wrk : QuadDobl_Complex_VecVecs.VecVec(1..nt);
      e,f,b : Standard_Natural_VecVecs.VecVec(1..nt);
    begin
      Coefficients_and_Supports(p,c,e);
      Split_Common_Factors(e,f,b);
      for i in b'range loop
        wrk(i) := new QuadDobl_Complex_Vectors.Vector(0..integer32(n));
      end loop;
      Gradient_with_Inverse_Condition
        (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,gznrc,gzdrc,gzrco);
      put_line("Comparing evaluated values ...");
      Compare(fz,g,fgz,1.0E-8,true);
      put_line("Condition number of evaluation :");
      put("  inverse condition number : "); put(fzrco,3); new_line;
      put("    numerator : "); put(numfz,3); new_line;
      put("  denominator : "); put(denfz,3); new_line;
      put_line("Condition of gradient :");
      put("    numerator : "); put(gznrc,3); new_line;
      put("  denominator : "); put(gzdrc,3); new_line;
      put("  inverse condition number : "); put(gzrco,3); new_line;
    end;
  end QuadDobl_Conditioning_Test;

  procedure Multprec_Conditioning_Test
              ( n,d,m,c,sz : natural32;
                cffsz,pntsz,close : in double_float ) is

  -- DESCRIPTION :
  --   Generates a polynomial, a point, and its gradient of
  --   prescribed conditioning number.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   size     the size of the numbers;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root.

    p : Multprec_Complex_Polynomials.Poly;
    x : Multprec_Complex_Vectors.Vector(1..integer32(n));
    g : Multprec_Complex_Vectors.Vector(1..integer32(n))
      := Multprec_Random_Vectors.Random_Vector(1,integer32(n),sz);
    fz : Multprec_Complex_Numbers.Complex_Number;
    fgz : Multprec_Complex_Vectors.Vector(0..x'last);
    nt : integer32;
    numfz,denfz,fzrco,gznrc,gzdrc,gzrco : Floating_Number;

    use Multprec_Gradient_Evaluations;

  begin
    Random_Conditioned_Gradient_Evaluation(n,d,m,c,sz,cffsz,pntsz,close,g,p,x);
    nt := integer32(Multprec_Complex_Polynomials.Number_of_Terms(p));
    declare
      c : Multprec_Complex_Vectors.Vector(1..nt);
      wrk : Multprec_Complex_VecVecs.VecVec(1..nt);
      e,f,b : Standard_Natural_VecVecs.VecVec(1..nt);
    begin
      Coefficients_and_Supports(p,c,e);
      Split_Common_Factors(e,f,b);
      for i in b'range loop
        wrk(i) := new Multprec_Complex_Vectors.Vector(0..integer32(n));
      end loop;
      Gradient_with_Inverse_Condition
        (f,b,c,x,wrk,fgz,numfz,denfz,fzrco,gznrc,gzdrc,gzrco);
      put_line("Comparing evaluated values ...");
      Compare(fz,g,fgz,1.0E-8,true);
      put_line("Condition number of evaluation :");
      put("  inverse condition number : "); put(fzrco,3); new_line;
      put("    numerator : "); put(numfz,3); new_line;
      put("  denominator : "); put(denfz,3); new_line;
      put_line("Condition of gradient :");
      put("    numerator : "); put(gznrc,3); new_line;
      put("  denominator : "); put(gzdrc,3); new_line;
      put("  inverse condition number : "); put(gzrco,3); new_line;
    end;
  end Multprec_Conditioning_Test;

  procedure Ask_for_Parameters
              ( n,d,m,c : out natural32;
                cff,pnt,cls : out double_float;
                prc : out character ) is

  -- DESCRIPTIN :
  --   Prompts the user for the parameters of the problem.

  -- ON RETURN :
  --   n        the dimension is the number of variables;
  --   d        the largest degree of the monomials;
  --   m        the number of monomials (0 for dense);
  --   c        coefficient type, complex, one, or real (1, 2, 3);
  --   cff      magnitude of the coefficients;
  --   pnt      magnitude of the points;
  --   cls      closeness to a root;
  --   prc      working precision, either s, d, q, or m for double,
  --            double double, quad double, or arbitrary multiprecision.

    cond : double_float;

  begin
    n := 0; d := 0; m := 0; c := 0;
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of monomials (0 for dense) : "); get(m);
    put("Give the type of coefficient (0 cmplx, 1 one, 2 real) : ");
    get(c);
    cff := 0.0; pnt := 0.0; cls := 0.0;
    put("Give magnitude of the coefficients : "); get(cff);
    put("Give magnitude of the coordinates of the point : "); get(pnt);
    put("Give closeness to a root : "); get(cls);
    cond := cff*(pnt**integer(d))/cls;
    new_line;
    put("Predicted condition number : "); put(cond,3); new_line;
    new_line;
    put("double, double double, quad double, or multiprecision? (s/d/q/m) ");
    Ask_Alternative(prc,"sdqm");
  end Ask_for_Parameters;

  procedure Gradient_Conditioning_Test is

  -- DESCRIPTION :
  --   For various levels of precision, a polynomial f is generated
  --   and a point z with prescribed condition number.

    n,d,m,c,deci,size : natural32 := 0;
    cff,pnt,cls,cond : double_float := 0.0;
    ans : character;

  begin
    new_line;
    put("Give the dimension : "); get(n);
    put("Give the largest degree : "); get(d);
    put("Give the number of monomials (0 for dense) : "); get(m);
    put("Give the type of coefficient (0 cmplx, 1 one, 2, real) : ");
    get(c);
    put("Give magnitude of the coefficients : "); get(cff);
    put("Give magnitude of the coordinates of the point : "); get(pnt);
    put("Give closeness to a root : "); get(cls);
    cond := cff*(pnt**integer(d))/cls;
    new_line;
    put("Predicted condition number : "); put(cond,3); new_line;
    new_line;
    put("double, double double, quad double, or multiprecision? (s/d/q/m) ");
    Ask_Alternative(ans,"sdqm");
    new_line;
    case ans is
      when 's' => Standard_Conditioning_Test(n,d,m,c,cff,pnt,cls);
      when 'd' => DoblDobl_Conditioning_Test(n,d,m,c,cff,pnt,cls);
      when 'q' => QuadDobl_Conditioning_Test(n,d,m,c,cff,pnt,cls);
      when 'm' =>
        put("Give the number of decimal places : "); get(deci);
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Multprec_Conditioning_Test(n,d,m,c,size,cff,pnt,cls);
      when others => null;
    end case;
  end Gradient_Conditioning_Test;

  procedure Compare ( A : in Standard_Complex_Matrices.Matrix;
                      B : in Standard_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in the Jacobian matrix A
  --   with the evaluated gradients in B.
  --   If the corresponding values in A and B differ by more than
  --   the tolerance tol, then an error message is written.
  --   If the flag output is true, then all numbers are written.

    use Standard_Complex_Numbers;

    dff : Complex_Number;
    val : double_float;

  begin
    for col in A'range(2) loop
      for row in A'range(1) loop
        dff := A(row,col) - B(row)(col);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference in Jacobian and gradient at");
          put(" row = "); put(row,1);
          put(", column = "); put(col,1); put_line(" :");
          put(A(row,col)); new_line;
          put(B(row)(col)); new_line;
          put("error :"); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( A : in DoblDobl_Complex_Matrices.Matrix;
                      B : in DoblDobl_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in the Jacobian matrix A
  --   with the evaluated gradients in B.
  --   If the corresponding values in A and B differ by more than
  --   the tolerance tol, then an error message is written.
  --   If the flag output is true, then all numbers are written.

    use DoblDobl_Complex_Numbers;

    dff : Complex_Number;
    val : double_double;

  begin
    for col in A'range(2) loop
      for row in A'range(1) loop
        dff := A(row,col) - B(row)(col);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference in Jacobian and gradient at");
          put(" row = "); put(row,1);
          put(", column = "); put(col,1); put_line(" :");
          put(A(row,col)); new_line;
          put(B(row)(col)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( A : in QuadDobl_Complex_Matrices.Matrix;
                      B : in QuadDobl_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in the Jacobian matrix A
  --   with the evaluated gradients in B.
  --   If the corresponding values in A and B differ by more than
  --   the tolerance tol, then an error message is written.
  --   If the flag output is true, then all numbers are written.

    use QuadDobl_Complex_Numbers;

    dff : Complex_Number;
    val : quad_double;

  begin
    for col in A'range(2) loop
      for row in A'range(1) loop
        dff := A(row,col) - B(row)(col);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference in Jacobian and gradient at");
          put(" row = "); put(row,1);
          put(", column = "); put(col,1); put_line(" :");
          put(A(row,col)); new_line;
          put(B(row)(col)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
      end loop;
    end loop;
  end Compare;

  procedure Compare ( A : in Multprec_Complex_Matrices.Matrix;
                      B : in Multprec_Complex_VecVecs.VecVec;
                      tol : in double_float; output : in boolean ) is

  -- DESCRIPTION :
  --   Compares the values in the Jacobian matrix A
  --   with the evaluated gradients in B.
  --   If the corresponding values in A and B differ by more than
  --   the tolerance tol, then an error message is written.
  --   If the flag output is true, then all numbers are written.

    use Multprec_Complex_Numbers;

    dff : Complex_Number;
    val : Floating_Number;

  begin
    for col in A'range(2) loop
      for row in A'range(1) loop
        dff := A(row,col) - B(row)(col);
        val := AbsVal(dff);
        if val > tol or output then
          put("difference in Jacobian and gradient at");
          put(" row = "); put(row,1);
          put(", column = "); put(col,1); put_line(" :");
          put(A(row,col)); new_line;
          put(B(row)(col)); new_line;
          put("error : "); put(val,3); new_line;
        end if;
        Clear(dff); Clear(val);
      end loop;
    end loop;
  end Compare;

  function Jacobian ( p : Standard_Complex_Poly_Systems.Poly_Sys;
                      x : Standard_Complex_Vectors.Vector )
                    return Standard_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the straightforward computation and evaluation
  --   of the Jacobian matrix of the system p at the point x,
  --   computed in standard double precision, for testing purposes.

    res : Standard_Complex_Matrices.Matrix(p'range,x'range);
    jac : Standard_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
        := Standard_Complex_Jaco_Matrices.Create(p);

  begin
    res := Standard_Complex_Jaco_Matrices.Eval(jac,x);
    Standard_Complex_Jaco_Matrices.Clear(jac);
    return res;
  end Jacobian;

  function Jacobian ( p : DoblDobl_Complex_Poly_Systems.Poly_Sys;
                      x : DoblDobl_Complex_Vectors.Vector )
                    return DoblDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the straightforward computation and evaluation
  --   of the Jacobian matrix of the system p at the point x,
  --   computed in DoblDobl double precision, for testing purposes.

    res : DoblDobl_Complex_Matrices.Matrix(p'range,x'range);
    jac : DoblDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
        := DoblDobl_Complex_Jaco_Matrices.Create(p);

  begin
    res := DoblDobl_Complex_Jaco_Matrices.Eval(jac,x);
    DoblDobl_Complex_Jaco_Matrices.Clear(jac);
    return res;
  end Jacobian;

  function Jacobian ( p : QuadDobl_Complex_Poly_Systems.Poly_Sys;
                      x : QuadDobl_Complex_Vectors.Vector )
                    return QuadDobl_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the straightforward computation and evaluation
  --   of the Jacobian matrix of the system p at the point x,
  --   computed in QuadDobl double precision, for testing purposes.

    res : QuadDobl_Complex_Matrices.Matrix(p'range,x'range);
    jac : QuadDobl_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
        := QuadDobl_Complex_Jaco_Matrices.Create(p);

  begin
    res := QuadDobl_Complex_Jaco_Matrices.Eval(jac,x);
    QuadDobl_Complex_Jaco_Matrices.Clear(jac);
    return res;
  end Jacobian;

  function Jacobian ( p : Multprec_Complex_Poly_Systems.Poly_Sys;
                      x : Multprec_Complex_Vectors.Vector )
                    return Multprec_Complex_Matrices.Matrix is

  -- DESCRIPTION :
  --   Returns the straightforward computation and evaluation
  --   of the Jacobian matrix of the system p at the point x,
  --   computed in Multprec double precision, for testing purposes.

    res : Multprec_Complex_Matrices.Matrix(p'range,x'range);
    jac : Multprec_Complex_Jaco_Matrices.Jaco_Mat(p'range,x'range)
        := Multprec_Complex_Jaco_Matrices.Create(p);

  begin
    res := Multprec_Complex_Jaco_Matrices.Eval(jac,x);
    Multprec_Complex_Jaco_Matrices.Clear(jac);
    return res;
  end Jacobian;
 
  procedure Standard_Jacobian_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with standard complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+16.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : Standard_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Standard_Complex_Vectors.Vector(p'range);
    jm : Standard_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    jmpx : Standard_Complex_Matrices.Matrix(p'range,x'range);
    rco,fxnrc,fxdrc,fxrco,maxng,mindg,rcogd : double_float;
    loss_jac,loss_eva : integer32;

    use Standard_Gradient_Evaluations;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if rco = 0.0
     then loss_eva := -2**30;
     else loss_eva := integer32(log10(rco));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    declare
      cff : Standard_Complex_VecVecs.VecVec(p'range);
      ydx : Standard_Complex_VecVecs.VecVec(p'range);
      e,f,b : Standard_Natural_VecVecs.Array_of_VecVecs(p'range);
      wrk : Standard_Complex_VecVecs.Array_of_VecVecs(p'range);
      nt : integer32;
    begin
      for i in p'range loop
        nt := integer32(Standard_Complex_Polynomials.Number_of_Terms(p(i)));
        cff(i) := new Standard_Complex_Vectors.Vector(1..nt);
        e(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Coefficients_and_Supports(p(i),cff(i).all,e(i).all);
        f(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        b(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Split_Common_Factors(e(i).all,f(i).all,b(i).all);
        wrk(i) := new Standard_Complex_VecVecs.VecVec(1..nt);
        for j in b(i)'range loop
          wrk(i)(j) := new Standard_Complex_Vectors.Vector(0..x'last);
        end loop;
        ydx(i) := new Standard_Complex_Vectors.Vector(0..x'last);
      end loop;
      Jacobian_with_Inverse_Condition
        (f,b,cff,x,wrk,ydx,fxnrc,fxdrc,fxrco,maxng,mindg,rcogd);
      new_line;
      put_line("Comparing values in Jacobian matrix with gradients ...");
      Compare(jm,ydx,1.0E-8,true);
      new_line;
      put_line("Comparing values with straightforward computation ...");
      jmpx := Jacobian(p,x);
      Compare(jmpx,ydx,1.0E-8,true);
      put_line("Condition of the polynomial evaluation :");
      put("  numerator : "); put(fxnrc,3); new_line;
      put("denominator : "); put(fxdrc,3); new_line;
      put("inverse condition number : "); put(fxrco,3); new_line;
      put_line("Condition of the Jacobian evaluation :");
      put("  numerator : "); put(maxng,3); new_line;
      put("denominator : "); put(mindg,3); new_line;
      put("inverse condition number : "); put(rcogd,3); new_line;
    end;
  end Standard_Jacobian_Test;

  procedure DoblDobl_Jacobian_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with double double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+32.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : DoblDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : DoblDobl_Complex_Vectors.Vector(p'range);
    jm : DoblDobl_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    jmpx : DoblDobl_Complex_Matrices.Matrix(p'range,x'range);
    rco,fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco : double_double;
    loss_jac,loss_eva : integer32;

    use DoblDobl_Gradient_Evaluations;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Is_Zero(rco)
     then loss_eva := -2**30;
     else loss_eva := integer32(to_double(log10(rco)));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    declare
      cff : DoblDobl_Complex_VecVecs.VecVec(p'range);
      ydx : DoblDobl_Complex_VecVecs.VecVec(p'range);
      e,f,b : Standard_Natural_VecVecs.Array_of_VecVecs(p'range);
      wrk : DoblDobl_Complex_VecVecs.Array_of_VecVecs(p'range);
      nt : integer32;
    begin
      for i in p'range loop
        nt := integer32(DoblDobl_Complex_Polynomials.Number_of_Terms(p(i)));
        cff(i) := new DoblDobl_Complex_Vectors.Vector(1..nt);
        e(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Coefficients_and_Supports(p(i),cff(i).all,e(i).all);
        f(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        b(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Split_Common_Factors(e(i).all,f(i).all,b(i).all);
        wrk(i) := new DoblDobl_Complex_VecVecs.VecVec(1..nt);
        for j in b(i)'range loop
          wrk(i)(j) := new DoblDobl_Complex_Vectors.Vector(0..x'last);
        end loop;
        ydx(i) := new DoblDobl_Complex_Vectors.Vector(0..x'last);
      end loop;
      Jacobian_with_Inverse_Condition
        (f,b,cff,x,wrk,ydx,fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco);
      new_line;
      put_line("Comparing values in Jacobian matrix with gradients ...");
      Compare(jm,ydx,1.0E-8,true);
      new_line;
      put_line("Comparing values with straightforward computation ...");
      jmpx := Jacobian(p,x);
      Compare(jmpx,ydx,1.0E-8,true);
      put_line("Condition of the polynomial evaluation :");
      put("  numerator : "); put(fxnrc,3); new_line;
      put("denominator : "); put(fxdrc,3); new_line;
      put("inverse condition number : "); put(fxrco,3); new_line;
      put_line("Condition of the Jacobian evaluation :");
      put("  numerator : "); put(gxnrc,3); new_line;
      put("denominator : "); put(gxdrc,3); new_line;
      put("inverse condition number : "); put(gxrco,3); new_line;
    end;
  end DoblDobl_Jacobian_Test;

  procedure QuadDobl_Jacobian_Test
              ( n,d,m,c : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with quad double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided condjm < 1.0E+64.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : QuadDobl_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : QuadDobl_Complex_Vectors.Vector(p'range);
    jm : QuadDobl_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    jmpx : QuadDobl_Complex_Matrices.Matrix(p'range,x'range);
    rco,fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco : quad_double;
    loss_jac,loss_eva : integer32;

    use QuadDobl_Gradient_Evaluations;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation(n,d,m,c,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Is_Zero(rco)
     then loss_eva := -2**30;
     else loss_eva := integer32(to_double(log10(rco)));
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    declare
      cff : QuadDobl_Complex_VecVecs.VecVec(p'range);
      ydx : QuadDobl_Complex_VecVecs.VecVec(p'range);
      e,f,b : Standard_Natural_VecVecs.Array_of_VecVecs(p'range);
      wrk : QuadDobl_Complex_VecVecs.Array_of_VecVecs(p'range);
      nt : integer32;
    begin
      for i in p'range loop
        nt := integer32(QuadDobl_Complex_Polynomials.Number_of_Terms(p(i)));
        cff(i) := new QuadDobl_Complex_Vectors.Vector(1..nt);
        e(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Coefficients_and_Supports(p(i),cff(i).all,e(i).all);
        f(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        b(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Split_Common_Factors(e(i).all,f(i).all,b(i).all);
        wrk(i) := new QuadDobl_Complex_VecVecs.VecVec(1..nt);
        for j in b(i)'range loop
          wrk(i)(j) := new QuadDobl_Complex_Vectors.Vector(0..x'last);
        end loop;
        ydx(i) := new QuadDobl_Complex_Vectors.Vector(0..x'last);
      end loop;
      Jacobian_with_Inverse_Condition
        (f,b,cff,x,wrk,ydx,fxnrc,fxdrc,fxrco,gxnrc,gxdrc,gxrco);
      new_line;
      put_line("Comparing values in Jacobian matrix with gradients ...");
      Compare(jm,ydx,1.0E-8,true);
      new_line;
      put_line("Comparing values with straightforward computation ...");
      jmpx := Jacobian(p,x);
      Compare(jmpx,ydx,1.0E-8,true);
      put_line("Condition of the polynomial evaluation :");
      put("  numerator : "); put(fxnrc,3); new_line;
      put("denominator : "); put(fxdrc,3); new_line;
      put("inverse condition number : "); put(fxrco,3); new_line;
      put_line("Condition of the Jacobian evaluation :");
      put("  numerator : "); put(gxnrc,3); new_line;
      put("denominator : "); put(gxdrc,3); new_line;
      put("inverse condition number : "); put(gxrco,3); new_line;
    end;
  end QuadDobl_Jacobian_Test;

  procedure Multprec_Jacobian_Test
              ( n,d,m,c,size : in natural32;
                cffsz,pntsz,close,condjm : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial system with quad double complex coefficients.
  --   The condition number of the numerical evaluation problem is
  --   determined by the degree, the coefficient size, the size of the
  --   coordinates of the point, and the distance of the point to the
  --   closest root.  The condition number of the Jacobian number
  --   will be as prescribed provided the value for the size of
  --   the working precision is sufficiently high.

  -- ON ENTRY :
  --   n        number of variables;
  --   d        largest degree of the monomials;
  --   m        number of monomials (0 for a dense polynomial);
  --   c        type of coefficient, 0 is random complex, 1 is one,
  --            and 2 is random real;
  --   size     size of the working precision;
  --   cffsz    size of the coefficients;
  --   pntsz    size of the coordinates of the point where to evaluate;
  --   close    distance of the point to a root;
  --   condjm   condition number of the Jacobian matrix.

    p : Multprec_Complex_Poly_Systems.Poly_Sys(1..integer32(n));
    x : Multprec_Complex_Vectors.Vector(p'range);
    jm : Multprec_Complex_Matrices.Matrix(p'range,x'range)
       := Random_Conditioned_Matrix(integer32(n),condjm);
    jmpx : Multprec_Complex_Matrices.Matrix(p'range,x'range);
    rco,fxnrc,fxdrc,fxrco,maxng,mindg,rcogd : Floating_Number;
    loss_jac,loss_eva : integer32;

    use Multprec_Gradient_Evaluations;

  begin
    put("The given condition number of the Jacobian :");
    put(condjm,3); new_line;
    loss_jac := Estimated_Loss_of_Decimal_Places(jm);
    put("-> Estimated loss of decimal places : "); put(loss_jac,1); new_line;
    Random_Conditioned_Jacobian_Evaluation
      (n,d,m,c,size,cffsz,pntsz,close,jm,p,x);
    rco := Inverse_Condition_Number(p,x);
    if Equal(rco,0.0) then
      loss_eva := -2**30;
    else
      declare
        mp_log10rco : Floating_Number := log10(rco);
        st_log10rco : double_float := Round(mp_log10rco);
      begin
        loss_eva := integer32(st_log10rco);
      end;
    end if;
    new_line;
    put("Inverse condition of polynomial evaluation :");
    put(rco,3); new_line;
    put("-> Estimated loss of decimal places : "); put(loss_eva,1); new_line;
    declare
      cff : Multprec_Complex_VecVecs.VecVec(p'range);
      ydx : Multprec_Complex_VecVecs.VecVec(p'range);
      e,f,b : Standard_Natural_VecVecs.Array_of_VecVecs(p'range);
      wrk : Multprec_Complex_VecVecs.Array_of_VecVecs(p'range);
      nt : integer32;
    begin
      for i in p'range loop
        nt := integer32(Multprec_Complex_Polynomials.Number_of_Terms(p(i)));
        cff(i) := new Multprec_Complex_Vectors.Vector(1..nt);
        e(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Coefficients_and_Supports(p(i),cff(i).all,e(i).all);
        f(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        b(i) := new Standard_Natural_VecVecs.VecVec(1..nt);
        Split_Common_Factors(e(i).all,f(i).all,b(i).all);
        wrk(i) := new Multprec_Complex_VecVecs.VecVec(1..nt);
        for j in b(i)'range loop
          wrk(i)(j) := new Multprec_Complex_Vectors.Vector(0..x'last);
        end loop;
        ydx(i) := new Multprec_Complex_Vectors.Vector(0..x'last);
      end loop;
      Jacobian_with_Inverse_Condition
        (f,b,cff,x,wrk,ydx,fxnrc,fxdrc,fxrco,maxng,mindg,rcogd);
      new_line;
      put_line("Comparing values with those in the gradients ...");
      Compare(jm,ydx,1.0E-8,true);
      new_line;
      put_line("Comparing values with straightforward computation ...");
      jmpx := Jacobian(p,x);
      Compare(jmpx,ydx,1.0E-8,true);
      put_line("Condition of the polynomial evaluation :");
      put("  numerator : "); put(fxnrc,3); new_line;
      put("denominator : "); put(fxdrc,3); new_line;
      put("inverse condition number : "); put(fxrco,3); new_line;
      put_line("Condition of the Jacobian evaluation :");
      put("  numerator : "); put(maxng,3); new_line;
      put("denominator : "); put(mindg,3); new_line;
      put("inverse condition number : "); put(rcogd,3); new_line;
    end;
  end Multprec_Jacobian_Test;

  procedure Jacobian_Conditioning_Test is

    n,d,m,c,deci,size : natural32 := 0;
    cff,pnt,cls,cndjm : double_float := 0.0;
    ans : character;

  begin
    Ask_for_Parameters(n,d,m,c,cff,pnt,cls,ans);
    new_line;
    put("Condition number of the Jacobian matrix : "); get(cndjm);
    new_line;
    case ans is
      when 's' => Standard_Jacobian_Test(n,d,m,c,cff,pnt,cls,cndjm);
      when 'd' => DoblDobl_Jacobian_Test(n,d,m,c,cff,pnt,cls,cndjm);
      when 'q' => QuadDobl_Jacobian_Test(n,d,m,c,cff,pnt,cls,cndjm);
      when 'm' =>
        put("Give the number of decimal places : "); get(deci);
        size := Multprec_Floating_Numbers.Decimal_to_Size(deci);
        Multprec_Jacobian_Test(n,d,m,c,size,cff,pnt,cls,cndjm);
      when others => null;
    end case;
  end Jacobian_Conditioning_Test;

  procedure Main is

    ans : character;

  begin
    new_line;
    put_line("Evaluation of a gradient and its condition number ...");
    put_line("  1. generate a random polynomial and select the precision;");
    put_line("  2. do a performance test on a cyclic n-roots polynomial;");
    put_line("  3. generate a polynomial with prescribed conditioning;");
    put_line("  4. generate a system with prescribed conditioning.");
    put("Type 1, 2, 3, or 4 to choose : ");
    Ask_Alternative(ans,"1234");
    case ans is
      when '1' => Evaluate_Gradient;
      when '2' => Performance_Test;
      when '3' => Gradient_Conditioning_Test;
      when '4' => Jacobian_Conditioning_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_vmpdiff;
