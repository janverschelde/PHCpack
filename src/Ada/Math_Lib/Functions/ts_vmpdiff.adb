with text_io;                             use text_io;
with Communications_with_User;            use Communications_with_User;
with Standard_Natural_Numbers;            use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;         use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;            use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;         use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;           use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;        use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;         use Standard_Complex_Numbers_io;
with Standard_Natural_Vectors;
with Standard_Natural_Vectors_io;         use Standard_Natural_Vectors_io;
with Standard_Natural_VecVecs;
with Standard_Integer_Vectors;
with Standard_Floating_Vectors;
with Standard_Floating_Vectors_io;        use Standard_Floating_Vectors_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;         use Standard_Complex_Vectors_io;
with Standard_Complex_VecVecs;
with Standard_Complex_Norms_Equals;       use Standard_Complex_Norms_Equals;
with Standard_Random_Numbers;
with Standard_Random_Vectors;             use Standard_Random_Vectors;
with Standard_Speelpenning_Products;
with Double_Double_Numbers;               use Double_Double_Numbers;
with Double_Double_Numbers_io;            use Double_Double_Numbers_io;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;         use DoblDobl_Complex_Numbers_io;
with Double_Double_Vectors;
with Double_Double_Vectors_io;            use Double_Double_Vectors_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;         use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_VecVecs;
with DoblDobl_Complex_Vector_Norms;       use DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with DoblDobl_Speelpenning_Products;
with Quad_Double_Numbers;                 use Quad_Double_Numbers;
with Quad_Double_Numbers_io;              use Quad_Double_Numbers_io;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;         use QuadDobl_Complex_Numbers_io;
with Quad_Double_Vectors;
with Quad_Double_Vectors_io;              use Quad_Double_Vectors_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;         use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_VecVecs;
with QuadDobl_Complex_Vector_Norms;       use QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with QuadDobl_Speelpenning_Products;
with Multprec_Floating_Numbers;           use Multprec_Floating_Numbers;
with Multprec_Floating_Numbers_io;        use Multprec_Floating_Numbers_io;
with Multprec_Complex_Numbers;
with Multprec_Complex_Numbers_io;         use Multprec_Complex_Numbers_io;
with Multprec_Floating_Vectors;
with Multprec_Floating_Vectors_io;        use Multprec_Floating_Vectors_io;
with Multprec_Complex_Vectors;
with Multprec_Complex_Vectors_io;         use Multprec_Complex_Vectors_io;
with Multprec_Complex_Norms_Equals;       use Multprec_Complex_Norms_Equals;
with Multprec_Random_Vectors;
with Multprec_Complex_VecVecs;
with Multprec_Speelpenning_Products;
with Standard_Monomial_Evaluations;
with DoblDobl_Monomial_Evaluations;
with QuadDobl_Monomial_Evaluations;
with Multprec_Monomial_Evaluations;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Poly_Functions;
with Standard_Complex_Polynomials_io;     use Standard_Complex_Polynomials_io;
with DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Poly_Functions;
with QuadDobl_Complex_Polynomials;
with QuadDobl_Complex_Poly_Functions;
with Multprec_Complex_Polynomials;
with Multprec_Complex_Poly_Functions;
with Coefficient_Supported_Polynomials;   use Coefficient_Supported_Polynomials;
with Standard_Gradient_Evaluations;
with DoblDobl_Gradient_Evaluations;
with QuadDobl_Gradient_Evaluations;
with Multprec_Gradient_Evaluations;
with Varbprec_Polynomial_Evaluations;     use Varbprec_Polynomial_Evaluations;
with Varbprec_Gradient_Evaluations;       use Varbprec_Gradient_Evaluations;

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

  procedure Main is

   -- ans : character;

  begin
    new_line;
    put_line("Evaluation of a gradient and its condition number ...");
   -- put_line("  1. generate a random polynomial and select the precision.");
   -- put("Type 1 to choose : ");
   -- Ask_Alternative(ans,"1");
   -- case ans is
   --   when '1' => 
    Evaluate_Gradient;
   --   when others => null;
   -- end case;
  end Main;

begin
  Main;
end ts_vmpdiff;
