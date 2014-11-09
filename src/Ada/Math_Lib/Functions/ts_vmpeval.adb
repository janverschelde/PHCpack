with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Standard_Complex_Numbers;          use Standard_Complex_Numbers;
with Standard_Natural_Vectors;
with Standard_Complex_Vectors;
with Standard_Random_Vectors;
with Symbol_Table;
with Standard_Complex_Polynomials;
with Standard_Complex_Polynomials_io;   use Standard_Complex_Polynomials_io;
with Standard_Complex_Poly_Functions;   use Standard_Complex_Poly_Functions;
with Standard_Random_Polynomials;
with VarbPrec_Polynomial_Evaluations;   use VarbPrec_Polynomial_Evaluations;

procedure ts_vmpeval is

-- DESCRIPTION :
--   Tests condition number computation of polynomial evaluation.

  procedure Standard_Test ( n,d,m,c : in natural32; cond : in double_float ) is

  -- DESCRIPTION :
  --   Given in n the number of variables, in d the largest degree,
  --   in m the maximum number of monomials (or 0 for dense), and in c
  --   the type of coefficients, then this procedure will generate a
  --   random polynomial with standard complex coefficients and show it.
  --   Generates coordinates for the a random point and then adjusts the
  --   constant term of the polynomial to match the given condition cond.

    p : Standard_Complex_Polynomials.Poly;
    x : Standard_Complex_Vectors.Vector(1..integer32(n))
      := Standard_Random_Vectors.Random_Vector(1,integer32(n));
    t : Standard_Complex_Polynomials.Term;
    rco : double_float;

  begin
    if m = 0
     then p := Standard_Random_Polynomials.Random_Dense_Poly(n,d,c);
     else p := Standard_Random_Polynomials.Random_Sparse_Poly(n,d,m,c);
    end if;
    put_line("The generated random polynomial :"); put(p); new_line;
    t.cf := Eval(p,x);
    t.dg := new Standard_Natural_Vectors.Vector'(1..integer32(n) => 0);
    Standard_Complex_Polynomials.Sub(p,t);
    t.cf := Create(1.0/cond);
    Standard_Complex_Polynomials.Add(p,t);
    rco := Inverse_Condition_Number(p,x);
    put("The inverse condition number :"); put(rco,3); new_line;
  end Standard_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the dimensions of a random polynomial.

    n,d,m : natural32 := 0;
    cond : double_float := 0.0;

  begin
    new_line;
    put_line("Condition numbers for polynomial evaluation ...");
    new_line;
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give the maximal degree : "); get(d);
    put("Give number of monomials (0 for dense): "); get(m);
    put("Give the condition number : "); get(cond);
    Standard_Test(n,d,m,2,cond);
  end Main;

begin
  Main;
end ts_vmpeval;
