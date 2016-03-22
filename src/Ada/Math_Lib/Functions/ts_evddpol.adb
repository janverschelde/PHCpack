with text_io;                           use text_io;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Symbol_Table;
with DoblDobl_Complex_Numbers;          use DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;          use DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Polynomials;      use DoblDobl_Complex_Polynomials;
with DoblDobl_Complex_Polynomials_io;   use DoblDobl_Complex_Polynomials_io;
with DoblDobl_Complex_Poly_Functions;   use DoblDobl_Complex_Poly_Functions;

procedure ts_evddpol is

-- DESCRIPTION :
--   Test on evaluating a polynomial with double double complex coefficients.

  procedure Test ( n : in natural32; p : in Poly ) is

  -- DESCRIPTION :
  --   Test on the evaluation of the polynomial p.
  --   The number of varables in p equals n.

    x : Vector(1..integer32(n));
    y,z : Complex_Number;
    f : Eval_Poly := Create(p);

  begin
    put("Give "); put(n,1); put_line(" complex numbers :");
    for i in x'range loop
      get(x(i));
    end loop;
    put_line("x ="); put_line(x);
    y := Eval(p,x);
    put_line("y = "); put(y); new_line;
    z := Eval(f,x);
    put_line("z = "); put(z); new_line;
  end Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts the user for the number of variables and a polynomial.
  --   Then calls the evaluation test.

    n : natural32 := 0;
    p : Poly;

  begin
    put("Give the number of variables : "); get(n);
    Symbol_Table.Init(n);
    put("Give a polynomial : "); get(p);
    put("Your polynomial : "); put(p); new_line;
    Test(n,p);
  end Main;

begin
  Main;
end ts_evddpol;
