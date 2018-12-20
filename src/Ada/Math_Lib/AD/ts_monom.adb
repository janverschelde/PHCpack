with text_io;                           use text_io;
with Communications_with_User;          use Communications_with_User;
with Standard_Natural_Numbers;          use Standard_Natural_Numbers;
with Standard_Natural_Numbers_io;       use Standard_Natural_Numbers_io;
with Standard_Integer_Numbers;          use Standard_Integer_Numbers;
with Standard_Integer_Numbers_io;       use Standard_Integer_Numbers_io;
with Standard_Floating_Numbers;         use Standard_Floating_Numbers;
with Standard_Floating_Numbers_io;      use Standard_Floating_Numbers_io;
with Double_Double_Numbers;             use Double_Double_Numbers;
with Quad_Double_Numbers;               use Quad_Double_Numbers;
with Standard_Complex_Numbers;
with Standard_Complex_Numbers_io;       use Standard_Complex_Numbers_io;
with Standard_Complex_Vectors;
with Standard_Complex_Vectors_io;       use Standard_Complex_Vectors_io;
with Standard_Complex_Vector_Norms;
with Standard_Random_Vectors;
with DoblDobl_Complex_Numbers;
with DoblDobl_Complex_Numbers_io;       use DoblDobl_Complex_Numbers_io;
with DoblDobl_Complex_Vectors;
with DoblDobl_Complex_Vectors_io;       use DoblDobl_Complex_Vectors_io;
with DoblDobl_Complex_Vector_Norms;
with DoblDobl_Random_Vectors;
with QuadDobl_Complex_Numbers;
with QuadDobl_Complex_Numbers_io;       use QuadDobl_Complex_Numbers_io;
with QuadDobl_Complex_Vectors;
with QuadDobl_Complex_Vectors_io;       use QuadDobl_Complex_Vectors_io;
with QuadDobl_Complex_Vector_Norms;
with QuadDobl_Random_Vectors;
with Standard_Complex_Monomials;
with Standard_Complex_Monomials_io;     use Standard_Complex_Monomials_io;
with DoblDobl_Complex_Monomials;
with DoblDobl_Complex_Monomials_io;     use DoblDobl_Complex_Monomials_io;
with QuadDobl_Complex_Monomials;
with QuadDobl_Complex_Monomials_io;     use QuadDobl_Complex_Monomials_io;
with Random_Monomials;

procedure ts_monom is

-- DESCRIPTION :
--   Tests the operations on monomials.

  procedure Standard_Eval ( m : in Standard_Complex_Monomials.Monomial ) is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in double precision.

    x : constant Standard_Complex_Vectors.Vector(1..integer32(m.dim))
      := Standard_Random_Vectors.Random_Vector(1,integer32(m.dim));
    y,z : Standard_Complex_Numbers.Complex_Number;
    yd,zd : Standard_Complex_Vectors.Vector(x'range)
          := (x'range => Standard_Complex_Numbers.Create(0.0));
    nrm : double_float;

  begin
    z := Standard_Complex_Monomials.Eval(m,x);
    put("z : "); put(z); new_line;
    Standard_Complex_Monomials.Diff(m,x,zd);
    Standard_Complex_Monomials.Speel(m,x,y,yd);
    put("y : "); put(y); new_line;
    Standard_Complex_Numbers.Sub(y,z);
    nrm := Standard_Complex_Numbers.AbsVal(y);
    put("Norm of the difference :"); put(nrm,3); new_line;
    put_line("The derivatives computed with the Speelpenning algorithm :");
    put_line(yd(1..integer32(m.nvr)));
    put_line("The derivatives computed with the straightforward algorithm :");
    put_line(zd(1..integer32(m.nvr)));
    Standard_Complex_Vectors.Sub(yd,zd);
    nrm := Standard_Complex_Vector_Norms.Max_Norm(yd);
    put("Max norm of the difference :"); put(nrm,3); new_line;
  end Standard_Eval;

  procedure DoblDobl_Eval ( m : in DoblDobl_Complex_Monomials.Monomial ) is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in double double precision.

    x : constant DoblDobl_Complex_Vectors.Vector(1..integer32(m.dim))
      := DoblDobl_Random_Vectors.Random_Vector(1,integer32(m.dim));
    y,z : DoblDobl_Complex_Numbers.Complex_Number;
    yd,zd : DoblDobl_Complex_Vectors.Vector(x'range)
          := (x'range => DoblDobl_Complex_Numbers.Create(integer32(0)));
    nrm : double_float;

  begin
    z := DoblDobl_Complex_Monomials.Eval(m,x);
    put("z : "); put(z); new_line;
    DoblDobl_Complex_Monomials.Diff(m,x,zd);
    DoblDobl_Complex_Monomials.Speel(m,x,y,yd);
    put("y : "); put(y); new_line;
    DoblDobl_Complex_Numbers.Sub(y,z);
    nrm := hi_part(DoblDobl_Complex_Numbers.AbsVal(y));
    put("Norm of the difference :"); put(nrm,3); new_line;
    put_line("The derivatives computed with the Speelpenning algorithm :");
    put_line(yd(1..integer32(m.nvr)));
    put_line("The derivatives computed with the straightforward algorithm :");
    put_line(zd(1..integer32(m.nvr)));
    DoblDobl_Complex_Vectors.Sub(yd,zd);
    nrm := hi_part(DoblDobl_Complex_Vector_Norms.Max_Norm(yd));
    put("Max norm of the difference :"); put(nrm,3); new_line;
  end DoblDobl_Eval;

  procedure QuadDobl_Eval ( m : in QuadDobl_Complex_Monomials.Monomial ) is

  -- DESCRIPTION :
  --   Evaluates m at a random vector, in quad double precision.

    x : constant QuadDobl_Complex_Vectors.Vector(1..integer32(m.dim))
      := QuadDobl_Random_Vectors.Random_Vector(1,integer32(m.dim));
    y,z : QuadDobl_Complex_Numbers.Complex_Number;
    yd,zd : QuadDobl_Complex_Vectors.Vector(x'range)
          := (x'range => QuadDobl_Complex_Numbers.Create(integer32(0)));
    nrm : double_float;

  begin
    z := QuadDobl_Complex_Monomials.Eval(m,x);
    put("z : "); put(z); new_line;
    QuadDobl_Complex_Monomials.Diff(m,x,zd);
    QuadDobl_Complex_Monomials.Speel(m,x,y,yd);
    put("y : "); put(y); new_line;
    QuadDobl_Complex_Numbers.Sub(y,z);
    nrm := hihi_part(QuadDobl_Complex_Numbers.AbsVal(y));
    put("Norm of the difference :"); put(nrm,3); new_line;
    put_line("The derivatives computed with the Speelpenning algorithm :");
    put_line(yd(1..integer32(m.nvr)));
    put_line("The derivatives computed with the straightforward algorithm :");
    put_line(zd(1..integer32(m.nvr)));
    QuadDobl_Complex_Vectors.Sub(yd,zd);
    nrm := hihi_part(QuadDobl_Complex_Vector_Norms.Max_Norm(yd));
    put("Max norm of the difference :"); put(nrm,3); new_line;
  end QuadDobl_Eval;

  procedure Standard_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial with coefficient in standard double precision.

    dim : integer32 := 0;
    expmax : natural32 := 0;
    m : Standard_Complex_Monomials.Monomial;

  begin
    put_line("Testing monomial operations in double precision ...");
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    m := Random_Monomials.Standard_Random_Monomial(dim,expmax,true);
    put_line("A random monomial :"); put(m);
    Standard_Eval(m);
  end Standard_Test;

  procedure DoblDobl_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial with coefficient in double double precision.

    dim : integer32 := 0;
    expmax : natural32 := 0;
    m : DoblDobl_Complex_Monomials.Monomial;

  begin
    put_line("Testing monomial operations in double double precision ...");
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    m := Random_Monomials.DoblDobl_Random_Monomial(dim,expmax,true);
    put_line("A random monomial :"); put(m);
    DoblDobl_Eval(m);
  end DoblDobl_Test;

  procedure QuadDobl_Test is

  -- DESCRIPTION :
  --   Prompts for dimension, largest exponent, and then generates
  --   a random monomial with coefficient in quad double precision.

    dim : integer32 := 0;
    expmax : natural32 := 0;
    m : QuadDobl_Complex_Monomials.Monomial;

  begin
    put_line("Testing monomial operations in quad double precision ...");
    put("Give the dimension : "); get(dim);
    put("Give the largest exponent : "); get(expmax);
    m := Random_Monomials.QuadDobl_Random_Monomial(dim,expmax,true);
    put_line("A random monomial :"); put(m);
    QuadDobl_Eval(m);
  end QuadDobl_Test;

  procedure Main is

  -- DESCRIPTION :
  --   Prompts for the precision and then lauches the test.

    ans : character;

  begin
    new_line;
    put_line("MENU to select the precision :");
    put_line("  0. double precision");
    put_line("  1. double double precision");
    put_line("  2. quad double precision");
    put("Type 0, 1, or 2 to select the precision : ");
    Ask_Alternative(ans,"012");
    new_line;
    case ans is
      when '0' => Standard_Test;
      when '1' => DoblDobl_Test;
      when '2' => QuadDobl_Test;
      when others => null;
    end case;
  end Main;

begin
  Main;
end ts_monom;
